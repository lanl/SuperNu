* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module gasmod
c
      implicit none
c***********************************************************************
c gas grid structure
c***********************************************************************
      integer,parameter :: gas_nelem=111
      integer,parameter :: gas_nchain=3  !number of hardwired decay chains: ni56, fe52, cr48 
c
c-- available isotopes
      integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays
      integer,parameter :: gas_ife52=-3, gas_imn52=-4 !positions in mass0fr and natom1fr arrays
      integer,parameter :: gas_icr48=-5, gas_iv48=-6  !positions in mass0fr and natom1fr arrays
c
c-- input parameters
      real*8 :: gas_gastempinit
      real*8 :: gas_radtempinit
c
c-- wavelength grid (gridmod has a copy as well)
      integer,private :: ng=0
c
c-- domain decomposed grid variables used to calculate the state of the material (gas)
      integer :: gas_ncell=0
      integer :: gas_icell1=0
      real*8,allocatable :: gas_temp(:)       !(ncell)
      real*8,allocatable :: gas_eraddens(:)
      real*8,allocatable :: gas_ur(:)
      real*8,allocatable :: gas_rho(:)
      real*8,allocatable :: gas_bcoef(:)
      real*8,allocatable :: gas_decaygamma(:)
      real*8,allocatable :: gas_decaybeta(:)
      real*8,allocatable :: gas_vol(:)        !cell volume [cm^3]
      real*8,allocatable :: gas_mass(:)       !cell mass [g]
      real*8,allocatable :: gas_ye(:)         !electron fraction
      real*8,allocatable :: gas_ye0(:)        !initial electron fraction
      real*8,allocatable :: gas_dynfr(:)      !dynamical ejecta fraction
      real*8,allocatable :: gas_natom(:)      !cell number of atoms
      real*8,allocatable :: gas_nelec(:)      !cell number of electrons per atom
      real*8,allocatable :: gas_natom1fr(:,:) !(-2*gas_nchain:gas_nelem,ncell)  !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
      real*8,allocatable :: gas_natom0fr(:,:,:) !(-2:2,ncell,nchain) !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
c-- mate,allocatablerial energy (temperature) source (may be manufactured), rev>244
      real*8,allocatable :: gas_matsrc(:)     !-- material energy (temperature) source (may be manufactured)
c
      real*8,allocatable :: gas_edep(:)

c== DD copies
c-- Line+Cont extinction coeff
      real*8,allocatable :: gas_capcoef(:) !(ncell)
      real*4,allocatable :: gas_cap(:,:) !(ng,ncell)
      real*4,allocatable :: gas_em_cap(:,:) !(ng,ncell) -- emission opacity
c-- leakage opacities
c     real*8,allocatable :: dd_opacleak(:,:) !(6,ncell)
c-- scattering coefficient
      real*8,allocatable :: gas_sig(:) !(ncell)
c-- Gamma ray gray opacity
      real*8,allocatable :: gas_capgam(:) !(ncell)
c-- Planck opacity (gray)
      real*8,allocatable :: gas_capgrey(:)!(ncell)
      real*8,allocatable :: gas_em_capgrey(:)!(ncell) -- Planck-integrated emission opacity
c-- Fleck factor
      real*8,allocatable :: gas_fcoef(:)  !(ncell)
c
      real*8,allocatable :: gas_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
      real*8,allocatable :: gas_emitex(:) !(ncell) amount of external energy emitted per cell per group in a time step
      real*8,allocatable :: gas_evolinit(:) !(ncell) amount of initial energy per cell per group


c-- temperature structure history (only allocated when used)
      real*8,allocatable :: gas_temppreset(:,:) !(ncell,tim_nt)
c
      save
c
      contains
c
c
      subroutine gasmod_init(ltalk,icell1,ncell,ngin,lemiss)
c----------------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk
      integer,intent(in) :: icell1,ncell,ngin
      logical,intent(in) :: lemiss
************************************************************************
* Allocate gas variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      ng = ngin
c
      gas_icell1 = icell1
      gas_ncell = ncell
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = gas_ncell*(8*(21 + (gas_nelem + 1 + 2*gas_nchain) +
     &   (gas_nchain*5)) + 4*(2 + ng))/1024 !kB
       write(6,*) 'ALLOC gas      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif !ltalk
c
c-- ndim=1 alloc
      allocate(gas_temp(gas_ncell))
      allocate(gas_ur(gas_ncell))
      allocate(gas_rho(gas_ncell))
      allocate(gas_bcoef(gas_ncell))
      allocate(gas_decaygamma(gas_ncell))
      gas_decaygamma = 0d0
      allocate(gas_decaybeta(gas_ncell))
      gas_decaybeta = 0d0
      allocate(gas_vol(gas_ncell))
      allocate(gas_mass(gas_ncell))
      allocate(gas_ye(gas_ncell))
      allocate(gas_ye0(gas_ncell))
      allocate(gas_dynfr(gas_ncell))
      allocate(gas_natom(gas_ncell))
      gas_natom = 0d0
      allocate(gas_nelec(gas_ncell))
      gas_nelec = 1d0
      allocate(gas_matsrc(gas_ncell))
      gas_matsrc = 0d0
c   allocate(dd_opacleak(6,gas_ncell))
      allocate(gas_sig(gas_ncell))
      allocate(gas_capgam(gas_ncell))
      allocate(gas_capgrey(gas_ncell))
      allocate(gas_capcoef(gas_ncell))
      allocate(gas_fcoef(gas_ncell))
c
      allocate(gas_eraddens(gas_ncell))
      allocate(gas_edep(gas_ncell))
c
      allocate(gas_emit(gas_ncell))
      allocate(gas_emitex(gas_ncell))
      allocate(gas_evolinit(gas_ncell))
c
c-- ndim=2 alloc small
      allocate(gas_natom1fr(-2*gas_nchain:gas_nelem,gas_ncell))
      gas_natom1fr = 0d0
      allocate(gas_natom0fr(-2:2,gas_ncell,gas_nchain))
      gas_natom0fr = 0d0
c
c-- ndim=2 alloc big
      allocate(gas_cap(ng,gas_ncell))
c-- optional emissivity arrays
      if (lemiss) then
        allocate(gas_em_cap(ng,gas_ncell))
        allocate(gas_em_capgrey(gas_ncell))
      endif
!}}}
      end subroutine gasmod_init
c
c
      subroutine gas_dealloc
      deallocate(gas_temp)!{{{
      deallocate(gas_ur)
      deallocate(gas_rho)
      deallocate(gas_bcoef)
      deallocate(gas_decaygamma)
      deallocate(gas_decaybeta)
      deallocate(gas_vol)
      deallocate(gas_mass)
      deallocate(gas_ye)
      deallocate(gas_ye0)
      deallocate(gas_dynfr)
      deallocate(gas_natom)
      deallocate(gas_nelec)
      deallocate(gas_matsrc)
      deallocate(gas_sig)
      deallocate(gas_capgam)
      deallocate(gas_capgrey)
      deallocate(gas_capcoef)
      deallocate(gas_fcoef)
      deallocate(gas_eraddens)
      deallocate(gas_edep)
      deallocate(gas_emit)
      deallocate(gas_emitex)
      deallocate(gas_evolinit)
      deallocate(gas_natom1fr)
      deallocate(gas_natom0fr)
      deallocate(gas_cap)
      if (allocated(gas_em_cap)) then
        deallocate(gas_em_cap)
        if (.not.allocated(gas_em_capgrey)) then
          stop 'gas_em_cap alloc-ed but not gas_em_capgrey'
        endif
        deallocate(gas_em_capgrey)
      endif
!}}}
      end subroutine gas_dealloc
c
      end module gasmod
c vim: fdm=marker
