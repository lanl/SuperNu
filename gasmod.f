      module gasmod
c
      implicit none
c***********************************************************************
c gas grid structure
c***********************************************************************
      integer,parameter :: gas_nelem=30
      integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays
c
c-- wavelength grid (gridmod has a copy as well)
      integer :: gas_ng=0
      real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
c
c-- domain decomposed grid variables used to calculate the state of the material (gas)
      integer :: gas_ncell=0
      real*8,allocatable :: gas_temp(:)       !(ncell)
      real*8,allocatable :: gas_eraddens(:)
      real*8,allocatable :: gas_ur(:)
      real*8,allocatable :: gas_rho(:)
      real*8,allocatable :: gas_bcoef(:)
      real*8,allocatable :: gas_nisource(:)
      real*8,allocatable :: gas_vol(:)        !cell volume [cm^3]
      real*8,allocatable :: gas_mass(:)       !cell mass [g]
      real*8,allocatable :: gas_ye(:)         !electron fraction
      real*8,allocatable :: gas_natom(:)      !cell number of atoms
      real*8,allocatable :: gas_nelec(:)      !cell number of electrons per atom
      real*8,allocatable :: gas_natom1fr(:,:) !(-2:gas_nelem,ncell)  !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
      real*8,allocatable :: gas_natom0fr(:,:) !(-2:2,ncell)  !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
c-- mate,allocatablerial energy (temperature) source (may be manufactured), rev>244
      real*8,allocatable :: gas_matsrc(:)     !-- material energy (temperature) source (may be manufactured)
c
      real*8,allocatable :: gas_edep(:)

c== DD copies
c-- Probability of emission in a given zone and group
      real*8,allocatable :: gas_emitprob(:,:) !(gas_ng,ncell)
c-- Line+Cont extinction coeff
      real*8,allocatable :: gas_cap(:,:) !(gas_ng,ncell)
c-- leakage opacities
c     real*8,allocatable :: dd_opacleak(:,:) !(6,ncell)
c-- scattering coefficient
      real*8,allocatable :: gas_sig(:) !(ncell)
c-- Gamma ray gray opacity
      real*8,allocatable :: gas_capgam(:) !(ncell)
c-- Planck opacity (gray)
      real*8,allocatable :: gas_capgrey(:)!(ncell)
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
      subroutine gas_init(ltalk,ncell)
c-------------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk
      integer,intent(in) :: ncell
************************************************************************
* Allocate gas variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      gas_ncell = ncell
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = gas_ncell*(20 + (gas_nelem+3+5) + 2*gas_ng)/1024 !kB
       write(6,*) 'ALLOC gas    :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif !ltalk
c
c-- ndim=1 alloc
      allocate(gas_temp(gas_ncell))
      allocate(gas_ur(gas_ncell))
      allocate(gas_rho(gas_ncell))
      allocate(gas_bcoef(gas_ncell))
      allocate(gas_nisource(gas_ncell))
      allocate(gas_vol(gas_ncell))
      allocate(gas_mass(gas_ncell))
      allocate(gas_ye(gas_ncell))
      allocate(gas_natom(gas_ncell))
      allocate(gas_nelec(gas_ncell))
      allocate(gas_matsrc(gas_ncell))
      gas_ye = .5d0
      gas_nelec = 1d0
      gas_matsrc = 0d0
c   allocate(dd_opacleak(6,gas_ncell))
      allocate(gas_sig(gas_ncell))
      allocate(gas_capgam(gas_ncell))
      allocate(gas_capgrey(gas_ncell))
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
      allocate(gas_natom1fr(-2:gas_nelem,gas_ncell))
      allocate(gas_natom0fr(-2:2,gas_ncell))
      gas_natom1fr = 0d0
      gas_natom0fr = 0d0
c
c-- ndim=2 alloc big
      allocate(gas_emitprob(gas_ng,gas_ncell))
      allocate(gas_cap(gas_ng,gas_ncell))
!}}}
      end subroutine gas_init
c
      end module gasmod
