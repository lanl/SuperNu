*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module gridmod
c     --------------
      implicit none
c
      logical :: grd_isvelocity=.false.
c
      integer :: grd_igeom=0
c
      integer,private :: ng=0
c
      integer :: grd_nep=0 !number of emission probability bins
      integer :: grd_nepg=0 !number of groups per emission probability bin
c
c-- complete domain
      integer :: grd_nx=0
      integer :: grd_ny=0
      integer :: grd_nz=0
c
      real*8,allocatable :: grd_xarr(:)  !(nx+1), left cell edge values
      real*8,allocatable :: grd_yarr(:)  !(ny+1), left cell edge values
      real*8,allocatable :: grd_zarr(:)  !(nz+1), left cell edge values
c
c-- maximum radial grid velocity
      real*8 :: grd_rout=0d0   !particle flux edge radius
      real*8 :: grd_rvoid=0d0  !void-corner radius, used as cell criterium
c
c-- polar angles
      real*8,allocatable :: grd_yacos(:)   !(ny+1)
c-- pointer into compressed domain
      integer,allocatable :: grd_icell(:,:,:) !(nx,ny,nz)
c
c-- domain decomposition
      integer :: grd_idd1=0
      integer :: grd_ndd=0
c
c-- compressed domain
      integer :: grd_ncell=0  !number of cells
      integer :: grd_ivoid=0  !the void cell id
c
c-- Probability of emission in a given zone and group
      real*8,allocatable :: grd_emitprob(:,:) !(nep,ncell)

c-- Line+Cont extinction coeff
      real*4,allocatable :: grd_cap(:,:) !(ng,ncell)

c-- leakage opacities
      real*8,allocatable :: grd_opaclump(:,:) !(10,ncell) leak(6),speclump,caplump,igemitmax,doplump
      real*8,allocatable :: grd_tempinv(:) !(ncell)
c-- scattering coefficient
      real*8,allocatable :: grd_sig(:) !(ncell) !grey scattering opacity
c-- Planck opacity (gray)
      real*8,allocatable :: grd_capgrey(:) !(ncell)
c-- Fleck factor
      real*8,allocatable :: grd_fcoef(:)  !(ncell)



      real*8,allocatable :: grd_tally(:,:)   !(2,ncell) (edep,eraddens)
c-- amplification factor excess
      real*8,allocatable :: grd_eamp(:)   !(ncell)


c-- number of IMC-DDMC method changes per cell per time step
      integer,allocatable :: grd_methodswap(:) !(ncell)
c-- number of census prt_particles per cell
      integer,allocatable :: grd_numcensimc(:) !(ncell)
      integer,allocatable :: grd_numcensddmc(:) !(ncell)

c
c-- packet number and energy distribution
c========================================
      real*8,allocatable :: grd_capgam(:)   !(ncell) gray gamma opacity
c
      real*8,allocatable :: grd_vol(:)  !(ncell)
c
      integer,allocatable :: grd_nvol(:) !(ncell) number of thermal source particles generated per cell
      integer,allocatable :: grd_nvolinit(:) !(ncell) number of initial (t=tfirst) particles per cell
c
      real*8,allocatable :: grd_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
      real*8,allocatable :: grd_emitex(:) !(ncell) amount of external energy emitted per cell in a time step
      real*8,allocatable :: grd_evolinit(:) !(ncell) amount of initial energy per cell per group
c
c-- temperature structure history (allocated only if used)
      real*8,allocatable :: grd_temppreset(:,:) !(ncell,tim_nt)
c
      interface
      pure function emitgroup(r,ic) result(ig)
      integer :: ig
      real*8,intent(in) :: r
      integer,intent(in) :: ic
      end function emitgroup
      end interface
c
      save
c
      contains
c
      subroutine gridmod_init(ltalk,ngin,ncell,lvoid,idd1,ndd)
c     --------------------------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk,lvoid
      integer,intent(in) :: ngin
      integer,intent(in) :: ncell,idd1,ndd
************************************************************************
* Allocate grd variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      ng = ngin
c-- void cell exists
      if(lvoid) grd_ivoid = ncell
c
c-- emission probability
      grd_nep = nint(sqrt(dble(ng)))
      grd_nepg = ceiling(ng/(grd_nep + 1d0))
c
c-- number of non-void cells, plus one optional dummy cell if void cells exist
      grd_ncell = ncell
      grd_idd1 = idd1
      grd_ndd = ndd
c
      allocate(grd_xarr(grd_nx+1))
      allocate(grd_yarr(grd_ny+1))
      allocate(grd_zarr(grd_nz+1))
c-- polar
      if(grd_igeom==1) allocate(grd_yacos(grd_ny+1))
c
c-- complete domain
      allocate(grd_icell(grd_nx,grd_ny,grd_nz))
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = int((int(grd_ncell,8)*(8*(12+6) + 5*4))/1024) !kB
       write(6,*) 'ALLOC grd      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
       n = int((int(grd_ncell,8)*4*ng)/1024) !kB
       write(6,*) 'ALLOC grd_cap  :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif
c
c-- ndim=3 alloc
      allocate(grd_tally(2,grd_ncell))
      allocate(grd_eamp(grd_ncell))
      allocate(grd_capgrey(grd_ncell))
      allocate(grd_sig(grd_ncell))
      allocate(grd_fcoef(grd_ncell))
      allocate(grd_tempinv(grd_ncell))
      allocate(grd_vol(grd_ncell))
c
      allocate(grd_capgam(grd_ncell))
      allocate(grd_emit(grd_ncell))
      grd_emit = 0d0
      allocate(grd_emitex(grd_ncell))
      allocate(grd_evolinit(grd_ncell))
c
c-- ndim=3 integer
      allocate(grd_nvol(grd_ncell))
      grd_nvol = 0
      allocate(grd_nvolinit(grd_ncell))
      grd_nvolinit = 0
c
      allocate(grd_methodswap(grd_ncell))
      allocate(grd_numcensimc(grd_ncell))
      allocate(grd_numcensddmc(grd_ncell))
c
c-- ndim=4 alloc
      allocate(grd_opaclump(10,grd_ncell))
      allocate(grd_emitprob(grd_nep,grd_ncell))
c-- ndim=4 alloc
      allocate(grd_cap(ng,grd_ncell))
c!}}}
      end subroutine gridmod_init
c
c
      subroutine grid_dealloc
      deallocate(grd_xarr)!{{{
      deallocate(grd_yarr)
      deallocate(grd_zarr)
c-- polar
      if(grd_igeom==1) deallocate(grd_yacos)
c-- complete domain
      deallocate(grd_icell)
c-- gasmod
      deallocate(grd_tally)
      deallocate(grd_eamp)
      deallocate(grd_capgrey)
      deallocate(grd_sig)
      deallocate(grd_fcoef)
      deallocate(grd_tempinv)
      deallocate(grd_vol)
c
      deallocate(grd_capgam)
      deallocate(grd_emit)
      deallocate(grd_emitex)
      deallocate(grd_evolinit)
c-- ndim=3 integer
      deallocate(grd_nvol)
      deallocate(grd_nvolinit)
      deallocate(grd_methodswap)
      deallocate(grd_numcensimc)
      deallocate(grd_numcensddmc)
c-- ndim=4 alloc
      deallocate(grd_opaclump)
      deallocate(grd_emitprob)
c-- ndim=4 alloc
      deallocate(grd_cap)!}}}
      end subroutine grid_dealloc
c
      end module gridmod
c vim: fdm=marker
