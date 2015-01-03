      module gridmod
c     --------------
      implicit none
c
      logical :: grd_isvelocity = .false.
c
      integer :: grd_igeom = 0
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
      real*8,allocatable :: grd_xarr(:)   !(nx+1), left cell edge values
      real*8,allocatable :: grd_yarr(:)   !(ny+1), left cell edge values
      real*8,allocatable :: grd_zarr(:)   !(nz+1), left cell edge values
c-- polar angles
      real*8,allocatable :: grd_yacos(:)   !(ny+1)
c-- pointer into compressed domain
      integer,allocatable :: grd_icell(:,:,:) !(nx,ny,nz)
c
c
c-- compressed domain
      integer :: grd_nc=0  !number of cells
      integer :: grd_ncp=0 !with padding
c
c-- Probability of emission in a given zone and group
      real*8,allocatable :: grd_emitprob(:,:) !(nep,ncp)
c-- Line+Cont extinction coeff
      real*4,allocatable :: grd_cap(:,:) !(ng,ncp)
c-- leakage opacities
      real*8,allocatable :: grd_opacleak(:,:) !(6,ncp)

c
      real*8,allocatable :: grd_vol(:)  !(ncp)
      real*8,allocatable :: grd_temp(:) !(ncp)

c-- scattering coefficient
      real*8,allocatable :: grd_sig(:) !(ncp) !grey scattering opacity
c-- Planck opacity (gray)
      real*8,allocatable :: grd_capgrey(:) !(ncp)
c-- Fleck factor
      real*8,allocatable :: grd_fcoef(:)  !(ncp)


c-- energy absorbed by material
      real*8,allocatable :: grd_edep(:)   !(ncp)
      real*8,allocatable :: grd_eamp(:)   !(ncp)
c-- radiation energy density in tsp_dt
      real*8,allocatable :: grd_eraddens(:) !(ncp)


c-- number of IMC-DDMC method changes per cell per time step
      integer,allocatable :: grd_methodswap(:) !(ncp)
c-- number of census prt_particles per cell
      integer,allocatable :: grd_numcensus(:) !(ncp)

c
c-- packet number and energy distribution
c========================================
      integer,allocatable :: grd_nvol(:) !(ncp) number of thermal source particles generated per cell
      integer,allocatable :: grd_nvolinit(:) !(ncp) number of initial (t=tfirst) particles per cell
c      
      real*8,allocatable :: grd_emit(:) !(ncp) amount of fictitious thermal energy emitted per cell in a time step
      real*8,allocatable :: grd_emitex(:) !(ncp) amount of external energy emitted per cell in a time step
      real*8,allocatable :: grd_evolinit(:) !(ncp) amount of initial energy per cell per group
c
c-- temperature structure history (allocated only if used)
      real*8,allocatable :: grd_temppreset(:,:) !(ncp,tim_nt)
c
      save
c
      contains
c
      subroutine grid_init(ltalk,ngin,igeom,ndim,nc,ncp,isvelocity)
c     -------------------------------------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk,isvelocity
      integer,intent(in) :: ngin,igeom
      integer,intent(in) :: ndim(3)
      integer,intent(in) :: nc,ncp
************************************************************************
* Allocate grd variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n
c
      grd_igeom = igeom
c
      ng = ngin
c
c-- emission probability
      grd_nep = nint(sqrt(dble(ng)))
      grd_nepg = ceiling(ng/(grd_nep + 1d0))
c
      grd_nx = ndim(1)
      grd_ny = ndim(2)
      grd_nz = ndim(3)
c
c-- ncell, ncell-with-padding
      grd_nc = nc
      grd_ncp = ncp
c
      grd_isvelocity = isvelocity
c
      allocate(grd_xarr(grd_nx+1))
      allocate(grd_yarr(grd_ny+1))
      allocate(grd_zarr(grd_nz+1))
c-- polar
      if(igeom==1) allocate(grd_yacos(grd_ny+1))
c
c-- complete domain
      allocate(grd_icell(grd_nx,grd_ny,grd_nz))
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = int((int(grd_ncp,8)*(8*(11+6) + 4*4))/1024) !kB
       write(6,*) 'ALLOC grd      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
       n = int((int(grd_ncp,8)*4*ng)/1024) !kB
       write(6,*) 'ALLOC grd_cap  :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif
c
c-- ndim=3 alloc
      allocate(grd_edep(grd_ncp))
      allocate(grd_eamp(grd_ncp))
      allocate(grd_capgrey(grd_ncp))
      allocate(grd_sig(grd_ncp))
      allocate(grd_fcoef(grd_ncp))
      allocate(grd_eraddens(grd_ncp))
      allocate(grd_temp(grd_ncp))
      allocate(grd_vol(grd_ncp))
c
      allocate(grd_emit(grd_ncp))
      allocate(grd_emitex(grd_ncp))
      allocate(grd_evolinit(grd_ncp))
c
c-- ndim=3 integer
      allocate(grd_nvol(grd_ncp))
      allocate(grd_nvolinit(grd_ncp))
c
      allocate(grd_methodswap(grd_ncp))
      allocate(grd_numcensus(grd_ncp))
c
c-- ndim=4 alloc
      allocate(grd_opacleak(6,grd_ncp))
      allocate(grd_emitprob(grd_nep,grd_ncp))
c-- ndim=4 alloc
      allocate(grd_cap(ng,grd_ncp))
c!}}}
      end subroutine grid_init
c
      end module gridmod
