      module gridmod
c     --------------
      implicit none
c
      integer :: grd_igeom = 0
c
      integer,private :: ng=0
c
      integer :: grd_nep=0 !number of emission probability bins
      integer :: grd_nepg=0 !number of groups per emission probability bin
c
      integer :: grd_nx=0
      integer :: grd_ny=0
      integer :: grd_nz=0

      logical :: grd_isvelocity = .false.

      real*8,allocatable :: grd_xarr(:)   !(nx+1), left cell edge values
      real*8,allocatable :: grd_yarr(:)   !(ny+1), left cell edge values
      real*8,allocatable :: grd_zarr(:)   !(nz+1), left cell edge values


c-- Probability of emission in a given zone and group
      real*8,allocatable :: grd_emitprob(:,:,:,:) !(nep,nx,ny,nz)
c-- Line+Cont extinction coeff
      real*4,allocatable :: grd_cap(:,:,:,:) !(ng,nx,ny,nz)
c-- leakage opacities
      real*8,allocatable :: grd_opacleak(:,:,:,:) !(6,nx,ny,nz)

c
      real*8,allocatable :: grd_temp(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: grd_vol(:,:,:) !(nx,ny,nz)

c-- scattering coefficient
      real*8,allocatable :: grd_sig(:,:,:) !(nx,ny,nz) !grey scattering opacity
c-- Gamma ray gray opacity
      real*8,allocatable :: grd_capgam(:,:,:) !(nx,ny,nz)
c-- Planck opacity (gray)
      real*8,allocatable :: grd_capgrey(:,:,:)!(nx,ny,nz)
c-- Fleck factor
      real*8,allocatable :: grd_fcoef(:,:,:)  !(nx,ny,nz)


c-- energy absorbed by material
      real*8,allocatable :: grd_edep(:,:,:)   !(nx,ny,nz)
      real*8,allocatable :: grd_eamp(:,:,:)   !(nx,ny,nz)
c-- radiation energy density in tsp_dt
      real*8,allocatable :: grd_eraddens(:,:,:) !(nx,ny,nz)


c-- number of IMC-DDMC method changes per cell per time step
      integer,allocatable :: grd_methodswap(:,:,:) !(nx,ny,nz)
c-- number of census prt_particles per cell
      integer,allocatable :: grd_numcensus(:,:,:) !(nx,ny,nz)

c
c-- packet number and energy distribution
c========================================
      integer,allocatable :: grd_nvol(:,:,:) !(nx,ny,nz) number of thermal source particles generated per cell
      integer,allocatable :: grd_nvolinit(:,:,:) !(nx,ny,nz) number of initial (t=tfirst) particles per cell
c      
      real*8,allocatable :: grd_emit(:,:,:) !(nx,ny,nz) amount of fictitious thermal energy emitted per cell in a time step
      real*8,allocatable :: grd_emitex(:,:,:) !(nx,ny,nz) amount of external energy emitted per cell in a time step
      real*8,allocatable :: grd_evolinit(:,:,:) !(nx,ny,nz) amount of initial energy per cell per group
c
c-- temperature structure history (allocated only if used)
      real*8,allocatable :: grd_temppreset(:,:,:,:) !(nx,ny,nz,tim_nt)
c
      save
c
      contains
c
      subroutine grid_init(ltalk,ngin,igeom,ndim,isvelocity)
c     --------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk,isvelocity
      integer,intent(in) :: ngin,igeom
      integer,intent(in) :: ndim(3)
************************************************************************
* Allocate grd variables.
*
* Don't forget to update the print statement if variables are added or
* removed
************************************************************************
      integer :: n,nx,ny,nz
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

      grd_isvelocity = isvelocity

      allocate(grd_xarr(grd_nx+1))
      allocate(grd_yarr(grd_ny+1))
      allocate(grd_zarr(grd_nz+1))
c
c-- shortcuts
      nx = grd_nx
      ny = grd_ny
      nz = grd_nz
c
c-- print alloc size (keep this updated)
c---------------------------------------
      if(ltalk) then
       n = nx*ny*nz
       n = int((int(n,8)*(8*(12+6) + 4*4))/1024) !kB
       write(6,*) 'ALLOC grd      :',n,"kB",n/1024,"MB",n/1024**2,"GB"
       n = nx*ny*nz
       n = int((int(n,8)*4*ng)/1024) !kB
       write(6,*) 'ALLOC grd_cap  :',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif
c
c-- ndim=3 alloc
      allocate(grd_edep(nx,ny,nz))
      allocate(grd_eamp(nx,ny,nz))
      allocate(grd_capgrey(nx,ny,nz))
      allocate(grd_capgam(nx,ny,nz))
      allocate(grd_sig(nx,ny,nz))
      allocate(grd_fcoef(nx,ny,nz))
      allocate(grd_eraddens(nx,ny,nz))
      allocate(grd_temp(nx,ny,nz))
      allocate(grd_vol(nx,ny,nz))
c
      allocate(grd_emit(nx,ny,nz))
      allocate(grd_emitex(nx,ny,nz))
      allocate(grd_evolinit(nx,ny,nz))
c
c-- ndim=3 integer
      allocate(grd_nvol(nx,ny,nz))
      allocate(grd_nvolinit(nx,ny,nz))
c
      allocate(grd_methodswap(nx,ny,nz))
      allocate(grd_numcensus(nx,ny,nz))
c
c-- ndim=4 alloc
      allocate(grd_opacleak(6,nx,ny,nz))
      allocate(grd_emitprob(grd_nep,nx,ny,nz))
c-- ndim=4 alloc
      allocate(grd_cap(ng,nx,ny,nz))
c!}}}
      end subroutine grid_init
c
      end module gridmod
