      module gridmod
c     --------------
      implicit none
c
      integer :: grd_igeom = 0
c
      integer :: grd_ng=0
      real*8,allocatable :: grd_wl(:)
c
      integer :: grd_nx=0
      integer :: grd_ny=0
      integer :: grd_nz=0

      logical :: grd_isvelocity = .false.

      real*8,allocatable :: grd_xarr(:)   !(nx+1), left cell edge values
      real*8,allocatable :: grd_yarr(:)   !(ny+1), left cell edge values
      real*8,allocatable :: grd_zarr(:)   !(nz+1), left cell edge values


!-- Probability of emission in a given zone and group
      real*8,allocatable :: grd_emitprob(:,:,:,:) !(ng,nx,ny,nz)
!-- Line+Cont extinction coeff
      real*8,allocatable :: grd_cap(:,:,:,:) !(ng,nx,ny,nz)
!-- leakage opacities
      real*8,allocatable :: grd_opacleak(:,:,:,:) !(6,nx,ny,nz)


!-- scattering coefficient
      real*8,allocatable :: grd_sig(:,:,:) !(nx,ny,nz)
!-- Gamma ray gray opacity
      real*8,allocatable :: grd_capgam(:,:,:) !(nx,ny,nz)
!-- Planck opacity (gray)
      real*8,allocatable :: grd_siggrey(:,:,:)!(nx,ny,nz)
!-- Fleck factor
      real*8,allocatable :: grd_fcoef(:,:,:)  !(nx,ny,nz)


!-- energy absorbed by material
      real*8,allocatable :: grd_edep(:,:,:)   !(nx,ny,nz)
!-- radiation energy density in tsp_dt
      real*8,allocatable :: grd_eraddens(:,:,:) !(nx,ny,nz)


!-- number of IMC-DDMC method changes per cell per time step
      integer,allocatable :: grd_methodswap(:,:,:) !(nx,ny,nz)
!-- number of census prt_particles per cell
      integer,allocatable :: grd_numcensus(:,:,:) !(nx,ny,nz)

!
!-- packet number and energy distribution
!========================================
      integer,allocatable :: grd_nvol(:,:,:) !(nx,ny,nz) number of thermal source particles generated per cell
      integer,allocatable :: grd_nvolex(:,:,:) !(nx,ny,nz) number of external source particles generated per cell
      integer,allocatable :: grd_nvolinit(:,:,:) !(nx,ny,nz) number of initial (t=tfirst) particles per cell
!      
      real*8,allocatable :: grd_emit(:,:,:) !(nx,ny,nz) amount of fictitious thermal energy emitted per cell in a time step
      real*8,allocatable :: grd_emitex(:,:,:) !(nx,ny,nz) amount of external energy emitted per cell per group in a time step
      real*8,allocatable :: grd_evolinit(:,:,:) !(nx,ny,nz) amount of initial energy per cell per group
!
      real*8,allocatable :: grd_temp(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: grd_vol(:,:,:) !(nx,ny,nz)
c
c-- temperature structure history (allocated only if used)
      real*8,allocatable :: grd_temppreset(:,:,:,:) !(nx,ny,nz,tim_nt)
c
      save
c
      contains
c
      subroutine grid_init(ltalk,ng,wlarr,igeom,ndim,isvelocity)
c     --------------------------------!{{{
      implicit none
      logical,intent(in) :: ltalk,isvelocity
      integer,intent(in) :: ng,igeom
      integer,intent(in) :: ndim(3)
      real*8,intent(in) :: wlarr(ng+1)
c
      integer :: n,nx,ny,nz
c
      grd_igeom = igeom
c
      grd_ng = ng
      allocate(grd_wl(ng+1))
      grd_wl = wlarr
c
      grd_nx = ndim(1)
      grd_ny = ndim(2)
      grd_nz = ndim(3)

      grd_isvelocity = isvelocity

      allocate(grd_xarr(grd_nx+1))
      allocate(grd_yarr(grd_ny+1))
      allocate(grd_zarr(grd_nz+1))

      nx = grd_nx !shortcut
      ny = grd_ny !shortcut
      nz = grd_nz !shortcut
      allocate(grd_edep(nx,ny,nz))
      allocate(grd_siggrey(nx,ny,nz))
      allocate(grd_capgam(nx,ny,nz))
!
!- Ryan W.: using power law to calculate grd_sig (similar to Planck opacity)
      allocate(grd_sig(nx,ny,nz))    !grey scattering opacity
!----------------------------------------------------------------
      allocate(grd_fcoef(nx,ny,nz))
       allocate(grd_emitprob(ng,nx,ny,nz))
      allocate(grd_opacleak(6,nx,ny,nz))
      allocate(grd_eraddens(nx,ny,nz))

!-Ryan W: gas_wl being allocated in gasgrid_setup now--
      !allocate(gas_wl(ng)) !wavelength grid
!------------------------------------------------------
       allocate(grd_cap(ng,nx,ny,nz)) !Line+Cont extinction coeff

!--Ryan W: values below were formerly secondary (rev 183)
      allocate(grd_temp(nx,ny,nz))  !cell average temperature
      allocate(grd_vol(nx,ny,nz))  !cell average temperature

      allocate(grd_nvol(nx,ny,nz))
      allocate(grd_nvolex(nx,ny,nz))
      allocate(grd_nvolinit(nx,ny,nz))
      allocate(grd_emit(nx,ny,nz))
      allocate(grd_emitex(nx,ny,nz))
      allocate(grd_evolinit(nx,ny,nz))

      allocate(grd_methodswap(nx,ny,nz))
      allocate(grd_numcensus(nx,ny,nz))

      if(ltalk) then
       n = nx*ny*nz
       n = int((int(n,8)*8*21)/1024) !kB
       write(6,*) 'ALLOC grid:',n,"kB",n/1024,"MB",n/1024**2,"GB"
       n = nx*ny*nz
       n = int(((8+8)*int(n,8)*ng)/1024) !kB
       write(6,*) 'ALLOC grd_cap:',n,"kB",n/1024,"MB",n/1024**2,"GB"
      endif
c!}}}
      end subroutine grid_init
c
      end module gridmod
