      module fluxmod
c     ------------------
      implicit none

      integer :: flx_ng=0  !number of flux groups (1D,2D,3D)
      integer :: flx_nmu=0 !number of polar bins (2D,3D)
      integer :: flx_nom=0 !number of azimuthal bins (3D)
c
c-- wavelength (wl), polar (mu), and azimuthal (om) bins
      real*8,allocatable :: flx_wl(:) !(flx_ng)
      real*8,allocatable :: flx_mu(:) !(flx_nmu)
      real*8,allocatable :: flx_om(:) !(flx_nom)
!
!-- radiative flux
!=================
!-- outbound grouped luminosity
      real*8,allocatable :: flx_luminos(:,:,:) !(flx_ng,flx_nmu,flx_nom)
!-- sampled devation of group luminosity
      real*8,allocatable :: flx_lumdev(:,:,:) !(flx_ng,flx_nmu,flx_nom)
!-- number of escaped particles per group
      integer,allocatable :: flx_lumnum(:,:,:) !(flx_ng,flx_nmu,flx_nom)
c
c-- grey gamma flux
      real*8,allocatable :: flx_gamluminos(:,:) !(flx_nmu,flx_nom)
      real*8,allocatable :: flx_gamlumdev(:,:)  !(flx_nmu,flx_nom)
      integer,allocatable :: flx_gamlumnum(:,:) !(flx_nmu,flx_nom)
c
      save
c
      contains
c
c
      subroutine fluxgrid_setup(nflx,wlmin,wlmax)
c     --------------------------------------------
      implicit none
      integer,intent(in) :: nflx(3)
      real*8,intent(in) :: wlmin,wlmax
*******************************************************
* calculate flux grid and tally arrays:
* -- if any dimensions of flux grid are negative, a
* file could be read for the corresponding grid using
* the absolute value of the negative dimension as the
* index.
* -- if a dimension number is positive, it could be
* treated as the grid size.
* -- for wavelength: if nflx(1)=0, then gas_wl could
* be used as the tally grid, flx_wl.
*******************************************************
      integer :: i
      character(12) :: fnames(3)
c
      fnames = (/'input.fluxwl','input.fluxmu',
     &     'input.fluxom'/)
c
c-- check if bins are read or generated
      do i = 1,3
         if(nflx(i)<0) then
            call read_fluxgrid(i,nflx(i),fnames(i))
         else
            call generate_fluxgrid(i,nflx(i),wlmin,wlmax)
         endif
      enddo
c
c-- allocate flux tally arrays
      allocate(flx_luminos(flx_ng,flx_nmu,flx_nom))
      allocate(flx_lumdev(flx_ng,flx_nmu,flx_nom))
      allocate(flx_lumnum(flx_ng,flx_nmu,flx_nom))
c
c-- grey gamma flux
      allocate(flx_gamluminos(flx_nmu,flx_nom))
      allocate(flx_gamlumdev(flx_nmu,flx_nom))
      allocate(flx_gamlumnum(flx_nmu,flx_nom))
c
      end subroutine fluxgrid_setup
c
c
      subroutine read_fluxgrid(iflx,idex,fname)
c     ------------------------------
      use physconstmod
      implicit none
      integer,intent(in) :: idex,iflx
      character(12),intent(in) :: fname
*************************************************************
* read lab wavelength grid for flux from file
*************************************************************
      logical :: lexists
      integer :: ihelp,l,n,irr
      real*8 :: help
      real*8,allocatable :: store(:)
c
c-- sanity check
      if(iflx>3.or.iflx<1) stop 'read_fluxgrid: invalid 1st arg'
      
      inquire(file=fname,exist=lexists)
      if(.not.lexists) stop 'read_fluxgrid: missing file'
c
c-- grid lookup (rtw: this routine is merely wlgrid_setup)
      ihelp = abs(idex)
      open(4,file=fname,status='old')
c
c-- strip header
      do l=1,3
       read(4,*)
      enddo
c-- read dimensions
      read(4,*) !not used
c
      do l=1,ihelp-1
       read(4,*)
      enddo
      read(4,*) irr, n
c
      allocate(store(n+3))
      rewind(4)
      do l=1,ihelp+3
       read(4,*)
      enddo
      read(4,*) store(:)
      close(4)

c
c-- check monotonicity
      help = -2d0
      do l=1,n+1
       if(store(l+2)<=help) stop 'read_fluxgrid: grid decreasing'
       help = store(l+2)
      enddo

      if(iflx==1) then
         flx_ng = n
         allocate(flx_wl(n+1))
         flx_wl = store(3:)
      elseif(iflx==2) then
         flx_nmu = n
         allocate(flx_mu(n+1))
         flx_mu = store(3:)
      else
         flx_nom = n
         allocate(flx_om(n+1))
         flx_om = pc_pi2*store(3:)
      endif
      deallocate(store)
c
      end subroutine read_fluxgrid
c
c
      subroutine generate_fluxgrid(iflx,idex,wlmin,wlmax)
c     ---------------------------------------
      use physconstmod
      use groupmod
      implicit none
      integer,intent(in) :: idex,iflx
      real*8,intent(in) :: wlmin,wlmax
*************************************************************
* generate lab wavelength grid for flux
*************************************************************
      integer :: i
      real*8 :: help
c
c-- wavelength
      if(iflx==1) then
         if(idex==0) then
c-- set wl to transport grid
            flx_ng = grp_ng
            allocate(flx_wl(flx_ng+1))
            flx_wl = grp_wl
c
         elseif(idex>0) then
c-- logarithmic wavelength
            flx_ng = idex
            allocate(flx_wl(flx_ng+1))
            help = wlmax/dble(wlmin)
            forall(i=1:flx_ng+1) flx_wl(i) =
     &        wlmin*help**((i-1d0)/flx_ng)
         else
            stop 'generate_fluxgrid: invalid nflx(1)'
         endif
c
c-- polar projection
      elseif(iflx==2) then
         if(idex>0) then
c-- uniform polar array
            flx_nmu = idex
            allocate(flx_mu(flx_nmu+1))
            help = 2d0/flx_nmu
            forall(i=1:flx_nmu+1) flx_mu(i) = -1d0+(i-1)*help
         else
            stop 'generate_fluxgrid: invalid nflx(2)'
         endif
c
c-- azimuthal angle
      else
         if(idex>0) then
c-- uniform azimuthal array
            flx_nom = idex
            allocate(flx_om(flx_nom+1))
            help = pc_pi2/flx_nom
            forall(i=1:flx_nom+1) flx_om(i) = (i-1)*help
         else
            stop 'generate_fluxgrid: invalid nflx(3)'
         endif
      endif
c
      end subroutine generate_fluxgrid
c
c
      subroutine flux_dealloc
      deallocate(flx_wl)!{{{
      deallocate(flx_mu)
      deallocate(flx_om)
      deallocate(flx_luminos)
      deallocate(flx_lumdev)
      deallocate(flx_lumnum)
c-- grey gamma flux
      deallocate(flx_gamluminos)
      deallocate(flx_gamlumdev)
      deallocate(flx_gamlumnum)!}}}
      end subroutine flux_dealloc
c
c
      end module fluxmod
