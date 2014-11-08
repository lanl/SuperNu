      subroutine wlgrid_setup(ng)
c     ---------------------------
      use gasmod
      use inputparmod
      implicit none
      integer,intent(out) :: ng
************************************************************************
* setup wavelength grid. Either read from file or use simple logarithmic
* spacing
************************************************************************
      integer :: ig
c
c-- read wavelength grid from file
      if(in_ng==0) then
       call read_wlgrid(ng)
      else
       ng = in_ng
       allocate(gas_wl(ng+1))
       forall(ig=1:ng+1) gas_wl(ig) =
     &   in_wlmin*(in_wlmax/dble(in_wlmin))**((ig-1d0)/ng)
      endif
c
      contains
c
      subroutine read_wlgrid(ng)
c     --------------------------!{{{
      use gasmod
      use inputparmod
      implicit none
      integer,intent(out) :: ng
************************************************************************
* read wavelength grid from file
*
* EXAMPLE FILE LAYOUT:
* --------------------
* #wavelength bin boundaries. units: [cm]
* #ncell ngroupmax
* #icell ngroup wlbound
* 10 5
*  1 1 1e-5 32e-5
*  2 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
* 
************************************************************************
      real*8 :: help
      real*8,allocatable :: wlstore(:)
      integer :: ngm,nrm,irr,l,ig
c
      open(4,file='input.wlgrid',status='old')
c
c-- strip header
      do l=1,3
       read(4,*)
      enddo
c-- read dimensions
      read(4,*) nrm,ngm !not used
c
      do l=1,in_wldex-1
       read(4,*)
      enddo
      read(4,*) irr, ng
c
      allocate(gas_wl(ng+1))
      allocate(wlstore(ng+3))
      rewind(4)
      do l=1,in_wldex+3
       read(4,*)
      enddo
      read(4,*) wlstore(:)
      close(4)
      gas_wl = wlstore(3:)
      deallocate(wlstore)
c
c-- check monotonicity
      help = 0d0
      do ig=1,ng+1
       if(gas_wl(ig)<=help) stop 'read_wlgrid: wlgrid not increasing'
       help = gas_wl(ig)
      enddo
c
      write(6,*)
      write(6,*) 'read_wlgrid: ng=',ng
!     write(6,*) 'wavelength grid in [cm]'
!     write(6,*) gas_wl
c!}}}
      end subroutine read_wlgrid
c
      end subroutine wlgrid_setup
