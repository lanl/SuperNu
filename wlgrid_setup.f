      subroutine wlgrid_setup(ngin,ng)
c     ---------------------------------
      use gasgridmod
      use inputparmod
      use physconstmod
      use miscmod
      implicit none
      integer,intent(in) :: ngin
      integer,intent(out) :: ng
************************************************************************
* setup wavelength grid. Either read from file or use simple exponential
* spacing
************************************************************************
      real*8, allocatable :: wlstore(:) 
      real*8 :: help
      integer :: i,ngm,nrm,irr
c
c-- read wavelength grid from file
      if(in_iswlread) then
       call read_wlgrid(ng)
      else
       ng = ngin
!--Ryan W: allocation moved here from gasgridmod--
       allocate(gas_wl(ng))
!----- --------------------------------------------
       forall(i=1:ng) gas_wl(i) =
     &   in_wlmin*(in_wlmax/dble(in_wlmin))**((i-1d0)/(ng-1d0))
      endif
c
      end subroutine wlgrid_setup
