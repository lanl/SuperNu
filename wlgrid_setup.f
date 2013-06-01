      subroutine wlgrid_setup
c     -----------------------
      use gasgridmod
      use inputparmod
      use physconstmod
      use miscmod
      implicit none
************************************************************************
* setup wavelength grid. Either read from file or use simple exponential
* spacing
************************************************************************
      real*8, allocatable :: wlstore(:) 
      real*8 :: help
      integer :: i,ng,ngm,nrm,irr
c
c-- read wavelength grid from file
      if(in_iswlread) then
       call read_wlgrid
      else
!--Ryan W: allocation moved here from gasgridmod--
       allocate(gas_wl(in_ng))
!----- --------------------------------------------
       forall(i=1:in_ng) gas_wl(i) =
     &   in_wlmin*(in_wlmax/dble(in_wlmin))**((i-1d0)/(in_ng-1d0))
       gas_dwl = gas_wl*log(in_wlmax/dble(in_wlmin)) /
     &   (in_ng-1d0)      !wl grid bin width
c-- sa nity test
       help = sum(gas_dwl)
       if(abs(help/(in_wlmax-in_wlmin) - 1d0) .gt. 1d-3) then
        call warn('gasgrid_setup','ggrid_dwl not accurate')
        write(6,*) help,in_wlmax-in_wlmin
       endif
      endif
c
      end subroutine wlgrid_setup
