      subroutine grid_update(t)
c     -------------------------
      use gridmod
      use gasgridmod
      implicit none
      real*8,intent(in) :: t
************************************************************************
* Setup the grid on the whole domain
************************************************************************
c-- volume 
      call grid_volume(grd_igeom,gas_isvelocity,t)

      end subroutine grid_update
