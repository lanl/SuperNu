      subroutine grid_update(t)
c     -------------------------
      use gridmod
      implicit none
      real*8,intent(in) :: t
************************************************************************
* Setup the grid on the whole domain
************************************************************************
c-- volume 
      call grid_volume(grd_igeom,grd_isvelocity,t)
c
      end subroutine grid_update
