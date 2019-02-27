*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
c vim: fdm=marker
