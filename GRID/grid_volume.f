*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine grid_volume(igeom,isvelocity,time)
c     ---------------------------------------------
      use gridmod
      use physconstmod
      implicit none
      integer,intent(in) :: igeom
      logical,intent(in) :: isvelocity
      real*8,intent(in) :: time
************************************************************************
* calculate volumes for all grid cells for expansion time t
************************************************************************
      integer :: i,j,k,l
      real*8 :: t
      real*8,allocatable :: vol(:,:,:) !too big for the stack
c
      allocate(vol(grd_nx,grd_ny,grd_nz))
c
      if(isvelocity) then
       t = time
      else
       t = 1d0
      endif
c
      select case(igeom)
      case(1,11)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
        vol(i,j,k) = t**3 *
     &    (grd_xarr(i+1)**3 - grd_xarr(i)**3) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_zarr(k+1) - grd_zarr(k))/3d0
       endforall
      case(2)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
        vol(i,j,k) = t**3 *
     &    (grd_xarr(i+1)**2 - grd_xarr(i)**2) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    .5d0*(grd_zarr(k+1) - grd_zarr(k))
       endforall
      case(3)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
        vol(i,j,k) = t**3 *
     &    (grd_xarr(i+1) - grd_xarr(i)) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_zarr(k+1) - grd_zarr(k))
       endforall
      case default
       stop 'grid_volume: invalid igeom'
      endselect !igeom
c
c-- map to compressed grid
      grd_vol = 0d0
      do k=1,grd_nz
      do j=1,grd_ny
      do i=1,grd_nx
       l = grd_icell(i,j,k)
       grd_vol(l) = grd_vol(l) + vol(i,j,k) !multiple void cells are linked to the dummy cell
      enddo
      enddo
      enddo
c
      deallocate(vol)
c
      end subroutine grid_volume
c vim: fdm=marker
