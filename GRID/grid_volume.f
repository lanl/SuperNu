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
      integer :: i,j,k
      real*8 :: t
c
      if(isvelocity) then
       t = time
      else
       t = 1d0
      endif
c
      select case(igeom)
      case(1,4)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
     &  grd_vol(i,j,k) = t**3 *
     &    (grd_xarr(i+1)**3 - grd_xarr(i)**3) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_zarr(k+1) - grd_zarr(k))/3d0
      case(2)
       forall(i=1:grd_nx,j=1:grd_ny)
        grd_vol(i,j,1) = pc_pi*t**3 *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_xarr(i+1)**2 - grd_xarr(i)**2)
       endforall
      case(3)
       forall(i=1:grd_nx,j=1:grd_ny,k=1:grd_nz)
        grd_vol(i,j,k) = t**3 *
     &    (grd_xarr(i+1) - grd_xarr(i)) *
     &    (grd_yarr(j+1) - grd_yarr(j)) *
     &    (grd_zarr(k+1) - grd_zarr(k))
       endforall
      case default
       stop 'grid_volume: invalid igeom'
      endselect !in_igeom
c
      end subroutine grid_volume
