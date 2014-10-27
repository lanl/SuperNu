      subroutine gridvolume(igeom,isvelocity,time)
c     --------------------------------------------
      use gasgridmod
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
      case(1)
       gas_vals2(:,1,1)%vol = pc_pi43*t**3 *
     &   (gas_xarr(2:)**3 - gas_xarr(:gas_nx)**3)
      case(2)
       forall(i=1:gas_nx,j=1:gas_ny)
        gas_vals2(i,j,1)%vol = pc_pi*t**3 *
     &    (gas_yarr(i+1) - gas_yarr(i)) *
     &    (gas_xarr(j+1)**2 - gas_xarr(j)**2)
       endforall
      case(3)
       forall(i=1:gas_nx,j=1:gas_ny,k=1:gas_nz)
        gas_vals2(i,j,k)%vol = t**3 *
     &    (gas_xarr(i+1) - gas_xarr(i)) *
     &    (gas_yarr(j+1) - gas_yarr(j)) *
     &    (gas_zarr(k+1) - gas_zarr(k))
       endforall
      endselect !in_igeom
c
      end subroutine gridvolume
