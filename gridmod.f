      module gridmod
c     --------------
      implicit none
c
      integer :: grd_igeom = 0
c
      integer :: grd_nx=0
      integer :: grd_ny=0
      integer :: grd_nz=0
c
      save
c
      contains
c
      subroutine grid_init(igeom,ndim)
c     --------------------------------!{{{
      use gasgridmod, only:gas_xarr,gas_yarr,gas_zarr
      implicit none
      integer,intent(in) :: igeom
      integer,intent(in) :: ndim(3)
c
      grd_igeom = igeom
c
      grd_nx = ndim(1)
      grd_ny = ndim(2)
      grd_nz = ndim(3)

      allocate(gas_xarr(grd_nx+1))
      allocate(gas_yarr(grd_ny+1))
      allocate(gas_zarr(grd_nz+1))
c!}}}
      end subroutine grid_init
c
      end module gridmod
