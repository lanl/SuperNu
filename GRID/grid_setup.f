      subroutine grid_setup
c     ---------------------
      use gridmod
      use inputstrmod
      use physconstmod
      implicit none
************************************************************************
* Setup the grid on the computational domain
************************************************************************
      logical :: lexist
c
c-- agnostic grid setup
      grd_xarr = str_xleft
      grd_yarr = str_yleft
      grd_zarr = str_zleft
c
c-- sanity check
      if(grd_isvelocity) then
       if(maxval(abs(grd_xarr))>pc_c) stop 'grid_setup: grd_xarr > pc_c'
       if(maxval(abs(grd_yarr))>pc_c) stop 'grid_setup: grd_yarr > pc_c'
       if(maxval(abs(grd_zarr))>pc_c) stop 'grid_setup: grd_zarr > pc_c'
      endif
c-- sanity check
      select case(grd_igeom)
      case(1)
       if(minval(grd_xarr)<0d0) stop 'grid_setup: grd_xarr < 0'
      case(2)
       if(minval(grd_xarr)<0d0) stop 'grid_setup: grd_xarr < 0'
      endselect
c
c-- read preset temperature profiles
      inquire(file='input.temp',exist=lexist)
      if(lexist) call read_temp_preset
c
      end subroutine grid_setup
