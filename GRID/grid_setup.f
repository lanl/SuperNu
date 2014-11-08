      subroutine grid_setup
c     ---------------------
      use gasgridmod
      use inputstrmod
      use timestepmod
      implicit none
************************************************************************
* Setup the grid on the computational domain
************************************************************************
c-- agnostic grid setup
      grd_xarr = str_xleft
      grd_yarr = str_yleft
      grd_zarr = str_zleft
c
      end subroutine grid_setup
