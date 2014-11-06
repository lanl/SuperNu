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
      gas_xarr = str_xleft
      gas_yarr = str_yleft
      gas_zarr = str_zleft
c
      end subroutine grid_setup
