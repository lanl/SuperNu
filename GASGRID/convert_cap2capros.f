      subroutine convert_cap2capros
c     -----------------------------
      use gasgridmod
      implicit none
************************************************************************
* Take the cell-centered planck-averaged group opacities and convert
* them to face-centered rosseland opacities.
*
* Note: The default value of gas_capros is huge(real*8), so that in the
* sum below its contribution vanishes.
************************************************************************
c
      gas_caprosl(:,2:) = 1d0/(1d0/gas_caprosl(:,2:) +
     &  .5d0*(1d0/gas_cap(:,2:) + 1d0/gas_cap(:,:gas_nr-1)))
      gas_caprosl(:,1) = 1d0/(1d0/gas_caprosl(:,1) +
     &  1d0/gas_cap(:,1))
      gas_caprosr(:,:gas_nr-1) = 1d0/(1d0/gas_caprosr(:,:gas_nr-1) +
     &  .5d0*(1d0/gas_cap(:,:gas_nr-1) + 1d0/gas_cap(:,2:)))
      gas_caprosr(:,gas_nr) = 1d0/(1d0/gas_caprosr(:,gas_nr) +
     &  1d0/gas_cap(:,gas_nr))
c
      end subroutine convert_cap2capros
