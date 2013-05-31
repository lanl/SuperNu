      subroutine convert_cap2capros
c     -----------------------------
      use gasgridmod
      implicit none
************************************************************************
* Take the cell-centered planck-averaged group opacities and convert
* them to face-centered rosseland opacities.
************************************************************************
      gas_caprosl(:,2:) = 1d0/(1d0/gas_caprosl(:,2:) +
     &  .5d0*(1d0/gas_cap(:,2:) + 1d0/gas_cap(:,:gas_nr-1)))
      gas_caprosl(:,1) = 1d0/(1d0/gas_caprosl(:,1) +
     &  1d0/gas_cap(:,1))
      gas_caprosr(:,:gas_nr-1) = 1d0/(1d0/gas_caprosr(:,:gas_nr-1) +
     &  .5d0*(1d0/gas_cap(:,:gas_nr-1) + 1d0/gas_cap(:,2:)))
      gas_caprosr(:,gas_nr) = 1d0/(1d0/gas_caprosr(:,gas_nr) +
     &  1d0/gas_cap(:,gas_nr))
c
      write(0,*) 'TEST',gas_caprosl(1,2),gas_cap(1,1),gas_cap(1,2)
c
      end subroutine convert_cap2capros
