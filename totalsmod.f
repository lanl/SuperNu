      module totalsmod
c     ----------------
      implicit none
************************************************************************
* energy conservation check quantities
************************************************************************
      real*8 :: tot_eext = 0d0 !time-integrated input energy from external source
      real*8 :: tot_emat = 0d0 !material energy
      real*8 :: tot_erad = 0d0 !census radiation energy
      real*8 :: tot_eout = 0d0 !energy escaped
      real*8 :: tot_evelo = 0d0 !total energy change to rad field from fluid
      real*8 :: tot_eerror = 0d0 !error in integral problem energy
      real*8 :: tot_esurf = 0d0
c
c-- initial external energy
      real*8 :: tot_eext0 = 0d0
c
      save
c
      contains
c
      subroutine totals_startup
c     -------------------------
      implicit none
************************************************************************
* The initial radiation field is estimated from the assumption
* that radiative work losses and energy deposition in the radiation
* field each recieve half of the total energy deposition.
* Once the radiation field is determined from this assumption, these
* values are adopted in the totals terms.
************************************************************************
c-- no radiation is assumed to have escaped
      tot_eout = 0d0
c-- energy decayed prior to first timestep
      tot_eext = tot_eext0
c-- exact balance
      tot_evelo = tot_eext-tot_erad-tot_emat !tot_emat is added to tot_eext in the first time step in gas_update
c
      end subroutine totals_startup
c
c
      subroutine totals_error
c     -----------------------
      implicit none
************************************************************************
* Check that all particle energy (weight) is accounted for from
* conservation in comoving quantities.
************************************************************************
      tot_eerror = (tot_eext-tot_evelo-tot_eout-tot_erad-tot_emat)/
     &  tot_eext
      end subroutine totals_error
c
      end module totalsmod
