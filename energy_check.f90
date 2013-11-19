subroutine energy_check

  use gasgridmod
  use timestepmod
  use physconstmod
  use manufacmod
  use inputparmod
  implicit none

!-----------------------------------------------------
!This subroutine checks that all particle energy 
!(weight) is accounted for from conservation in
!comoving quantities.
!-----------------------------------------------------


  gas_eext = gas_eext-gas_eleft-gas_eright
  gas_eerror = (gas_eext-gas_evelo-gas_erad-gas_emat)/&
       (gas_eext-gas_evelo)


end subroutine energy_check
