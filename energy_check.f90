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


  gas_eerror = (gas_eextav-gas_eveloav-gas_erad-gas_emat)/&
       gas_eextav



end subroutine energy_check
