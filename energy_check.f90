subroutine energy_check

  use totalsmod
  implicit none

!-----------------------------------------------------
!This subroutine checks that all particle energy 
!(weight) is accounted for from conservation in
!comoving quantities.
!-----------------------------------------------------


  tot_eerror = (tot_eext-tot_evelo-tot_eout - &
       tot_erad-tot_emat)/tot_eext



end subroutine energy_check
