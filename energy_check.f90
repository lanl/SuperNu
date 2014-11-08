subroutine energy_check

  use totalsmod
  implicit none

!-----------------------------------------------------
!This subroutine checks that all particle energy 
!(weight) is accounted for from conservation in
!comoving quantities.
!-----------------------------------------------------


  tot_eerror = (tot_eextav-tot_eveloav-tot_erad-tot_emat)/&
       tot_eextav



end subroutine energy_check
