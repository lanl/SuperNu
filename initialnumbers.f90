subroutine initialnumbers

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine computes the distribution of initial particle energy
  !before the first time step.  A fraction of the total initial particle
  !number is given to each cell based on the amount of inital radiative
  !energy profile.
!##################################################
  
  integer :: ir, ig
  real*8 :: etotinit
  
  prt_ninitnew=0
!
  gas_nvolinit = 0
  gas_evolinit = 0d0
  etotinit = 0d0
!
  gas_eext = 0d0
!
  call analytic_initial

  etotinit = sum(gas_evolinit)
!
  gas_eext = etotinit
!
  do ir=1,gas_nr
     gas_nvolinit(ir)=nint(sum(gas_evolinit(:,ir))*prt_ninit/etotinit) !+50
  enddo

  prt_ninitnew = sum(gas_nvolinit)


end subroutine initialnumbers
