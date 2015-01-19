subroutine initialnumbers

  use gridmod
  use sourcemod
  use totalsmod
  use particlemod
  implicit none

!##################################################
  !This subroutine computes the distribution of initial particle energy
  !before the first time step.  A fraction of the total initial particle
  !number is given to each cell based on the amount of inital radiative
  !energy profile.
!##################################################
  real*8 :: etotinit
  
  src_ninitnew = 0
!
  grd_nvolinit = 0
  grd_evolinit = 0d0
  etotinit = 0d0
!
  tot_eext = 0d0
!
  call analytic_initial

  etotinit = sum(grd_evolinit)
!
  tot_eext = etotinit
!
  if(etotinit > 0d0) then
     grd_nvolinit = nint(grd_evolinit*src_ninit/etotinit) !+50
  endif

  src_ninitnew = sum(grd_nvolinit)


end subroutine initialnumbers