subroutine analytic_initial

  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2, x3, x4

!###############################################
! This subroutines attributes radiation energy to
! each cell and group depeding on user specification
! of gas_srctype
!###############################################
!
  gas_evolinit = 0d0
!
!-- currently only supplying nonzero for gas_srctype=manu
  if(gas_srctype=='none') then
     return
  elseif(gas_srctype=='heav') then
     return
  elseif(gas_srctype=='strt') then
     return
  elseif(gas_srctype=='manu') then
     call init_manuprofile(tsp_texp)
  else
     stop 'analytic_initial'
  endif

end subroutine analytic_initial
