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
     if(gas_opacanaltype=='pick') then
!
!-- checking input validity for picket fence
        if(gas_ng/=2) stop 'analytic_initial: gas_ng/=2 and opacanaltype=pick'
        if(gas_isvelocity) stop 'analytic_initial: invalid gas_isvelocity'
!
!-- tstd initial energy profile currently approximation
        if(gas_suol=='tstd') then
           gas_evolinit(1,:)=0d0!gas_ppick(1)*gas_vals2%ur*&
                !gas_vals2%volr*(gas_l0+gas_lr)**3
           gas_evolinit(2,:)=gas_ppick(2)*gas_vals2%ur*&
                gas_vals2%volr*(gas_l0+gas_lr)**3
!           gas_evolinit(2,:)=gas_vals2%ur*&
!                gas_vals2%volr*(gas_l0+gas_lr)**3           
        else
           stop 'analytic_initial: gas_suol/=tstd not available'
        endif
     else
        return
     endif
  elseif(gas_srctype=='heav') then
     return
  elseif(gas_srctype=='strt') then
     return
  elseif(gas_srctype=='manu') then
     call init_manuprofile(tsp_texp)
  else
     stop 'analytic_initial: invalid gas_srctype'
  endif

end subroutine analytic_initial
