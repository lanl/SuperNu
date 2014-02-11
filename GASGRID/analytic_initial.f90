subroutine analytic_initial

  use gasgridmod
  use inputparmod
  use physconstmod
  use timestepmod
  use manufacmod
  use profiledatamod
  implicit none

  integer :: ir, ig
  real*8 :: specint
  real*8 :: trad(gas_nr)

!###############################################
! This subroutines attributes radiation energy to
! each cell and group depeding on user specification
! of gas_srctype
!###############################################
!
  gas_evolinit = 0d0
!
!-- initial radiation energy
  if(in_tradinittype=='prof') then
    if(prof_nttrad==0) stop 'analytic_initial: no trad profile data'
    trad = trad_profile(tsp_t)
    write(6,*) 'debug'
    write(6,'(1p,8e12.4)') trad
  elseif(in_tradinittype=='unif') then
    trad = in_tempradinit
  else
    stop 'analytic_initial: invalid in_tradinittype'
  endif
!
!-- map radiation temperature to gas_evolinit
  if(.not.gas_isvelocity) then
    do ig = 1, gas_ng
      gas_evolinit(ig,:)=pc_acoef*trad**4&
        *(1d0/gas_wl(ig)-1d0/gas_wl(ig+1))/&
        (1d0/gas_wl(1)-1d0/gas_wl(gas_ng+1))&
        *gas_vals2%volr*(gas_l0+gas_lr)**3
    enddo
  else
    do ig = 1, gas_ng
      gas_evolinit(ig,:)=pc_acoef*trad**4&
        *(1d0/gas_wl(ig)-1d0/gas_wl(ig+1))/&
        (1d0/gas_wl(1)-1d0/gas_wl(gas_ng+1))&
        *gas_vals2%volr*(tsp_t*gas_velout)**3
    enddo
  endif
!--
!
!-- source specific initial conditions (overrides gas_inittyp)
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
     call init_manuprofile(tsp_t)
  else
     stop 'analytic_initial: invalid gas_srctype'
  endif

end subroutine analytic_initial
