subroutine gasgrid_upold

  use mpimod, only:nmpi
  use gasgridmod
  use timestepmod
  use physconstmod
  use timingmod
  use inputparmod
  implicit none
  
!################################################
  !Simple version of gasgrid_update for debug purposes
!################################################
  integer :: ir
  real*8 :: help

!-- timing
  real :: t0,t1
!
!-- begin
!
      write(7,*)
      write(7,*) 'update gas grid:'
      write(7,*) '---------------------------'
      if(tsp_tn==1) then
       write(6,*)
       write(6,*) 'update gas grid:'
       write(6,*) '---------------------------'
      endif
!
!
      call time(t0)
!


!-- update volume and density
!============================
      if(gas_isvelocity) then
       help = gas_velout*tsp_texp
      else
       help = gas_lr
      endif
      !gas_vals2%vol = gas_vals2%volr*(gas_velout*tsp_tcenter)**3 !volume in cm^3!{{{
      gas_vals2%vol = gas_vals2%volr*help**3 !volume in cm^3
      gas_vals2%volcrp = gas_vals2%vol !effective volume in cm^3
!
!-- density
      gas_vals2%rho = gas_vals2%mass/gas_vals2%vol
      !gas_vals2%bcoef = 2.0*pc_acoef*gas_vals2%tempkev**3
      !gas_vals2%bcoef = 0.4*(1.e12*gas_vals2%rho)*580.25d0 !currently performed in xsection (power law)
!
!-- keep track of temperature evolution
      gas_temphist(:,tsp_tn) = gas_vals2%temp!}}}
!
!
!
!-- reset counters
!=================
      gas_erad = 0.0   !Total radiation energy
      gas_eint = 0.0   !Total internal energy

      gas_eraddensg =0d0 !radiation density field

      call time(t1)
      call timereg(t_gasupd,t1-t0)


end subroutine gasgrid_upold
