SUBROUTINE material_update

  USE gasgridmod
  USE timestepmod
  USE physconstmod
  IMPLICIT NONE

!##################################################
  !This subroutine updates the material state.  At the moment
  !it updates material temperature and the approximate amount
  !of gamma ray energy introduced in the time step.
!##################################################
  INTEGER :: ir
  REAL*8 :: dtemp, Um, expfact, tauNi, tauCo

  gas_emat = 0.0
  DO ir = 1, gas_nr
     dtemp = gas_edep(ir)*3.0/(4.0*pc_pi*gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**3))
     dtemp = (dtemp-tsp_dt*gas_fcoef(ir)*gas_sigmap(ir)*pc_c*gas_vals2(ir)%ur)/gas_vals2(ir)%bcoef
     !WRITE(*,*) dtemp
     gas_vals2(ir)%tempkev = gas_vals2(ir)%tempkev+dtemp
     !gas_vals2(ir)%ur=dtemp/(tsp_dt*pc_c*gas_sigmap(ir))
     !gas_vals2(ir)%tempkev = (gas_vals2(ir)%ur/pc_acoef)**(0.25d0)
     !gas_vals2(ir)%bcoef = 2.0*pc_acoef*gas_vals2(ir)%tempkev**3
     gas_vals2(ir)%ur = pc_acoef*gas_vals2(ir)%tempkev**4
     Um = gas_vals2(ir)%bcoef*gas_vals2(ir)%tempkev
     gas_emat = gas_emat + Um*4.0*pc_pi*gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**3)/3.0
     !Calculating expansion losses (if any)
     expfact = gas_velno*1.0+gas_velyes*tsp_texp/(tsp_texp+tsp_dt)
     gas_vals2(ir)%rho = gas_vals2(ir)%rho*expfact**3
     gas_vals2(ir)%bcoef = gas_vals2(ir)%bcoef*expfact**3
     !gas_edep(ir) = gas_edep(ir)*3.0/(4.0*pc_pi*gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**3))
  ENDDO
  tauCo = 111.3d0*86400.0d0
  tauNi = 8.8d0*86400.0d0
  
  gas_nidecay = (1.6022e-6)*1.87*(1.0d0-EXP(-(tsp_time+tsp_dt)/tauNi))
  gas_nidecay = gas_nidecay+(1.6022e-6)*1.87*tauCo*(1.0d0-EXP(-(tsp_time+tsp_dt)/tauCo))/(tauCo-tauNi)
  gas_nidecay = gas_nidecay-(1.6022e-6)*1.87*(1.0d0-EXP(-tsp_time/tauNi))
  gas_nidecay = gas_nidecay-(1.6022e-6)*1.87*tauCo*(1.0d0-EXP(-tsp_time/tauCo))/(tauCo-tauNi)
  gas_nidecay = gas_nidecay/tsp_dt
  !WRITE(*,*) gas_nidecay

END SUBROUTINE material_update
