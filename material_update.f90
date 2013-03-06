SUBROUTINE material_update

  USE data_mod
  IMPLICIT NONE

  INTEGER(iknd) :: ir
  REAL(rknd) :: dtemp, Um, expfact, tauNi, tauCo

  Emat = 0.0
  DO ir = 1, gas_nr
     dtemp = Edep(ir)*3.0/(4.0*pi*dr3arr(ir)*(velno*1.0+velyes*texp**3))
     dtemp = (dtemp-dt*fcoef(ir)*sigmap(ir)*lspeed*Ur(ir))/bcoef(ir)
     !WRITE(*,*) dtemp
     Temp(ir) = Temp(ir)+dtemp
     !Ur(ir)=dtemp/(dt*lspeed*sigmap(ir))
     !Temp(ir) = (Ur(ir)/a_coef)**(0.25_rknd)
     !bcoef(ir) = 2.0*a_coef*Temp(ir)**3
     Ur(ir) = a_coef*Temp(ir)**4
     Um = bcoef(ir)*Temp(ir)
     Emat = Emat + Um*4.0*pi*dr3arr(ir)*(velno*1.0+velyes*texp**3)/3.0
     !Calculating expansion losses (if any)
     expfact = velno*1.0+velyes*texp/(texp+dt) !(in_lr+rarr(gas_nr+1)*time)/(in_lr+rarr(gas_nr+1)*(time+dt))
     rhoarr(ir) = rhoarr(ir)*expfact**3
     bcoef(ir) = bcoef(ir)*expfact**3
     !Edep(ir) = Edep(ir)*3.0/(4.0*pi*dr3arr(ir)*(velno*1.0+velyes*texp**3))
  ENDDO
  tauCo = 111.3_rknd*86400.0_rknd
  tauNi = 8.8_rknd*86400.0_rknd
  
  nidecay = (1.6022e-6)*1.87*(1.0_rknd-EXP(-(time+dt)/tauNi))
  nidecay = nidecay+(1.6022e-6)*1.87*tauCo*(1.0_rknd-EXP(-(time+dt)/tauCo))/(tauCo-tauNi)
  nidecay = nidecay-(1.6022e-6)*1.87*(1.0_rknd-EXP(-time/tauNi))
  nidecay = nidecay-(1.6022e-6)*1.87*tauCo*(1.0_rknd-EXP(-time/tauCo))/(tauCo-tauNi)
  nidecay = nidecay/dt
  !WRITE(*,*) nidecay

END SUBROUTINE material_update
