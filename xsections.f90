SUBROUTINE xsections

  USE gasgridmod
  USE timestepmod
  USE physconstmod
  USE inputparmod
  IMPLICIT NONE

  INTEGER(iknd) :: ir, ig
  REAL(rknd) :: Um, beta, tt, gg, ggg, eps, bb
  ! Here: left=>toward r=0 and right=>outward

  !Interpolating cell boundary temperatures
  gas_tempb(1)=gas_temp(1)
  !gas_tempb(1) = 1.0
  DO ir = 2, gas_nr
     gas_tempb(ir) = (gas_temp(ir)**4+gas_temp(ir-1)**4)/2.0
     gas_tempb(ir) = gas_tempb(ir)**0.25
  ENDDO
  gas_tempb(gas_nr+1)=gas_temp(gas_nr)

  !Calculating (or loading) opacities (could be prt_done differently)
  !Picket fence (Planck):
  ! Picket-fence problem
  Ppick(1) = 1.0_rknd
  Ppick(2) = 0.0_rknd
  DO ig = 3, gas_ng
     Ppick(ig) = 0.0
  ENDDO
  DO ir = 1, gas_nr
     gas_sigmapg(1,ir) = 0.10*gas_rhoarr(ir) !/gas_temp(ir)**3
     gas_sigmapg(2,ir) = 0.10*gas_rhoarr(ir) !/gas_temp(ir)**3
     DO ig = 3, gas_ng
        gas_sigmapg(ig,ir) = 1.0 !/gas_temp(ir)**3
     ENDDO
     gas_sigmap(ir)=0.0
     DO ig = 1, gas_ng
        gas_sigmap(ir) = gas_sigmap(ir)+Ppick(ig)*gas_sigmapg(ig,ir)
     ENDDO
     Um = gas_bcoef(ir)*gas_temp(ir)
     beta = 4.0*gas_ur(ir)/Um
     gas_fcoef(ir) = 1.0/(1.0+alpha*beta*pc_c*tsp_dt*gas_sigmap(ir))
     DO ig = 1, gas_ng
        gas_emitprobg(ig,ir) = Ppick(ig)*gas_sigmapg(ig,ir)/gas_sigmap(ir)
     ENDDO
  ENDDO
  
  !Picket fence (Rosseland (same as Planck for P-fence)):
  DO ir = 1, gas_nr
     gas_sigmargleft(1,ir) = 0.10*gas_rhoarr(ir) !/gas_tempb(ir)**3
     gas_sigmargleft(2,ir) = 0.10*gas_rhoarr(ir) !/gas_tempb(ir)**3
     DO ig = 3, gas_ng
        gas_sigmargleft(ig,ir) = 1.0 !/gas_tempb(ir)**3
     ENDDO
     gas_sigmargright(1,ir) = 0.10*gas_rhoarr(ir) !/gas_tempb(ir+1)**3
     gas_sigmargright(2,ir) = 0.10*gas_rhoarr(ir) !/gas_tempb(ir+1)**3
     DO ig = 3, gas_ng
        gas_sigmargright(ig,ir) = 1.0 !/gas_tempb(ir+1)**3
     ENDDO
  ENDDO

  ! Calculating IMC-to-DDMC leakage probabilities/(angular polynomial)
  ! See Densmore, 2007
  DO ir = 1, gas_nr
     gg = (3.0*gas_fcoef(ir))**0.5
     eps = (4.0/3.0)*gg/(1.0+0.7104*gg)
     DO ig = 1, gas_ng
        !Calculating for leakage from left
        !tt = gas_sigmargleft(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        tt = gas_sigmapg(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        gas_ppl(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
        !Calculating for leakage from right
        !tt = gas_sigmargright(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        tt = gas_sigmapg(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        gas_ppr(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
     ENDDO
  ENDDO

  ! Calculating DDMC(-to-IMC) leakage opacities (vacuum outer bound)
  DO ir = 1, gas_nr
     DO ig = 1, gas_ng  
        !Computing left-leakage opacities
        IF (ir==1) THEN
           !gas_sigmal(ig,ir)=0.5*gas_ppl(ig,ir)/gas_drarr(ir)
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp))
        ELSEIF(gas_sigmapg(ig,ir-1)*gas_drarr(ir-1)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0_rknd) THEN
           !gas_sigmal(ig,ir)=0.5*gas_ppl(ig,ir)/gas_drarr(ir)
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp))
        ELSE
           tt = gas_sigmargleft(ig,ir)*gas_drarr(ir)+gas_sigmargright(ig,ir-1)*gas_drarr(ir-1)
           !gas_sigmal(ig,ir) = 2.0/(3.0*gas_drarr(ir)) 
           gas_sigmal(ig,ir) = (2.0*gas_rarr(ir)**2)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmal(ig,ir) = gas_sigmal(ig,ir)/tt
        ENDIF
        !Computing right-leakage opacities
        IF (ir==gas_nr) THEN
           !gas_sigmar(ig,ir)=0.5*gas_ppr(ig,ir)/gas_drarr(ir)
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp))
        ELSEIF(gas_sigmapg(ig,ir+1)*gas_drarr(ir+1)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0_rknd) THEN
           !gas_sigmar(ig,ir)=0.5*gas_ppr(ig,ir)/gas_drarr(ir)
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp))
        ELSE
           tt = gas_sigmargright(ig,ir)*gas_drarr(ir)+gas_sigmargleft(ig,ir+1)*gas_drarr(ir+1)
           !gas_sigmar(ig,ir) = 2.0/(3.0*gas_drarr(ir))
           gas_sigmar(ig,ir) = (2.0*gas_rarr(ir+1)**2)/(gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmar(ig,ir) = gas_sigmar(ig,ir)/tt
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE xsections
