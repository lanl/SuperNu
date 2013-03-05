SUBROUTINE xsections

  USE data_mod
  IMPLICIT NONE

  INTEGER(iknd) :: ir, ig
  REAL(rknd) :: Um, beta, tt, gg, ggg, eps, bb
  ! Here: left=>toward r=0 and right=>outward

  !Interpolating cell boundary temperatures
  Tempb(1)=Temp(1)
  !Tempb(1) = 1.0
  DO ir = 2, nr
     Tempb(ir) = (Temp(ir)**4+Temp(ir-1)**4)/2.0
     Tempb(ir) = Tempb(ir)**0.25
  ENDDO
  Tempb(nr+1)=Temp(nr)

  !Calculating (or loading) opacities (could be done differently)
  !Picket fence (Planck):
  ! Picket-fence problem
  Ppick(1) = 1.0_rknd
  Ppick(2) = 0.0_rknd
  DO ig = 3, ng
     Ppick(ig) = 0.0
  ENDDO
  DO ir = 1, nr
     sigmapg(1,ir) = 0.10*rhoarr(ir) !/Temp(ir)**3
     sigmapg(2,ir) = 0.10*rhoarr(ir) !/Temp(ir)**3
     DO ig = 3, ng
        sigmapg(ig,ir) = 1.0 !/Temp(ir)**3
     ENDDO
     sigmap(ir)=0.0
     DO ig = 1, ng
        sigmap(ir) = sigmap(ir)+Ppick(ig)*sigmapg(ig,ir)
     ENDDO
     Um = bcoef(ir)*Temp(ir)
     beta = 4.0*Ur(ir)/Um
     fcoef(ir) = 1.0/(1.0+alpha*beta*lspeed*dt*sigmap(ir))
     DO ig = 1, ng
        EmitProbg(ig,ir) = Ppick(ig)*sigmapg(ig,ir)/sigmap(ir)
     ENDDO
  ENDDO
  
  !Picket fence (Rosseland (same as Planck for P-fence)):
  DO ir = 1, nr
     sigmargleft(1,ir) = 0.10*rhoarr(ir) !/Tempb(ir)**3
     sigmargleft(2,ir) = 0.10*rhoarr(ir) !/Tempb(ir)**3
     DO ig = 3, ng
        sigmargleft(ig,ir) = 1.0 !/Tempb(ir)**3
     ENDDO
     sigmargright(1,ir) = 0.10*rhoarr(ir) !/Tempb(ir+1)**3
     sigmargright(2,ir) = 0.10*rhoarr(ir) !/Tempb(ir+1)**3
     DO ig = 3, ng
        sigmargright(ig,ir) = 1.0 !/Tempb(ir+1)**3
     ENDDO
  ENDDO

  ! Calculating IMC-to-DDMC leakage probabilities/(angular polynomial)
  ! See Densmore, 2007
  DO ir = 1, nr
     gg = (3.0*fcoef(ir))**0.5
     eps = (4.0/3.0)*gg/(1.0+0.7104*gg)
     DO ig = 1, ng
        !Calculating for leakage from left
        !tt = sigmargleft(ig,ir)*drarr(ir)*(velno*1.0+velyes*texp)
        tt = sigmapg(ig,ir)*drarr(ir)*(velno*1.0+velyes*texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        PPL(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
        !Calculating for leakage from right
        !tt = sigmargright(ig,ir)*drarr(ir)*(velno*1.0+velyes*texp)
        tt = sigmapg(ig,ir)*drarr(ir)*(velno*1.0+velyes*texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        PPR(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
     ENDDO
  ENDDO

  ! Calculating DDMC(-to-IMC) leakage opacities (vacuum outer bound)
  DO ir = 1, nr
     DO ig = 1, ng  
        !Computing left-leakage opacities
        IF (ir==1) THEN
           !sigmaL(ig,ir)=0.5*PPL(ig,ir)/drarr(ir)
           sigmaL(ig,ir)=1.5*PPL(ig,ir)*rarr(ir)**2
           sigmaL(ig,ir)=sigmaL(ig,ir)/(dr3arr(ir)*(velno*1.0+velyes*texp))
        ELSEIF(sigmapg(ig,ir-1)*drarr(ir-1)*(velno*1.0+velyes*texp)<5.0_rknd) THEN
           !sigmaL(ig,ir)=0.5*PPL(ig,ir)/drarr(ir)
           sigmaL(ig,ir)=1.5*PPL(ig,ir)*rarr(ir)**2
           sigmaL(ig,ir)=sigmaL(ig,ir)/(dr3arr(ir)*(velno*1.0+velyes*texp))
        ELSE
           tt = sigmargleft(ig,ir)*drarr(ir)+sigmargright(ig,ir-1)*drarr(ir-1)
           !sigmaL(ig,ir) = 2.0/(3.0*drarr(ir)) 
           sigmaL(ig,ir) = (2.0*rarr(ir)**2)/(dr3arr(ir)*(velno*1.0+velyes*texp**2))
           sigmaL(ig,ir) = sigmaL(ig,ir)/tt
        ENDIF
        !Computing right-leakage opacities
        IF (ir==nr) THEN
           !sigmaR(ig,ir)=0.5*PPR(ig,ir)/drarr(ir)
           sigmaR(ig,ir)=1.5*PPR(ig,ir)*rarr(ir+1)**2
           sigmaR(ig,ir)=sigmaR(ig,ir)/(dr3arr(ir)*(velno*1.0+velyes*texp))
        ELSEIF(sigmapg(ig,ir+1)*drarr(ir+1)*(velno*1.0+velyes*texp)<5.0_rknd) THEN
           !sigmaR(ig,ir)=0.5*PPR(ig,ir)/drarr(ir)
           sigmaR(ig,ir)=1.5*PPR(ig,ir)*rarr(ir+1)**2
           sigmaR(ig,ir)=sigmaR(ig,ir)/(dr3arr(ir)*(velno*1.0+velyes*texp))
        ELSE
           tt = sigmargright(ig,ir)*drarr(ir)+sigmargleft(ig,ir+1)*drarr(ir+1)
           !sigmaR(ig,ir) = 2.0/(3.0*drarr(ir))
           sigmaR(ig,ir) = (2.0*rarr(ir+1)**2)/(dr3arr(ir)*(velno*1.0+velyes*texp**2))
           sigmaR(ig,ir) = sigmaR(ig,ir)/tt
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE xsections
