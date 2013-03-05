!Pure diffusion routine

SUBROUTINE diffusion1alt(z,g,r,mu,t,E,E0,hyparam,vacnt)

  USE data_mod
  IMPLICIT NONE
  !
  INTEGER(iknd), INTENT(INOUT) :: z, g, hyparam
  REAL(rknd), INTENT(INOUT) :: r, mu, t, E, E0
  LOGICAL, INTENT(INOUT) :: vacnt
  !
  INTEGER(iknd) :: ig, iig
  REAL(rknd) :: r1, r2
  REAL(rknd) :: denom, denom2
  REAL(rknd) :: ddmct, tauA, tauL, tauR, tauS, tcensus
  REAL(rknd), DIMENSION(ng) :: PDFg

  !denom = sigmaL(g,z)+sigmaR(g,z)+fcoef(z)*sigmapg(g,z)
  !denom = denom+(1.0-EmitProbg(g,z))*(1.0-fcoef(z))*sigmapg(g,z)
  r1 = RAND()
  tauA = ABS(LOG(r1)/(lspeed*fcoef(z)*sigmapg(g,z)))
  r1 = RAND()
  tauS = ABS(LOG(r1)/(lspeed*(1.0-EmitProbg(g,z))*(1.0-fcoef(z))*sigmapg(g,z)))
  r1 = RAND()
  tauR = ABS(LOG(r1)/(lspeed*sigmaR(g,z)))
  r1 = RAND()
  tauL = ABS(LOG(r1)/(lspeed*sigmaL(g,z)))
  tcensus = time+dt-t

  ddmct = MIN(tauA,tauS,tauR,tauL,tcensus)
  E = E*(velno*1.0+velyes*EXP(-ddmct/texp))
  E0 = E0*(velno*1.0+velyes*EXP(-ddmct/texp))
  t = t+ddmct
  !WRITE(*,*) ddmct, tau, tcensus
  !IF (ddmct == tau) THEN
     !r1 = RAND()
     !PR = sigmaR(g,z)/denom
     !PL = sigmaL(g,z)/denom
     !PA = fcoef(z)*sigmapg(g,z)/denom
  IF (ddmct == tauL) THEN
     IF (z == 1) THEN
        !WRITE(*,*) 'Non-physical left leakage'
        vacnt = .TRUE.
        done = .TRUE.
        Eleft = Eleft+E
     !ELSEIF (sigmapg(g,z-1)*drarr(z-1)*(velno*1.0+velyes*texp)>=5.0_rknd) THEN
     !   z = z-1
     ELSE
     !   hyparam = 1
     !   r = rarr(z)
        z = z-1
     !   r1 = RAND()
     !   r2 = RAND()
     !   mu = -MAX(r1,r2)
     !   mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
     !   E = E/(1.0-velyes*r*mu/lspeed)
     !   E0 = E0/(1.0-velyes*r*mu/lspeed)
     ENDIF
  ELSEIF (ddmct == tauR) THEN
     IF (z == nr) THEN
        vacnt = .TRUE.
        done = .TRUE.
        Eright = Eright+E
     !ELSEIF (sigmapg(g,z+1)*drarr(z+1)*(velno*1.0+velyes*texp)>=5.0_rknd) THEN
     !   z = z+1
     ELSE
     !   hyparam = 1
     !   r = rarr(z+1)
        z = z+1
     !   r1 = RAND()
     !   r2 = RAND()
     !   mu = MAX(r1,r2)
     !   mu = (mu+velyes*r/lspeed)/(1.0+r*mu/lspeed)
     !   E = E/(1.0-velyes*r*mu/lspeed)
     !   E0 = E0/(1.0-velyes*r*mu/lspeed)
     ENDIF
  ELSEIF (ddmct == tauA) THEN
     vacnt = .TRUE.
     done = .TRUE.
     Edep(z) = Edep(z)+E
  ELSEIF (ddmct == tauS) THEN
     denom2 = sigmap(z)-Ppick(g)*sigmapg(g,z)
     DO ig = 1, ng
        PDFg(ig) = EmitProbg(ig,z)*sigmap(z)/denom2 
     ENDDO
     PDFg(g)=0.0
     denom2 = 0.0
     r1 = RAND()
     DO ig = 1, ng
        iig = ig
        IF (r1>=denom2.AND.r1<denom2+PDFg(ig)) EXIT
        denom2 = denom2+PDFg(ig)
     ENDDO
     g = iig
     !IF (sigmapg(g,z)*drarr(z)*(velno*1.0+velyes*texp)>=5.0_rknd) THEN
     !   hyparam = 2
     !ELSE
     !   hyparam = 1
     !   r1 = RAND()
     !   mu = 1.0-2.0*r1
     !   r1 = RAND()
     !   r = r1*rarr(z+1)+(1.0-r1)*rarr(z) !(r1*rarr(z+1)**3+(1.0-r1)*rarr(z)**3)**(1.0/3.0)
     !   mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
     !   E = E/(1.0-velyes*mu*r/lspeed)
     !   E0 = E0/(1.0-velyes*mu*r/lspeed)
     !ENDIF
  ELSE
     done = .TRUE.
     numcensus(z)=numcensus(z)+1
     Erad = Erad+E
  ENDIF

END SUBROUTINE diffusion1alt
