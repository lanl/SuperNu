!Pure transport routine

SUBROUTINE transport1(z,g,r,mu,t,E,E0,hyparam,vacnt)

  USE data_mod
  IMPLICIT NONE
  !
  INTEGER(iknd), INTENT(INOUT) :: z, g, hyparam
  REAL(rknd), INTENT(INOUT) :: r, mu, t, E, E0
  LOGICAL, INTENT(INOUT) :: vacnt
  !
  INTEGER(iknd) :: ig, iig
  REAL(rknd) :: r1, r2
  REAL(rknd) :: db, dcol, dcen, d
  REAL(rknd) :: siglabfact, dcollabfact, elabfact 
  REAL(rknd) :: rold, P, denom2, told

  siglabfact = 1.0_rknd - velyes*mu*r/lspeed
  dcollabfact = velno*1.0 + velyes*texp*(1.0_rknd-mu*r/lspeed)

  ! distance to boundary = db
  IF (z == 1) THEN
     db = ABS(SQRT(rarr(z+1)**2-(1.0-mu**2)*r**2)-mu*r)
  ELSEIF (mu < -SQRT(1.0_rknd-(rarr(z)/r)**2)) THEN
     db = ABS(SQRT(rarr(z)**2-(1.0_rknd-mu**2)*r**2)+mu*r)
  ELSE
     db = ABS(SQRT(rarr(z+1)**2-(1.0_rknd-mu**2)*r**2)-mu*r)
  ENDIF

  ! distance to collision = dcol
  IF((1.0_rknd-fcoef(z))*sigmapg(g,z)>0.0_rknd) THEN
     r1 = RAND()
     dcol = ABS(LOG(r1)/((1.0_rknd-fcoef(z))*sigmapg(g,z)*dcollabfact))
  ELSE
     dcol = 3.0*db
  ENDIF
  ! distance to census = dcen
  dcen = ABS(lspeed*(time+dt-t)/(velno*1.0+velyes*texp))
  ! minimum distance = d
  d = MIN(dcol,db,dcen)

  rold = r
  r = SQRT((1.0_rknd-mu**2)*r**2+(d+r*mu)**2)
  told = t
  t = t + (velno*1.0+velyes*texp)*d/lspeed
  mu = (rold*mu+d)/r
  elabfact = 1.0_rknd - velyes*mu*r/lspeed
  Edep(z)=Edep(z)+E*(1.0_rknd-EXP(-fcoef(z)*sigmapg(g,z)*d))*elabfact
  E = E*EXP(-fcoef(z)*sigmapg(g,z)*d)
  IF (E/E0<0.001_rknd) THEN
     vacnt = .TRUE.
     done = .TRUE.
     Edep(z) = Edep(z) + E*elabfact
  ENDIF
  
  IF (d == dcol) THEN
     !r1 = RAND()
     !IF (r1 < fcoef(z)) THEN
     !   vacnt = .TRUE.
     !   done = .TRUE.
     !   Edep(z) = Edep(z) + E*elabfact
     !ELSE
        r1 = RAND()
        mu = 1.0-2.0*r1
        mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
        E = E*elabfact/(1.0-velyes*mu*r/lspeed)
        denom2 = 0.0
        r1 = RAND()
        DO ig = 1, gas_ng
           iig = ig
           IF ((r1>=denom2).AND.(r1<denom2+EmitProbg(ig,z))) EXIT
           denom2 = denom2+EmitProbg(ig,z)
        ENDDO
        g = iig
        IF ((sigmapg(g,z)*drarr(z)*(velno*1.0+velyes*texp)>=5.0_rknd).AND.(puretran.EQV..FALSE.)) THEN
           hyparam = 2
           E = E*(1.0-velyes*r*mu/lspeed)
           E0 = E0*(1.0-velyes*r*mu/lspeed)
        ELSE
           hyparam = 1
        ENDIF
        !Edep(z)=Edep(z)+E*(1.0-velyes*mu*r/lspeed)/(sigmapg(g,z)*dt*lspeed)
     !ENDIF
  ELSEIF (d == db) THEN
     IF (mu>=0.0_rknd) THEN
        IF (z == gas_nr) THEN
           vacnt = .TRUE.
           done = .TRUE.
           Eright = Eright+E !*elabfact
        ! Checking if DDMC region right
        ELSEIF ((sigmapg(g,z+1)*drarr(z+1)*(velno*1.0+velyes*texp)>=5.0_rknd).AND.(puretran.EQV..FALSE.)) THEN
           r1 = RAND()
           mu = (mu-velyes*r/lspeed)/(1.0-velyes*r*mu/lspeed)
           P = PPL(g,z+1)*(1.0+1.5*ABS(mu))
           IF (r1 < P) THEN
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              z = z+1
           ELSE
              r1 = RAND()
              r2 = RAND()
              mu = -MAX(r1,r2)
              mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
           ENDIF
        ! End of check
        ELSE
           z = z+1   
        ENDIF
     ELSE
        IF (z==1) THEN
           IF ((sigmapg(g,z+1)*drarr(z+1)*(velno*1.0+velyes*texp)>=5.0_rknd).AND.(puretran.EQV..FALSE.)) THEN
              r1 = RAND()
              mu = (mu-velyes*r/lspeed)/(1.0-velyes*r*mu/lspeed)
              P = PPL(g,z+1)*(1.0+1.5*ABS(mu))
              IF (r1 < P) THEN
                 hyparam = 2
                 E = E*elabfact
                 E0 = E0*elabfact
                 z = z+1
              ELSE
                 r1 = RAND()
                 r2 = RAND()
                 mu = -MAX(r1,r2)
                 mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
              ENDIF
           ELSE
              z = z+1
           ENDIF
        !IF (z==1) THEN
        !   vacnt = .TRUE.
        !   done = .TRUE.
        !   Eleft = Eleft+E*elabfact
        ! Checking if DDMC region left   
        ELSEIF ((sigmapg(g,z-1)*drarr(z-1)*(velno*1.0+velyes*texp)>=5.0_rknd).AND.(puretran.EQV..FALSE.)) THEN
           r1 = RAND()
           mu = (mu-velyes*r/lspeed)/(1.0-velyes*r*mu/lspeed)
           P = PPR(g,z-1)*(1.0+1.5*ABS(mu))
           IF (r1 < P) THEN
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              z = z-1
           ELSE
              r1 = RAND()
              r2 = RAND()
              mu = MAX(r1,r2)
              mu = (mu+velyes*r/lspeed)/(1.0+velyes*r*mu/lspeed)
           ENDIF
        ! End of check
        ELSE
           z = z-1
        ENDIF
     ENDIF
  ELSEIF (d == dcen) THEN
     done = .TRUE.
     numcensus(z) = numcensus(z)+1
     Erad = Erad + E*elabfact
  ENDIF


END SUBROUTINE transport1
