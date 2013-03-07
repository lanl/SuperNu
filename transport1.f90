!Pure transport routine

SUBROUTINE transport1(z,g,r,mu,t,E,E0,hyparam,vacnt)

  USE gasgridmod
  USE timestepmod
  USE physconstmod
  USE particlemod
  USE inputparmod
  IMPLICIT NONE
  !
  INTEGER, INTENT(INOUT) :: z, g, hyparam
  REAL*8, INTENT(INOUT) :: r, mu, t, E, E0
  LOGICAL, INTENT(INOUT) :: vacnt
  !
  INTEGER :: ig, iig
  REAL*8 :: r1, r2
  REAL*8 :: db, dcol, dcen, d
  REAL*8 :: siglabfact, dcollabfact, elabfact 
  REAL*8 :: rold, P, denom2, told

  siglabfact = 1.0d0 - gas_velyes*mu*r/pc_c
  dcollabfact = gas_velno*1.0 + gas_velyes*tsp_texp*(1.0d0-mu*r/pc_c)

  ! distance to boundary = db
  IF (z == 1) THEN
     db = ABS(SQRT(gas_rarr(z+1)**2-(1.0-mu**2)*r**2)-mu*r)
  ELSEIF (mu < -SQRT(1.0d0-(gas_rarr(z)/r)**2)) THEN
     db = ABS(SQRT(gas_rarr(z)**2-(1.0d0-mu**2)*r**2)+mu*r)
  ELSE
     db = ABS(SQRT(gas_rarr(z+1)**2-(1.0d0-mu**2)*r**2)-mu*r)
  ENDIF

  ! distance to collision = dcol
  IF((1.0d0-gas_fcoef(z))*gas_sigmapg(g,z)>0.0d0) THEN
     r1 = RAND()
     dcol = ABS(LOG(r1)/((1.0d0-gas_fcoef(z))*gas_sigmapg(g,z)*dcollabfact))
  ELSE
     dcol = 3.0*db
  ENDIF
  ! distance to census = dcen
  dcen = ABS(pc_c*(tsp_time+tsp_dt-t)/(gas_velno*1.0+gas_velyes*tsp_texp))
  ! minimum distance = d
  d = MIN(dcol,db,dcen)

  rold = r
  r = SQRT((1.0d0-mu**2)*r**2+(d+r*mu)**2)
  told = t
  t = t + (gas_velno*1.0+gas_velyes*tsp_texp)*d/pc_c
  mu = (rold*mu+d)/r
  elabfact = 1.0d0 - gas_velyes*mu*r/pc_c
  gas_edep(z)=gas_edep(z)+E*(1.0d0-EXP(-gas_fcoef(z)*gas_sigmapg(g,z)*d))*elabfact
  E = E*EXP(-gas_fcoef(z)*gas_sigmapg(g,z)*d)
  IF (E/E0<0.001d0) THEN
     vacnt = .TRUE.
     prt_done = .TRUE.
     gas_edep(z) = gas_edep(z) + E*elabfact
  ENDIF
  
  IF (d == dcol) THEN
     !r1 = RAND()
     !IF (r1 < gas_fcoef(z)) THEN
     !   vacnt = .TRUE.
     !   prt_done = .TRUE.
     !   gas_edep(z) = gas_edep(z) + E*elabfact
     !ELSE
        r1 = RAND()
        mu = 1.0-2.0*r1
        mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
        E = E*elabfact/(1.0-gas_velyes*mu*r/pc_c)
        denom2 = 0.0
        r1 = RAND()
        DO ig = 1, gas_ng
           iig = ig
           IF ((r1>=denom2).AND.(r1<denom2+gas_emitprobg(ig,z))) EXIT
           denom2 = denom2+gas_emitprobg(ig,z)
        ENDDO
        g = iig
        IF ((gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0).AND.(in_puretran.EQV..FALSE.)) THEN
           hyparam = 2
           E = E*(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0*(1.0-gas_velyes*r*mu/pc_c)
        ELSE
           hyparam = 1
        ENDIF
        !gas_edep(z)=gas_edep(z)+E*(1.0-gas_velyes*mu*r/pc_c)/(gas_sigmapg(g,z)*tsp_dt*pc_c)
     !ENDIF
  ELSEIF (d == db) THEN
     IF (mu>=0.0d0) THEN
        IF (z == gas_nr) THEN
           vacnt = .TRUE.
           prt_done = .TRUE.
           prt_eright = prt_eright+E !*elabfact
        ! Checking if DDMC region right
        ELSEIF ((gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
                 .AND.(in_puretran.EQV..FALSE.)) THEN
           r1 = RAND()
           mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
           P = gas_ppl(g,z+1)*(1.0+1.5*ABS(mu))
           IF (r1 < P) THEN
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              z = z+1
           ELSE
              r1 = RAND()
              r2 = RAND()
              mu = -MAX(r1,r2)
              mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           ENDIF
        ! End of check
        ELSE
           z = z+1   
        ENDIF
     ELSE
        IF (z==1) THEN
           IF ((gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
                   .AND.(in_puretran.EQV..FALSE.)) THEN
              r1 = RAND()
              mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
              P = gas_ppl(g,z+1)*(1.0+1.5*ABS(mu))
              IF (r1 < P) THEN
                 hyparam = 2
                 E = E*elabfact
                 E0 = E0*elabfact
                 z = z+1
              ELSE
                 r1 = RAND()
                 r2 = RAND()
                 mu = -MAX(r1,r2)
                 mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
              ENDIF
           ELSE
              z = z+1
           ENDIF
        !IF (z==1) THEN
        !   vacnt = .TRUE.
        !   prt_done = .TRUE.
        !   prt_eleft = prt_eleft+E*elabfact
        ! Checking if DDMC region left   
        ELSEIF ((gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
             .AND.(in_puretran.EQV..FALSE.)) THEN
           r1 = RAND()
           mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
           P = gas_ppr(g,z-1)*(1.0+1.5*ABS(mu))
           IF (r1 < P) THEN
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              z = z-1
           ELSE
              r1 = RAND()
              r2 = RAND()
              mu = MAX(r1,r2)
              mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           ENDIF
        ! End of check
        ELSE
           z = z-1
        ENDIF
     ENDIF
  ELSEIF (d == dcen) THEN
     prt_done = .TRUE.
     gas_numcensus(z) = gas_numcensus(z)+1
     prt_erad = prt_erad + E*elabfact
  ENDIF


END SUBROUTINE transport1
