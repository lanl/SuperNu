!Pure diffusion routine

SUBROUTINE diffusion1alt(z,g,r,mu,t,E,E0,hyparam,vacnt)

  USE gasgridmod
  USE particlemod
  USE timestepmod
  USE inputparmod
  USE physconstmod
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
  REAL(rknd), DIMENSION(gas_ng) :: PDFg

  !denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  !denom = denom+(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)
  r1 = RAND()
  tauA = ABS(LOG(r1)/(pc_c*gas_fcoef(z)*gas_sigmapg(g,z)))
  r1 = RAND()
  tauS = ABS(LOG(r1)/(pc_c*(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)))
  r1 = RAND()
  tauR = ABS(LOG(r1)/(pc_c*gas_sigmar(g,z)))
  r1 = RAND()
  tauL = ABS(LOG(r1)/(pc_c*gas_sigmal(g,z)))
  tcensus = tsp_time+tsp_dt-t

  ddmct = MIN(tauA,tauS,tauR,tauL,tcensus)
  E = E*(gas_velno*1.0+gas_velyes*EXP(-ddmct/tsp_texp))
  E0 = E0*(gas_velno*1.0+gas_velyes*EXP(-ddmct/tsp_texp))
  t = t+ddmct
  !WRITE(*,*) ddmct, tau, tcensus
  !IF (ddmct == tau) THEN
     !r1 = RAND()
     !PR = gas_sigmar(g,z)/denom
     !PL = gas_sigmal(g,z)/denom
     !PA = gas_fcoef(z)*gas_sigmapg(g,z)/denom
  IF (ddmct == tauL) THEN
     IF (z == 1) THEN
        !WRITE(*,*) 'Non-physical left leakage'
        vacnt = .TRUE.
        prt_done = .TRUE.
        prt_eleft = prt_eleft+E
     !ELSEIF (gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0_rknd) THEN
     !   z = z-1
     ELSE
     !   hyparam = 1
     !   r = gas_rarr(z)
        z = z-1
     !   r1 = RAND()
     !   r2 = RAND()
     !   mu = -MAX(r1,r2)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*r*mu/pc_c)
     !   E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
     ENDIF
  ELSEIF (ddmct == tauR) THEN
     IF (z == gas_nr) THEN
        vacnt = .TRUE.
        prt_done = .TRUE.
        prt_eright = prt_eright+E
     !ELSEIF (gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0_rknd) THEN
     !   z = z+1
     ELSE
     !   hyparam = 1
     !   r = gas_rarr(z+1)
        z = z+1
     !   r1 = RAND()
     !   r2 = RAND()
     !   mu = MAX(r1,r2)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*r*mu/pc_c)
     !   E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
     ENDIF
  ELSEIF (ddmct == tauA) THEN
     vacnt = .TRUE.
     prt_done = .TRUE.
     gas_edep(z) = gas_edep(z)+E
  ELSEIF (ddmct == tauS) THEN
     denom2 = gas_sigmap(z)-gas_ppick(g)*gas_sigmapg(g,z)
     DO ig = 1, gas_ng
        PDFg(ig) = gas_emitprobg(ig,z)*gas_sigmap(z)/denom2 
     ENDDO
     PDFg(g)=0.0
     denom2 = 0.0
     r1 = RAND()
     DO ig = 1, gas_ng
        iig = ig
        IF (r1>=denom2.AND.r1<denom2+PDFg(ig)) EXIT
        denom2 = denom2+PDFg(ig)
     ENDDO
     g = iig
     !IF (gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0_rknd) THEN
     !   hyparam = 2
     !ELSE
     !   hyparam = 1
     !   r1 = RAND()
     !   mu = 1.0-2.0*r1
     !   r1 = RAND()
     !   r = r1*gas_rarr(z+1)+(1.0-r1)*gas_rarr(z) !(r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*mu*r/pc_c)
     !   E0 = E0/(1.0-gas_velyes*mu*r/pc_c)
     !ENDIF
  ELSE
     prt_done = .TRUE.
     gas_numcensus(z)=gas_numcensus(z)+1
     prt_erad = prt_erad+E
  ENDIF

END SUBROUTINE diffusion1alt
