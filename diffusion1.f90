!Pure diffusion routine

SUBROUTINE diffusion1(z,g,r,mu,t,E,E0,hyparam,vacnt)

  USE gasgridmod
  USE timestepmod
  USE physconstmod
  USE particlemod
  USE inputparmod
  IMPLICIT NONE

!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  !
  INTEGER, INTENT(INOUT) :: z, g, hyparam
  REAL*8, INTENT(INOUT) :: r, mu, t, E, E0
  LOGICAL, INTENT(INOUT) :: vacnt
  !
  INTEGER :: ig, iig
  REAL*8 :: r1, r2
  REAL*8 :: denom, denom2
  REAL*8 :: ddmct, tau, tcensus, PR, PL, PA
  REAL*8, DIMENSION(gas_ng) :: PDFg

  denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  denom = denom+(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)
  r1 = RAND()
  tau = ABS(LOG(r1)/(pc_c*denom))
  tcensus = tsp_time+tsp_dt-t
  ddmct = MIN(tau,tcensus)
  E = E*(gas_velno*1.0+gas_velyes*EXP(-ddmct/tsp_texp))
  E0 = E0*(gas_velno*1.0+gas_velyes*EXP(-ddmct/tsp_texp))
  t = t+ddmct
  !WRITE(*,*) ddmct, tau, tcensus
  !gas_edep(z) = gas_edep(z)+E
  IF (ddmct == tau) THEN
     r1 = RAND()
     PR = gas_sigmar(g,z)/denom
     PL = gas_sigmal(g,z)/denom
     PA = gas_fcoef(z)*gas_sigmapg(g,z)/denom
     IF (0.0d0<=r1 .AND. r1<PL) THEN
        IF (z == 1) THEN
           WRITE(*,*) 'Non-physical left leakage'
           !vacnt = .TRUE.
           !prt_done = .TRUE.
           !gas_eleft = gas_eleft+E
        ELSEIF (gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) THEN
           z = z-1
        ELSE
           hyparam = 1
           r = gas_rarr(z)
           z = z-1
           r1 = RAND()
           r2 = RAND()
           mu = -MAX(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
        ENDIF
     ELSEIF (PL<=r1 .AND. r1<PL+PR) THEN
        IF (z == gas_nr) THEN
           vacnt = .TRUE.
           prt_done = .TRUE.
           r1 = RAND()
           r2 = RAND()
           mu = MAX(r1,r2)
           gas_eright = gas_eright+E*(1.0+gas_velyes*gas_rarr(gas_nr+1)*mu/pc_c)
        ELSEIF (gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) THEN
           z = z+1
        ELSE
           hyparam = 1
           r = gas_rarr(z+1)
           z = z+1
           r1 = RAND()
           r2 = RAND()
           mu = MAX(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
        ENDIF
     ELSEIF (PL+PR<=r1 .AND. r1<PL+PR+PA) THEN
        vacnt = .TRUE.
        prt_done = .TRUE.
        gas_edep(z) = gas_edep(z)+E
     ELSE
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
        IF (gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) THEN
           hyparam = 2
        ELSE
           hyparam = 1
           r1 = RAND()
           mu = 1.0-2.0*r1
           r1 = RAND()
           r = (r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           E = E/(1.0-gas_velyes*mu*r/pc_c)
           E0 = E0/(1.0-gas_velyes*mu*r/pc_c)
        ENDIF
     ENDIF
  ELSE
     prt_done = .TRUE.
     gas_numcensus(z)=gas_numcensus(z)+1
     gas_erad = gas_erad+E
  ENDIF

END SUBROUTINE diffusion1
