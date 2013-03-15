!Pure diffusion routine

subroutine diffusion1(z,g,r,mu,t,E,E0,hyparam,vacnt)

  use gasgridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  implicit none

!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  !
  integer, intent(inout) :: z, g, hyparam
  real*8, intent(inout) :: r, mu, t, E, E0
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig
  real*8 :: r1, r2
  real*8 :: denom, denom2
  real*8 :: ddmct, tau, tcensus, PR, PL, PA
  real*8, dimension(gas_ng) :: PDFg

  denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  denom = denom+(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)
  r1 = rand()
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_time+tsp_dt-t
  ddmct = min(tau,tcensus)
  E = E*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  E0 = E0*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  t = t+ddmct
  !write(6,*) ddmct, tau, tcensus
  !gas_edep(z) = gas_edep(z)+E
  if (ddmct == tau) then
     r1 = rand()
     PR = gas_sigmar(g,z)/denom
     PL = gas_sigmal(g,z)/denom
     PA = gas_fcoef(z)*gas_sigmapg(g,z)/denom
     if (0.0d0<=r1 .and. r1<PL) then
        if (z == 1) then
           write(6,*) 'Non-physical left leakage'
           !vacnt = .true.
           !prt_done = .true.
           !gas_eleft = gas_eleft+E
        elseif (gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
           z = z-1
        else
           hyparam = 1
           r = gas_rarr(z)
           z = z-1
           r1 = rand()
           r2 = rand()
           mu = -max(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
        endif
     elseif (PL<=r1 .and. r1<PL+PR) then
        if (z == gas_nr) then
           vacnt = .true.
           prt_done = .true.
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
           gas_eright = gas_eright+E*(1.0+gas_velyes*gas_rarr(gas_nr+1)*mu/pc_c)
        elseif (gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
           z = z+1
        else
           hyparam = 1
           r = gas_rarr(z+1)
           z = z+1
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
        endif
     elseif (PL+PR<=r1 .and. r1<PL+PR+PA) then
        vacnt = .true.
        prt_done = .true.
        gas_edep(z) = gas_edep(z)+E
     else
        denom2 = gas_sigmap(z)-gas_ppick(g)*gas_sigmapg(g,z)
        do ig = 1, gas_ng
           PDFg(ig) = gas_emitprobg(ig,z)*gas_sigmap(z)/denom2
        enddo
        PDFg(g)=0.0
        denom2 = 0.0
        r1 = rand()
        do ig = 1, gas_ng
           iig = ig
           if (r1>=denom2.and.r1<denom2+PDFg(ig)) EXIT
           denom2 = denom2+PDFg(ig)
        enddo
        g = iig
        if (gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
           hyparam = 2
        else
           hyparam = 1
           r1 = rand()
           mu = 1.0-2.0*r1
           r1 = rand()
           r = (r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           E = E/(1.0-gas_velyes*mu*r/pc_c)
           E0 = E0/(1.0-gas_velyes*mu*r/pc_c)
        endif
     endif
  else
     prt_done = .true.
     gas_numcensus(z)=gas_numcensus(z)+1
     gas_erad = gas_erad+E
  endif

end subroutine diffusion1
