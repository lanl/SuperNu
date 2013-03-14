!Pure diffusion routine

subroutine diffusion1alt(z,g,r,mu,t,E,E0,hyparam,vacnt)

  use gasgridmod
  use particlemod
  use timestepmod
  use inputparmod
  use physconstmod
  implicit none
  !
  integer, intent(inout) :: z, g, hyparam
  real*8, intent(inout) :: r, mu, t, E, E0
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig
  real*8 :: r1, r2
  real*8 :: denom, denom2
  real*8 :: ddmct, tauA, tauL, tauR, tauS, tcensus
  real*8, dimension(gas_ng) :: PDFg

  !denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  !denom = denom+(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)
  r1 = rand()
  tauA = abs(log(r1)/(pc_c*gas_fcoef(z)*gas_sigmapg(g,z)))
  r1 = rand()
  tauS = abs(log(r1)/(pc_c*(1.0-gas_emitprobg(g,z))*(1.0-gas_fcoef(z))*gas_sigmapg(g,z)))
  r1 = rand()
  tauR = abs(log(r1)/(pc_c*gas_sigmar(g,z)))
  r1 = rand()
  tauL = abs(log(r1)/(pc_c*gas_sigmal(g,z)))
  tcensus = tsp_time+tsp_dt-t

  ddmct = min(tauA,tauS,tauR,tauL,tcensus)
  E = E*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  E0 = E0*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  t = t+ddmct
  !write(*,*) ddmct, tau, tcensus
  !if (ddmct == tau) then
     !r1 = rand()
     !PR = gas_sigmar(g,z)/denom
     !PL = gas_sigmal(g,z)/denom
     !PA = gas_fcoef(z)*gas_sigmapg(g,z)/denom
  if (ddmct == tauL) then
     if (z == 1) then
        !write(*,*) 'Non-physical left leakage'
        vacnt = .true.
        prt_done = .true.
        gas_eleft = gas_eleft+E
     !elseif (gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
     !   z = z-1
     else
     !   hyparam = 1
     !   r = gas_rarr(z)
        z = z-1
     !   r1 = rand()
     !   r2 = rand()
     !   mu = -max(r1,r2)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*r*mu/pc_c)
     !   E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
     endif
  elseif (ddmct == tauR) then
     if (z == gas_nr) then
        vacnt = .true.
        prt_done = .true.
        gas_eright = gas_eright+E
     !elseif (gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
     !   z = z+1
     else
     !   hyparam = 1
     !   r = gas_rarr(z+1)
        z = z+1
     !   r1 = rand()
     !   r2 = rand()
     !   mu = max(r1,r2)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*r*mu/pc_c)
     !   E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
     endif
  elseif (ddmct == tauA) then
     vacnt = .true.
     prt_done = .true.
     gas_edep(z) = gas_edep(z)+E
  elseif (ddmct == tauS) then
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
     !if (gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
     !   hyparam = 2
     !else
     !   hyparam = 1
     !   r1 = rand()
     !   mu = 1.0-2.0*r1
     !   r1 = rand()
     !   r = r1*gas_rarr(z+1)+(1.0-r1)*gas_rarr(z) !(r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
     !   mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
     !   E = E/(1.0-gas_velyes*mu*r/pc_c)
     !   E0 = E0/(1.0-gas_velyes*mu*r/pc_c)
     !endif
  else
     prt_done = .true.
     gas_numcensus(z)=gas_numcensus(z)+1
     gas_erad = gas_erad+E
  endif

end subroutine diffusion1alt
