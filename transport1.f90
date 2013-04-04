!Pure transport routine

subroutine transport1(z,wl,r,mu,t,E,E0,hyparam,vacnt)

  use gasgridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  implicit none

!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event (Fleck&Cummings, 1971).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous DDMC diffusion routine through the advance.
!##################################################
  !
  integer, intent(inout) :: z, hyparam !,g
  real*8, intent(inout) :: r, mu, t, E, E0, wl
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig, g
  real*8 :: r1, r2
  real*8 :: db, dcol, dcen, d
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, P, denom2, told

  siglabfact = 1.0d0 - gas_velyes*mu*r/pc_c
  dcollabfact = gas_velno*1.0 + gas_velyes*tsp_texp*(1.0d0-mu*r/pc_c)

  ! distance to boundary = db
  if(gas_isshell) then
     if (mu < -sqrt(1.0d0-(gas_rarr(z)/r)**2)) then
        db = abs(sqrt(gas_rarr(z)**2-(1.0d0-mu**2)*r**2)+mu*r)
     else
        db = abs(sqrt(gas_rarr(z+1)**2-(1.0d0-mu**2)*r**2)-mu*r)
     endif
  else
     if (z == 1) then
        db = abs(sqrt(gas_rarr(z+1)**2-(1.0-mu**2)*r**2)-mu*r)
     elseif (mu < -sqrt(1.0d0-(gas_rarr(z)/r)**2)) then
        db = abs(sqrt(gas_rarr(z)**2-(1.0d0-mu**2)*r**2)+mu*r)
     else
        db = abs(sqrt(gas_rarr(z+1)**2-(1.0d0-mu**2)*r**2)-mu*r)
     endif
  endif

  ! Calculating current group (rev. 120)
  g = minloc(abs(gas_wl-wl/(1.0d0-gas_velyes*r*mu/pc_c)),1)
  !write(*,*) 'g 1: ',g
  if(wl/(1.0d0-gas_velyes*r*mu/pc_c)-gas_wl(g)<0d0) then
     g = g-1
  endif
  !write(*,*) 'g 2: ',g, wl/(1.0d0-gas_velyes*r*mu/pc_c)
  ! distance to collision = dcol
  if((1.0d0-gas_fcoef(z))*gas_sigmapg(g,z)>0.0d0) then
     r1 = rand()
     dcol = abs(log(r1)/((1.0d0-gas_fcoef(z))*gas_sigmapg(g,z)*dcollabfact))
  else
     dcol = 3.0*db
  endif
  ! distance to census = dcen
  dcen = abs(pc_c*(tsp_time+tsp_dt-t)/(gas_velno*1.0+gas_velyes*tsp_texp))
  ! minimum distance = d
  d = min(dcol,db,dcen)

  rold = r
  r = sqrt((1.0d0-mu**2)*r**2+(d+r*mu)**2)
  told = t
  t = t + (gas_velno*1.0+gas_velyes*tsp_texp)*d/pc_c
  mu = (rold*mu+d)/r
  elabfact = 1.0d0 - gas_velyes*mu*r/pc_c
  gas_edep(z)=gas_edep(z)+E*(1.0d0-exp(-gas_fcoef(z)*gas_sigmapg(g,z)*d))*elabfact
  E = E*exp(-gas_fcoef(z)*gas_sigmapg(g,z)*d)
  if (E/E0<0.001d0) then
     vacnt = .true.
     prt_done = .true.
     gas_edep(z) = gas_edep(z) + E*elabfact
  endif
  ! Recalculating current group (rev. 120)
  g = minloc(abs(gas_wl-wl/(1.0d0-gas_velyes*r*mu/pc_c)),1)
  if(wl/(1.0d0-gas_velyes*r*mu/pc_c)-gas_wl(g)<0d0) then
     g = g-1
  endif
  !
  if (d == dcol) then  !fictitious scattering with implicit capture
        r1 = rand()
        mu = 1.0-2.0*r1
        mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
        E = E*elabfact/(1.0-gas_velyes*mu*r/pc_c)
        denom2 = 0.0
        r1 = rand()
        do ig = 1, gas_ng
           iig = ig
           if ((r1>=denom2).and.(r1<denom2+gas_emitprobg(ig,z))) exit
           denom2 = denom2+gas_emitprobg(ig,z)
        enddo
        g = iig
        ! uniformly sampling comoving wavelength in group
        r1 = rand()
        wl = (1d0-r1)*gas_wl(g)+r1*gas_wl(g+1)
        ! converting comoving wavelength to lab frame wavelength
        wl = wl*(1.0-gas_velyes*r*mu/pc_c)
        if ((gas_sigmapg(g,z)*gas_drarr(z)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0).and.(in_puretran.eqv..false.)) then
           hyparam = 2
           E = E*(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0*(1.0-gas_velyes*r*mu/pc_c)
           wl = wl/(1.0-gas_velyes*r*mu/pc_c)
        else
           hyparam = 1
        endif
  elseif (d == db) then
     if (mu>=0.0d0) then
        if (z == gas_nr) then
           vacnt = .true.
           prt_done = .true.
           gas_eright = gas_eright+E !*elabfact
        ! Checking if DDMC region right
        elseif ((gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
                 .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
           P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
           if (r1 < P) then
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              wl = wl/(1.0-gas_velyes*r*mu/pc_c)
              z = z+1
           else
              r1 = rand()
              r2 = rand()
              mu = -max(r1,r2)
              mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           endif
        ! End of check
        else
           z = z+1
        endif
     else
        if (z==1) then
           if(gas_isshell) then
              vacnt = .true.
              prt_done = .true.
              gas_eleft = gas_eleft+E*elabfact
           else
              if ((gas_sigmapg(g,z+1)*gas_drarr(z+1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
                   .and.(in_puretran.eqv..false.)) then
                 r1 = rand()
                 mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
                 P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
                 if (r1 < P) then
                    hyparam = 2
                    E = E*elabfact
                    E0 = E0*elabfact
                    wl = wl/(1.0-gas_velyes*r*mu/pc_c)
                    z = z+1
                 else
                    r1 = rand()
                    r2 = rand()
                    mu = -max(r1,r2)
                    mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
                 endif
              else
                 z = z+1
              endif
           endif
        elseif ((gas_sigmapg(g,z-1)*gas_drarr(z-1)*(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           mu = (mu-gas_velyes*r/pc_c)/(1.0-gas_velyes*r*mu/pc_c)
           P = gas_ppr(g,z-1)*(1.0+1.5*abs(mu))
           if (r1 < P) then
              hyparam = 2
              E = E*elabfact
              E0 = E0*elabfact
              wl = wl/(1.0-gas_velyes*r*mu/pc_c)
              z = z-1
           else
              r1 = rand()
              r2 = rand()
              mu = max(r1,r2)
              mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           endif
        ! End of check
        else
           z = z-1
        endif
     endif
  elseif (d == dcen) then
     prt_done = .true.
     gas_numcensus(z) = gas_numcensus(z)+1
     gas_erad = gas_erad + E*elabfact
  endif


end subroutine transport1
