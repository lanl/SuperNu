subroutine transport1(z,wl,r,mu,t,E,E0,hyparam,vacnt,trndx)

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
  real*8,parameter :: cinv = 1d0/pc_c
  !
  integer, intent(inout) :: z, hyparam !,g
  integer, intent(in) :: trndx
  real*8, intent(inout) :: r, mu, t, E, E0, wl
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig, g, binsrch
  real*8 :: r1, r2, thelp,thelpinv, mu0, wl0, vdiff, rdop1, rdop2
  real*8 :: db, dcol, dcen, dthm, ddop, d
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, P, denom2, told, zholder, muold
  real*8 :: bmax, x1, x2, xx0, ddop1, ddop2
  real*8 :: dtinv
  real*8 :: help

!--------------------------------------------------------------
!
!-- shortcut
  dtinv = 1d0/tsp_dt

  if(gas_isvelocity) then
     siglabfact = 1.0d0 - mu*r*cinv
     dcollabfact = tsp_t*(1d0-mu*r*cinv)
     thelp = tsp_t
  else
     siglabfact = 1d0
     dcollabfact = 1d0
     thelp = 1d0
  endif
  thelpinv = 1d0/thelp

!
!-- calculating current group (rev. 120)
  if(gas_isvelocity) then
     g = binsrch(wl/(1.0d0-r*mu*cinv),gas_wl,gas_ng+1,in_ng)
  else
     g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
  endif
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid energy bound
     if(g>gas_ng) then
        g=gas_ng
        if(gas_isvelocity) then
           wl=gas_wl(gas_ng+1)*(1.0d0-r*mu*cinv)
        else
           wl=gas_wl(gas_ng+1)
        endif
     elseif(g<1) then
        g=1
        if(gas_isvelocity) then
           wl=gas_wl(1)*(1.0d0-r*mu*cinv)
        else
           wl=gas_wl(1)
        endif
     else
        write(*,*) 'domain leak!!'
        prt_done = .true.
        vacnt = .true.
     endif
  endif
!
!== DISTANCE CALCULATIONS
!
!-- distance to boundary = db
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
!
!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(gas_cap(g,z)>0d0) then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        dcol = abs(log(r1)/(gas_cap(g,z)*dcollabfact))
     else
        dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
     endif
  else
     if((1.0d0-gas_fcoef(z))*gas_cap(g,z)>0.0d0) then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        dcol = abs(log(r1)/((1.0d0-gas_fcoef(z))*gas_cap(g,z)*dcollabfact))
     else
        dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
     endif
  endif
!
!-- distance to Thomson-type collision = dthm
  if(gas_sig(z)>0.0d0) then
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     dthm = abs(log(r1)/(gas_sig(z)*dcollabfact))
  else
     dthm = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- distance to census = dcen
  dcen = abs(pc_c*(tsp_t+tsp_dt-t)*thelpinv)
!
!-- distance to Doppler shift = ddop
   if(gas_isvelocity.and.g<gas_ng) then
!      r1 = rand()
!      prt_tlyrand=prt_tlyrand+1
!      ddop = pc_c*tsp_t*(gas_wl(g+1)-gas_wl(g))*abs(log(r1))/(gas_wl(g)*dcollabfact)
!     wl = r1*gas_wl(g)+(1d0-r1)*gas_wl(g+1) !uniform sample
!      wl=1d0/(r1/gas_wl(g+1) + (1d0-r1)/gas_wl(g))  !reciprocal sample
!      wl=wl*(1d0-mu*r*cinv)
!      ddop = pc_c*(1d0-mu*r*cinv)*(1d0-wl/(1d0-mu*r*cinv*gas_wl(g+1)))
!     ddop = pc_c*(1d0-mu*r*cinv)*(1d0-&
!          gas_wl(g)*log(gas_wl(g+1)/gas_wl(g))/(gas_wl(g+1)-gas_wl(g)))
!     write(*,*) pc_c*(wl/gas_wl(g+1)-1d0)+r*mu
      ddop = abs(pc_c*(1d0-wl/gas_wl(g+1))-r*mu)
  else
     ddop = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- minimum distance = d
!  if(tsp_it==29) write(*,*) dcol,dthm,db,dcen,ddop
  d = min(dcol,dthm,db,dcen,ddop)
!
!== END OF DISTANCE CALCULATIONS
!
!-- position, angle, time update  
  rold = r
  r = sqrt((1.0d0-mu**2)*r**2+(d+r*mu)**2)
!  r = sqrt(r**2+d**2+2d0*d*r*mu)
  told = t
  t = t + thelp*d*cinv
  muold = mu
  mu = (rold*mu+d)/r

!-- transformation factor set
  if(gas_isvelocity) then
     elabfact = 1.0d0 - muold*rold*cinv
  else
     elabfact = 1d0
  endif
  !calculating energy deposition and density
  !
  if(.not.prt_isimcanlog) then
        gas_edep(z)=gas_edep(z)+E*(1.0d0-exp(-gas_fcoef(z) &
             *gas_cap(g,z)*siglabfact*d*thelp))*elabfact
     !--
     if(gas_fcoef(z)*gas_cap(g,z)*gas_drarr(z)*thelp>1d-6) then     
        gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
             (1.0d0-exp(-gas_fcoef(z)*siglabfact*gas_cap(g,z)*d*thelp))* &
             elabfact/(gas_fcoef(z)*siglabfact*gas_cap(g,z)*pc_c*tsp_dt)
     else
        gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
             elabfact*d*dcollabfact*cinv*dtinv
     endif
     !--
!     E = E*exp(-gas_fcoef(z)*gas_cap(g,z)*d*dcollabfact)
     E = E*exp(-gas_fcoef(z)*gas_cap(g,z)*siglabfact*d*thelp)

  else
     !
     gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
          elabfact*d*dcollabfact*cinv*dtinv
  endif

!-- transformation factor reset
  if(gas_isvelocity) then
     elabfact = 1.0d0 - mu*r*cinv
  else
     elabfact = 1d0
  endif

  !
  if(d == ddop) then !group shift
!     r1 = rand()
!     prt_tlyrand=prt_tlyrand+1
!-- redshifting
     if(g<gas_ng) then
        g = g+1
!-- lab frame wavelength
!     wl = r1*gas_wl(g)+(1d0-r1)*gas_wl(g+1) !uniform sample
!        wl=1d0/(r1/gas_wl(g+1) + (1d0-r1)/gas_wl(g))  !reciprocal sample
!        wl = wl*(1d0-mu*r*cinv)
        wl = (gas_wl(g)+1d-6*(gas_wl(g+1)-gas_wl(g)))*(1d0-mu*r*cinv)
     else
        r1 = rand()
        prt_tlyrand=prt_tlyrand+1
!     wl = r1*gas_wl(gas_ng)+(1d0-r1)*gas_wl(gas_ng+1) !uniform sample
        wl=1d0/(r1/gas_wl(g+1) + (1d0-r1)/gas_wl(gas_ng))  !reciprocal sample
        wl = wl*(1d0-mu*r*cinv)
!        wl = gas_wl(gas_ng+1)*(1d0-mu*r*cinv)
     endif
!-- check if ddmc region
     if (((gas_sig(z)+gas_cap(g,z))*gas_drarr(z)* &
          thelp >= prt_tauddmc*gas_curvcent(z)) &
          .and.(in_puretran.eqv..false.)) then
        hyparam = 2
        gas_methodswap(z)=gas_methodswap(z)+1
        if(gas_isvelocity) then
!-- velocity effects accounting
           gas_evelo=gas_evelo+E*r*mu*cinv
!
           E = E*(1.0-r*mu*cinv)
           E0 = E0*(1.0-r*mu*cinv)
           wl = wl/(1.0-r*mu*cinv)
        endif
     else
        hyparam = 1
     endif
!
  elseif (d == dthm) then  !physical scattering (Thomson-type)
     !
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu = 1.0-2.0*r1
     if(abs(mu)<0.0000001d0) then
        mu = 0.0000001d0
     endif
     if(gas_isvelocity) then
        mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
        help = 1d0/(1.0-mu*r*cinv)
        gas_evelo=gas_evelo+E*(1d0-elabfact*help)
!
        E = E*elabfact*help
!        E0 = E0*elabfact/(1.0-mu*r*cinv)
        wl = wl*(1.0-mu*r*cinv)/elabfact
     endif
     !
     !
  elseif (d == dcol) then  !fictitious scattering with implicit capture
     !
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     if(r1<=gas_fcoef(z).and.prt_isimcanlog) then
        vacnt=.true.
        prt_done=.true.
        gas_edep(z) = gas_edep(z) + E*elabfact
!-- velocity effects accounting
        gas_evelo = gas_evelo+E*(1d0-elabfact)
!
     else
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        mu = 1.0-2.0*r1
        if(abs(mu)<0.0000001d0) then
           mu = 0.0000001d0
        endif
        if(gas_isvelocity) then
           mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1.0-mu*r*cinv)
           gas_evelo=gas_evelo+E*(1d0-elabfact*help)
!
           E = E*elabfact*help
!           wl = wl*(1.0-mu*r*cinv)/elabfact
           
        endif
!
        denom2 = 0.0
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        do ig = 1, gas_ng
           iig = ig
           if ((r1>=denom2).and.(r1<denom2+gas_emitprob(ig,z))) exit
           denom2 = denom2+gas_emitprob(ig,z)
        enddo
        g = iig
        !(rev 121): calculating radiation energy tally per group
        !gas_eraddens(g,z)=gas_eraddens(g,z)+E*elabfact
        !-------------------------------------------------------
        ! sampling comoving wavelength in group
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        wl = 1d0/((1d0-r1)/gas_wl(g)+r1/gas_wl(g+1))
        !wl = (1d0-r1)*gas_wl(g)+r1*gas_wl(g+1)
        !wl = 0.5d0*(gas_wl(g)+gas_wl(g+1))
        !
        ! sampling sub-group Planck function:
!         x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(z))
!         x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(z))
!         if (x2<pc_plkpk) then
!            bmax = x2**3/(exp(x2)-1d0)
!         elseif (x1>pc_plkpk) then
!            bmax = x1**3/(exp(x1)-1d0)
!         else
!            bmax = pc_plkpk
!         endif
!         r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!         r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!         xx0 = (1d0-r1)*x1+r1*x2
!         do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!            r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            xx0 = (1d0-r1)*x1+r1*x2
!         enddo
!         wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
        !
        !
        if(gas_isvelocity) then
!-- converting comoving wavelength to lab frame wavelength
           wl = wl*(1.0-r*mu*cinv)
        endif
        if (((gas_sig(z)+gas_cap(g,z))*gas_drarr(z)* &
             thelp >= prt_tauddmc*gas_curvcent(z)) &
             .and.(in_puretran.eqv..false.)) then
           hyparam = 2
           gas_methodswap(z)=gas_methodswap(z)+1
           if(gas_isvelocity) then
!-- velocity effects accounting
              gas_evelo = gas_evelo+E*r*mu*cinv
!
              E = E*(1.0-r*mu*cinv)
              E0 = E0*(1.0-r*mu*cinv)
              wl = wl/(1.0-r*mu*cinv)
           endif
        else
           hyparam = 1
        endif
     endif
     !
  elseif (d == db) then   !------boundary crossing ----
     if (mu>=0.0d0) then
        if (z == gas_nr) then
!           if(g/=1) then
              vacnt = .true.
              prt_done = .true.
!
!-- retrieve lab frame group
              g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
!
!-- check group bounds
              if(g>gas_ng.or.g<1) then
                 if(g>gas_ng) then
                    g=gas_ng
                    wl=gas_wl(gas_ng+1)
                 else
                    g=1
                    wl=gas_wl(1)
                 endif
              endif
!
!-- outbound luminosity tally
!-- velocity effects accounting
              gas_evelo = gas_evelo+E*(1d0-elabfact)
!
              gas_eright = gas_eright+E*elabfact
              gas_luminos(g) = gas_luminos(g)+E*dtinv
              gas_lumdev(g) = gas_lumdev(g)+(E*dtinv)**2
              gas_lumnum(g) = gas_lumnum(g)+1
!              gas_luminos(g) = gas_luminos(g)+mu*E*dtinv
!            else
!               r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!               r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!               mu = -max(r1,r2)
!               mu = -mu
!               mu = -r1
!            endif
        ! Checking if DDMC region right
        elseif (((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
             *thelp >= prt_tauddmc*gas_curvcent(z+1)) &
                 .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(gas_isvelocity) then
              mu = (mu-r*cinv)/(1.0-r*mu*cinv)
           endif
           P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
!-- new albedo test
!            P = P+gas_ppl(g,z+1)* &
!                 (2d0/((gas_sig(z+1)+gas_cap(g,z+1))*gas_rarr(z+1)*thelp))*&
!                 (7.09752*abs(mu)**4-15.5997*abs(mu)**3+9.62187*abs(mu)**2&
!                 -3.9273*abs(mu)+1.48571)
!            if(P<0d0) then
!               P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
!            endif
!--
           if (r1 < P) then
              hyparam = 2
              gas_methodswap(z)=gas_methodswap(z)+1
              if(gas_isvelocity) then
!-- velocity effects accounting
                 gas_evelo=gas_evelo+E*(1d0-elabfact)
!
                 E = E*elabfact
                 E0 = E0*elabfact
                 wl = wl/elabfact
              endif
              z = z+1
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = -max(r1,r2)
              if(gas_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
           endif
        ! End of check
        else
           z = z+1
           r = gas_rarr(z)
        endif
     else
        if (z==1) then
           if(gas_isshell) then
              vacnt = .true.
              prt_done = .true.
!-- velocity effects accounting
              gas_evelo = gas_evelo+E*(1d0-elabfact)
!
              gas_eleft = gas_eleft+E*elabfact
           else
              if (((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
                   *thelp >= prt_tauddmc*gas_curvcent(z+1)) &
                   .and.(in_puretran.eqv..false.)) then
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 if(gas_isvelocity) then
                    mu = (mu-r*cinv)/(1.0-r*mu*cinv)
                 endif
                 P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
!-- new albedo test
!                  P = P+gas_ppl(g,z+1)* &
!                       (2d0/((gas_sig(z+1)+gas_cap(g,z+1))*gas_rarr(z+1)*thelp))*&
!                       (7.09752*abs(mu)**4-15.5997*abs(mu)**3+9.62187*abs(mu)**2-&
!                       3.9273*abs(mu)+1.48571)
!                  if(P<0d0) then
!                     P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
!                  endif
!--
                 if (r1 < P) then
                    hyparam = 2
                    gas_methodswap(z)=gas_methodswap(z)+1
                    if(gas_isvelocity) then
!-- velocity effects accounting
                       gas_evelo=gas_evelo+E*(1d0-elabfact)
!
                       E = E*elabfact
                       E0 = E0*elabfact
                       wl = wl/elabfact
                    endif
                    z = z+1
                 else
                    r1 = rand()
                    prt_tlyrand = prt_tlyrand+1
                    r2 = rand()
                    prt_tlyrand = prt_tlyrand+1
                    mu = -max(r1,r2)
                    if(gas_isvelocity) then
                       mu = (mu+r*cinv)/(1.0+r*mu*cinv)
                    endif
                 endif
              else
                 z = z+1
              endif
           endif
        elseif (((gas_sig(z-1)+gas_cap(g,z-1))*gas_drarr(z-1) &
             *thelp >= prt_tauddmc*gas_curvcent(z-1)) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(gas_isvelocity) then
              mu = (mu-r*cinv)/(1.0-r*mu*cinv)
!-- amplification
!
!              E0=E0*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
!              E = E*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
              if(mu<0d0) then
!-- velocity effects accounting
                 help = 1d0/abs(mu)
                 gas_evelo = gas_evelo-E*2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv
!
                 E0 = E0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
                 E = E*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
              endif
               
!--
           endif
           P = gas_ppr(g,z-1)*(1.0+1.5*abs(mu))
!-- new albedo test
!            P = P+gas_ppr(g,z-1)* &
!                 (2d0/((gas_sig(z-1)+gas_cap(g,z-1))*gas_rarr(z)*thelp))* &
!                 (7.09752*abs(mu)**4-15.5997*abs(mu)**3+9.62187*abs(mu)**2-&
!                 3.9273*abs(mu)+1.48571)
!            if(P<0d0) then
!               P = gas_ppr(g,z-1)*(1.0+1.5*abs(mu))
!            endif
!--
           if (r1 < P) then
              hyparam = 2
              gas_methodswap(z)=gas_methodswap(z)+1
              if(gas_isvelocity) then
!-- velocity effects accounting
                 gas_evelo = gas_evelo+E*(1d0-elabfact)
!
                 E = E*elabfact
                 E0 = E0*elabfact
                 wl = wl/elabfact
              endif
              z = z-1
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = max(r1,r2)
              if(gas_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
           endif
        ! End of check
        else
           z = z-1
        endif
     endif
  elseif (d == dcen) then
     prt_done = .true.
     gas_numcensus(z) = gas_numcensus(z)+1
!     gas_erad = gas_erad + E*elabfact
!
  endif

  if (E<1d-6*E0.and..not.vacnt) then
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     if(r1<0.5d0) then
        vacnt = .true.
        prt_done = .true.
        gas_edep(z) = gas_edep(z) + E*elabfact
!-- velocity effects accounting
        gas_evelo=gas_evelo+E*(1d0-elabfact)
!
     else
!-- weight addition accounted for in external source
        gas_eext=gas_eext+E
!
        E = 2d0*E
        E0 = 2d0*E0
     endif
  endif

end subroutine transport1
