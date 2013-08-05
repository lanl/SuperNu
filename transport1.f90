!Pure transport routine

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
  !
  integer, intent(inout) :: z, hyparam !,g
  integer, intent(inout) :: trndx
  real*8, intent(inout) :: r, mu, t, E, E0, wl
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig, g, binsrch
  real*8 :: r1, r2, help
  real*8 :: db, dcol, dcen, dthm, d
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, P, denom2, told, zholder
  real*8 :: bmax, x1, x2, xx0
  logical :: reflectin=.true. !temporary

  if(gas_isvelocity) then
     siglabfact = 1.0d0 - mu*r/pc_c
     dcollabfact = tsp_texp*(1d0-mu*r/pc_c)
     help = tsp_texp
  else
     siglabfact = 1d0
     dcollabfact = 1d0
     help = 1d0
  endif

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
  if(gas_isvelocity) then
     g = binsrch(wl/(1.0d0-r*mu/pc_c),gas_wl,gas_ng+1)
  else
     g = binsrch(wl,gas_wl,gas_ng+1)
  endif
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid energy bound
     if(g>gas_ng) then
        g=gas_ng
        if(gas_isvelocity) then
           wl=gas_wl(gas_ng+1)*(1.0d0-r*mu/pc_c)
        else
           wl=gas_wl(gas_ng+1)
        endif
     elseif(g<1) then
        g=1
        if(gas_isvelocity) then
           wl=gas_wl(1)*(1.0d0-r*mu/pc_c)
        else
           wl=gas_wl(1)
        endif
     else
        write(*,*) 'domain leak!!'
        prt_done = .true.
        vacnt = .true.
     endif
  endif
  ! distance to fictitious collision = dcol
  !
  if(prt_isimcanlog) then
     if(gas_cap(g,z)>0d0) then
        r1 = rand()
        dcol = abs(log(r1)/(gas_cap(g,z)*dcollabfact))
     else
        dcol = 3.0*db
     endif
  else
     if((1.0d0-gas_fcoef(z))*gas_cap(g,z)>0.0d0) then
        r1 = rand()
        dcol = abs(log(r1)/((1.0d0-gas_fcoef(z))*gas_cap(g,z)*dcollabfact))
     else
        dcol = 3.0*db
     endif
  endif
  !
  ! distance to physical collision = dthm
  if(gas_sig(z)>0.0d0) then
     r1 = rand()
     dthm = abs(log(r1)/(gas_sig(z)*dcollabfact))
  else
     dthm = 3.0*db
  endif
  ! distance to census = dcen
  dcen = abs(pc_c*(tsp_time+tsp_dt-t)/help)

  ! minimum distance = d
  d = min(dcol,dthm,db,dcen)
  
  rold = r
  r = sqrt((1.0d0-mu**2)*r**2+(d+r*mu)**2)
  told = t
  t = t + help*d/pc_c
  mu = (rold*mu+d)/r
  if(gas_isvelocity) then
     elabfact = 1.0d0 - mu*r/pc_c
  else
     elabfact = 1d0
  endif
  !calculating energy deposition and density
  !
  if(.not.prt_isimcanlog) then
     gas_edep(z)=gas_edep(z)+E*(1.0d0-exp(-gas_fcoef(z) &
          *gas_cap(g,z)*d*help))*elabfact
     !--
     if(gas_fcoef(z)*gas_cap(g,z)>0d0) then     
        gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
             (1.0d0-exp(-gas_fcoef(z)*gas_cap(g,z)*d*help))* &
             elabfact/(gas_fcoef(z)*gas_cap(g,z)*pc_c*tsp_dt)
     else
        gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
             elabfact*d*dcollabfact/(pc_c*tsp_dt)
     endif
     !--
     E = E*exp(-gas_fcoef(z)*gas_cap(g,z)*d*dcollabfact)
     if (E/E0<0.001d0) then
        vacnt = .true.
        prt_done = .true.
        gas_edep(z) = gas_edep(z) + E*elabfact
     endif
  else
     !
     gas_eraddens(g,z) = gas_eraddens(g,z)+E* &
          elabfact*d*dcollabfact/(pc_c*tsp_dt)
  endif

  ! Recalculating current group (rev. 120)
  if(gas_isvelocity) then
     g = binsrch(wl/(1.0d0-r*mu/pc_c),gas_wl,gas_ng+1)
  else
     g = binsrch(wl,gas_wl,gas_ng+1)
  endif
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid energy bound
     if(g>gas_ng) then
        g=gas_ng
        if(gas_isvelocity) then
           wl=gas_wl(gas_ng+1)*(1.0d0-r*mu/pc_c)
        else
           wl=gas_wl(gas_ng+1)
        endif
     elseif(g<1) then
        g=1
        if(gas_isvelocity) then
           wl=gas_wl(1)*(1.0d0-r*mu/pc_c)
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

  if (d == dthm) then  !physical scattering (Thomson-type)
     !
     r1 = rand()
     mu = 1.0-2.0*r1
     if(gas_isvelocity) then
        mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
        E = E*elabfact/(1.0-mu*r/pc_c)
        wl=wl*(1.0-r*mu/pc_c)/elabfact
     endif
     !
     !
  elseif (d == dcol) then  !fictitious scattering with implicit capture
     !
     r1 = rand()
     if(r1<=gas_fcoef(z).and.prt_isimcanlog) then
        vacnt=.true.
        prt_done=.true.
        gas_edep(z) = gas_edep(z) + E*elabfact
     else
        r1 = rand()
        mu = 1.0-2.0*r1
        if(gas_isvelocity) then
           mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
           E = E*elabfact/(1.0-mu*r/pc_c)
        endif
        denom2 = 0.0
        r1 = rand()
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
        !r1 = rand()
        !wl = (1d0-r1)*gas_wl(g)+r1*gas_wl(g+1)
        !wl = 0.5d0*(gas_wl(g)+gas_wl(g+1))
        !
        ! sampling sub-group Planck function:
        x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(z))
        x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(z))
        if (x2<pc_plkpk) then
           bmax = x2**3/(exp(x2)-1d0)
        elseif (x1>pc_plkpk) then
           bmax = x1**3/(exp(x1)-1d0)
        else
           bmax = pc_plkpk
        endif
        r1 = rand()
        r2 = rand()
        xx0 = (1d0-r1)*x1+r1*x2
        do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
           r1 = rand()
           r2 = rand()
           xx0 = (1d0-r1)*x1+r1*x2
        enddo
        wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
        !
        !
        if(gas_isvelocity) then
           ! converting comoving wavelength to lab frame wavelength
           wl = wl*(1.0-r*mu/pc_c)
        endif
        if (((gas_sig(z)+gas_cap(g,z))*gas_drarr(z)* &
             help >= prt_tauddmc*gas_curvcent(z)) &
             .and.(in_puretran.eqv..false.)) then
           hyparam = 2
           if(gas_isvelocity) then
              E = E*(1.0-r*mu/pc_c)
              E0 = E0*(1.0-r*mu/pc_c)
              wl = wl/(1.0-r*mu/pc_c)
           endif
        else
           hyparam = 1
        endif
     endif
     !
  elseif (d == db) then   !------boundary crossing ----
     if (mu>=0.0d0) then
        if (z == gas_nr) then
           if(.not.reflectin) then
              vacnt = .true.
              prt_done = .true.
!-- outbound luminosity tally
              gas_eright = gas_eright+E*elabfact
              gas_luminos(g) = gas_luminos(g)+mu*E/tsp_dt
           else
              mu=-mu
           endif
        ! Checking if DDMC region right
        elseif (((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
             *help >= prt_tauddmc*gas_curvcent(z+1)) &
                 .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           if(gas_isvelocity) then
              mu = (mu-r/pc_c)/(1.0-r*mu/pc_c)
           endif
           P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
           if (r1 < P) then
              hyparam = 2
              if(gas_isvelocity) then
                 E = E*elabfact
                 E0 = E0*elabfact
                 wl = wl/(1.0-r*mu/pc_c)
              endif
              z = z+1
           else
              r1 = rand()
              r2 = rand()
              mu = -max(r1,r2)
              if(gas_isvelocity) then
                 mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
              endif
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
              if (((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
                   *help >= prt_tauddmc*gas_curvcent(z+1)) &
                   .and.(in_puretran.eqv..false.)) then
                 r1 = rand()
                 if(gas_isvelocity) then
                    mu = (mu-r/pc_c)/(1.0-r*mu/pc_c)
                 endif
                 P = gas_ppl(g,z+1)*(1.0+1.5*abs(mu))
                 if (r1 < P) then
                    hyparam = 2
                    if(gas_isvelocity) then
                       E = E*elabfact
                       E0 = E0*elabfact
                       wl = wl/(1.0-r*mu/pc_c)
                    endif
                    z = z+1
                 else
                    r1 = rand()
                    r2 = rand()
                    mu = -max(r1,r2)
                    if(gas_isvelocity) then
                       mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
                    endif
                 endif
              else
                 z = z+1
              endif
           endif
        elseif (((gas_sig(z-1)+gas_cap(g,z-1))*gas_drarr(z-1) &
             *help >= prt_tauddmc*gas_curvcent(z-1)) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           if(gas_isvelocity) then
              mu = (mu-r/pc_c)/(1.0-r*mu/pc_c)
           endif
           P = gas_ppr(g,z-1)*(1.0+1.5*abs(mu))
           if (r1 < P) then
              hyparam = 2
              if(gas_isvelocity) then
                 E = E*elabfact
                 E0 = E0*elabfact
                 wl = wl/(1.0-r*mu/pc_c)
              endif
              z = z-1
           else
              r1 = rand()
              r2 = rand()
              mu = max(r1,r2)
              if(gas_isvelocity) then
                 mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
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
     gas_erad = gas_erad + E*elabfact
  endif


end subroutine transport1
