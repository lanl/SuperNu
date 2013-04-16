!Pure diffusion routine

subroutine diffusion1(z,wl,r,mu,t,E,E0,hyparam,vacnt)

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
  integer, intent(inout) :: z, hyparam !,g
  real*8, intent(inout) :: r, mu, t, E, E0, wl
  logical, intent(inout) :: vacnt
  !
  integer :: ig, iig, g
  real*8 :: r1, r2
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus, PR, PL, PA
  real*8, dimension(gas_ng) :: PDFg

  ! Calculating current group (rev. 120)
  g = minloc(abs(gas_wl-wl),1)
  if(wl-gas_wl(g)<0d0) then
     g = g-1
  endif
  !
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid bound
     if(g>gas_ng) then
        g=gas_ng
        wl=gas_wl(gas_ng+1)
     elseif(g<1) then
        g=1
        wl=gas_wl(1)
     else
        write(*,*) 'domain leak!!'
        prt_done = .true.
        vacnt = .true.
     endif
  endif
  !
  denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  denom = denom+(1d0-gas_emitprobg(g,z))*(1d0-gas_fcoef(z))*gas_sigmapg(g,z)
  !write(*,*) gas_emitprobg(g,z),g
  r1 = rand()
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_time+tsp_dt-t
  ddmct = min(tau,tcensus)
  E = E*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  E0 = E0*(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  ! Recalculating comoving wavelength (rev. 120)
  wl = wl/(gas_velno*1.0+gas_velyes*exp(-ddmct/tsp_texp))
  ! Recalculating current group (rev. 120)
  g = minloc(abs(gas_wl-wl),1)
  if(wl-gas_wl(g)<0d0) then
     g = g-1
  endif
  !
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid bound
     if(g>gas_ng) then
        g=gas_ng
        wl=gas_wl(gas_ng+1)
     elseif(g<1) then
        g=1
        wl=gas_wl(1)
     else
        write(*,*) 'domain leak!!'
        prt_done = .true.
        vacnt = .true.
     endif
  endif
  !
  t = t+ddmct
  !
  !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E
  !
  !Recalculating histogram sum (rev. 120)
  denom = gas_sigmal(g,z)+gas_sigmar(g,z)+gas_fcoef(z)*gas_sigmapg(g,z)
  denom = denom+(1d0-gas_emitprobg(g,z))*(1d0-gas_fcoef(z))*gas_sigmapg(g,z)
  
  if (ddmct == tau) then
     r1 = rand()
     PR = gas_sigmar(g,z)/denom
     PL = gas_sigmal(g,z)/denom
     PA = gas_fcoef(z)*gas_sigmapg(g,z)/denom
     !tallying radiation energy density
     !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E/(tsp_dt*pc_c*denom)
     gas_eraddensg(g,z)=gas_eraddensg(g,z)+E*ddmct/tsp_dt
     !
     !
     if (0.0d0<=r1 .and. r1<PL) then
        if (z == 1) then
           if(gas_isshell) then
              vacnt = .true.
              prt_done = .true.
              gas_eleft = gas_eleft+E
           else
              write(6,*) 'Non-physical left leakage', g, gas_wl(g+1), wl
           endif
        elseif ((gas_sig(z-1)+gas_sigmapg(g,z-1))*gas_drarr(z-1) &
             *(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
           z = z-1
           !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E
        else
           hyparam = 1
           r = gas_rarr(z)
           z = z-1
           !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E
           r1 = rand()
           r2 = rand()
           mu = -max(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+gas_velyes*r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
           wl = wl*(1.0-gas_velyes*r*mu/pc_c)
        endif
     elseif (PL<=r1 .and. r1<PL+PR) then
        if (z == gas_nr) then
           vacnt = .true.
           prt_done = .true.
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
           gas_eright = gas_eright+E*(1.0+gas_velyes*gas_rarr(gas_nr+1)*mu/pc_c)
        elseif ((gas_sig(z+1)+gas_sigmapg(g,z+1))*gas_drarr(z+1) &
             *(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
           z = z+1
           !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E
        else
           hyparam = 1
           r = gas_rarr(z+1)
           z = z+1
           !gas_eraddensg(g,z)=gas_eraddensg(g,z)+E
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
           mu = (mu+gas_velyes*r/pc_c)/(1.0+r*mu/pc_c)
           E = E/(1.0-gas_velyes*r*mu/pc_c)
           E0 = E0/(1.0-gas_velyes*r*mu/pc_c)
           wl = wl*(1.0-gas_velyes*r*mu/pc_c)
        endif
     elseif (PL+PR<=r1 .and. r1<PL+PR+PA) then
        vacnt = .true.
        prt_done = .true.
        !gas_edep(z) = gas_edep(z)+E
        !if(tsp_tn==17) then
        !   write(*,*) PL,PR,PA, r1
        !   write(*,*) 'here: ', E, hyparam
        !endif        
     else
        !write(*,*) 'scatter'
        denom2 = 0d0
        do ig = 1, gas_ng
           if(ig.ne.g) then
              denom2 = denom2+gas_emitprobg(ig,z)
           endif
        enddo
        denom3 = 0d0
        r1 = rand()
        do ig = 1, gas_ng
           if(ig.ne.g) then
              iig = ig
              if((r1>=denom3).and.(r1<denom3+gas_emitprobg(ig,z)/denom2)) exit
              denom3 = denom3+gas_emitprobg(ig,z)/denom2
           endif
        enddo
        g = iig
        r1 = rand()
        wl = (1d0-r1)*gas_wl(g)+r1*gas_wl(g+1)
        !
        if ((gas_sig(z)+gas_sigmapg(g,z))*gas_drarr(z) &
             *(gas_velno*1d0+gas_velyes*tsp_texp)>=5.0d0) then
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
           wl = wl*(1.0-gas_velyes*r*mu/pc_c)
        endif
     !else
     !   stop 'diffusion1: invalid histogram sample'
     endif
  else
     gas_eraddensg(g,z)=gas_eraddensg(g,z)+E*ddmct/tsp_dt
     prt_done = .true.
     gas_numcensus(z)=gas_numcensus(z)+1
     gas_erad = gas_erad+E
  endif

end subroutine diffusion1
