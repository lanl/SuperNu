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
  integer :: ig, iig, g, binsrch
  real*8 :: r1, r2, help, x1, x2, r3, uur, uul, uumax, r0
  real*8 :: denom, denom2, denom3, xx0, bmax
  real*8 :: ddmct, tau, tcensus, PR, PL, PA, PD
  !real*8, dimension(gas_ng) :: PDFg
  real*8 :: deleff=0.38
  real*8 :: alpeff
  !
  alpeff= 0d0 !gas_fcoef(z)**(deleff/(1-deleff))
  !
  if(gas_isvelocity) then
     help = tsp_texp
  else
     help = 1d0
  endif
  !
  ! Calculating current group (rev 206)
  g = binsrch(wl,gas_wl,gas_ng+1)
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
  denom = gas_opacleakl(g,z)+gas_opacleakr(g,z) !+gas_fcoef(z)*gas_cap(g,z)
  denom = denom+(1d0-alpeff)*(1d0-gas_emitprob(g,z))*&
       (1d0-gas_fcoef(z))*gas_cap(g,z)
!--add doppler term
  if(gas_isvelocity.and.g<gas_ng) then
     denom=denom+(gas_wl(g)/(gas_wl(g+1)-gas_wl(g)))/(pc_c*tsp_texp)
  endif
!--add analog term
  if(prt_isddmcanlog) then
     denom = denom+gas_fcoef(z)*gas_cap(g,z)
  endif
  !write(*,*) gas_emitprob(g,z),g
  r1 = rand()
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_time+tsp_dt-t
  ddmct = min(tau,tcensus)
!
!-- redshift weight
  E=E*exp(-ddmct/tsp_texp)
  E0=E0*exp(-ddmct/tsp_texp)
!--
!
  !calculating energy depostion and density
  !
  if(.not.prt_isddmcanlog) then
     gas_edep(z)=gas_edep(z)+E*(1d0-exp(-gas_fcoef(z) &
          *gas_cap(g,z)*pc_c*ddmct))
     !--
     if(gas_fcoef(z)*gas_cap(g,z)*gas_drarr(z)*help>1d-6) then
        gas_eraddens(g,z)=gas_eraddens(g,z)+E* &
             (1d0-exp(-gas_fcoef(z)*gas_cap(g,z)*pc_c*ddmct))/ &
             (gas_fcoef(z)*gas_cap(g,z)*pc_c*tsp_dt)
!        gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt        
     else
        gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt
     endif
     !--
     E=E*exp(-gas_fcoef(z)*gas_cap(g,z)*pc_c*ddmct)
     if(E/E0<0.001d0) then
        vacnt=.true.
        prt_done=.true.
        gas_edep(z)=gas_edep(z)+E
     endif
  else
     !
     gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt
  endif
  !
  t = t+ddmct
  !
  !gas_eraddens(g,z)=gas_eraddens(g,z)+E
  !
  !Recalculating histogram sum (rev. 120)
  denom = gas_opacleakl(g,z)+gas_opacleakr(g,z) !+gas_fcoef(z)*gas_cap(g,z)
  denom = denom+(1d0-alpeff)*(1d0-gas_emitprob(g,z))*&
       (1d0-gas_fcoef(z))*gas_cap(g,z)
!--add doppler term
  if(gas_isvelocity.and.g<gas_ng) then
     denom=denom+(gas_wl(g)/(gas_wl(g+1)-gas_wl(g)))/(pc_c*tsp_texp)
  endif
!--add analog term
  if(prt_isddmcanlog) then
     denom=denom+gas_fcoef(z)*gas_cap(g,z)
  endif
  if (ddmct == tau) then
     r1 = rand()
!-- right leak probability
     PR = gas_opacleakr(g,z)/denom
!-- left leak probability
     PL = gas_opacleakl(g,z)/denom
!-- absorption probability
     if(prt_isddmcanlog) then
        PA = gas_fcoef(z)*gas_cap(g,z)/denom
     else
        PA = 0d0
     endif
!-- group Doppler shift probability
     if(gas_isvelocity.and.g<gas_ng) then
        PD = (gas_wl(g)/(gas_wl(g+1)-gas_wl(g)))/(pc_c*tsp_texp*denom)
     else
        PD = 0d0
     endif
!
     !tallying radiation energy density
     !gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt
     !gas_eraddens(g,z)=gas_eraddens(g,z)+E/(denom*pc_c*tsp_dt)
     !
     !
     if (0.0d0<=r1 .and. r1<PL) then
        !gas_eraddens(g,z)=gas_eraddens(g,z)+E/(gas_opacleakl(g,z)*pc_c*tsp_dt)
        if (z == 1) then
           if(gas_isshell) then
              vacnt = .true.
              prt_done = .true.
              gas_eleft = gas_eleft+E
           else
              write(6,*) 'Non-physical left leakage', g, gas_wl(g+1), wl
           endif
        elseif ((gas_sig(z-1)+gas_cap(g,z-1))*gas_drarr(z-1) &
             *help >= prt_tauddmc*gas_curvcent(z-1)) then
           
           z = z-1
           !gas_eraddens(g,z)=gas_eraddens(g,z)+E
        else
!
!-- wl thermal resample
!            x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(z))
!            x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(z))
!            if (x2<pc_plkpk) then
!               bmax = x2**3/(exp(x2)-1d0)
!            elseif (x1>pc_plkpk) then
!               bmax = x1**3/(exp(x1)-1d0)
!            else
!               bmax = pc_plkpk
!            endif
!            r1 = rand()
!            r2 = rand()
!            xx0 = (1d0-r1)*x1+r1*x2
!            do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!               r1 = rand()
!               r2 = rand()
!               xx0 = (1d0-r1)*x1+r1*x2
!            enddo
!            wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
            r1 = rand()
            wl = 1d0/(r1/gas_wl(g+1)+(1d0-r1)/gas_wl(g))
!
!-- method changed to IMC
           hyparam = 1
!
!-- location set right bound of left cell
           r = gas_rarr(z)
!-- current particle cell set to 1 left
           z = z-1
!
!-- particl angle sampled from isotropic b.c. inward
           r1 = rand()
           r2 = rand()
           mu = -max(r1,r2)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1d0+r*mu/pc_c)
              E = E/(1.0-r*mu/pc_c)
              E0 = E0/(1.0-r*mu/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif
     elseif (PL<=r1 .and. r1<PL+PR) then
        !gas_eraddens(g,z)=gas_eraddens(g,z)+E/(gas_opacleakr(g,z)*pc_c*tsp_dt)
        if (z == gas_nr) then
           vacnt = .true.
           prt_done = .true.
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
!-- outbound luminosity tally
           if(gas_isvelocity) then
              gas_eright = gas_eright+E*(1.0+gas_rarr(gas_nr+1)*mu/pc_c)
              gas_luminos(g)=gas_luminos(g)+(E/tsp_dt)* &
                   (mu+gas_rarr(gas_nr+1)/pc_c)/ &
                   (1.0+gas_rarr(gas_nr+1)*mu/pc_c)
           else
              gas_eright = gas_eright+E
              gas_luminos(g)=gas_luminos(g)+(E/tsp_dt)*mu
           endif
        elseif ((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
             *help >= prt_tauddmc*gas_curvcent(z+1)) then
           !gas_eraddens(g,z)=gas_eraddens(g,z)+E
           z = z+1
           !gas_eraddens(g,z)=gas_eraddens(g,z)+E
        else
!
!-- wavelength thermal resample
!            x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(z))
!            x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(z))
!            if (x2<pc_plkpk) then
!               bmax = x2**3/(exp(x2)-1d0)
!            elseif (x1>pc_plkpk) then
!               bmax = x1**3/(exp(x1)-1d0)
!            else
!               bmax = pc_plkpk
!            endif
!            r1 = rand()
!            r2 = rand()
!            xx0 = (1d0-r1)*x1+r1*x2
!            do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!               r1 = rand()
!               r2 = rand()
!               xx0 = (1d0-r1)*x1+r1*x2
!            enddo
!            wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
            r1 = rand()
            wl = 1d0/(r1/gas_wl(g+1)+(1d0-r1)/gas_wl(g))
!
!-- method changed to IMC
           hyparam = 1
!
!-- location set left bound of right cell
           r = gas_rarr(z+1)
!-- current particle cell set 1 right
           z = z+1
!
!--  particl angle sampled from isotropic b.c. outward
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
              E = E/(1.0-r*mu/pc_c)
              E0 = E0/(1.0-r*mu/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif
     elseif (PL+PR<=r1 .and. r1<PL+PR+PA) then
        vacnt = .true.
        prt_done = .true.
        gas_edep(z) = gas_edep(z)+E
     elseif(PL+PR+PA<=r1.and.r1<PL+PR+PA+PD) then
!-- group shift
        if(g<gas_ng) then
           g = g+1
           !r1 = rand()
           !wl = 1d0/((1d0-r1)/gas_wl(g)+r1/gas_wl(g+1))
           wl = gas_wl(g)
        else
           g = gas_ng
           wl = gas_wl(gas_ng+1)
        endif
        if((gas_sig(z)+gas_cap(g,z))*gas_drarr(z) &
             *help >= prt_tauddmc*gas_curvcent(z)) then
           hyparam = 2
        else
           hyparam = 1
!-- direction sampled isotropically           
           r1 = rand()
           mu = 1.0-2.0*r1
!-- position sampled uniformly
           r1 = rand()
           r = (r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
              E = E/(1.0-r*mu/pc_c)
              E0 = E0/(1.0-r*mu/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif
!
     else
        denom2 = 0d0
        do ig = 1, gas_ng
           if(ig.ne.g) then
              denom2 = denom2+gas_emitprob(ig,z)
           endif
        enddo
        denom3 = 0d0
        r1 = rand()
        do ig = 1, gas_ng
           if(ig.ne.g) then
              iig = ig
              if((r1>=denom3).and.(r1<denom3+gas_emitprob(ig,z)/denom2)) exit
              denom3 = denom3+gas_emitprob(ig,z)/denom2
           endif
        enddo
        !write(*,*) 'Scatter: ',g,'to ',iig
        g = iig
        r1 = rand()
        wl = 1d0/((1d0-r1)/gas_wl(g)+r1/gas_wl(g+1))
        ! during DDMC phase, wavelength is only a placeholder for group 
        !wl = 0.5d0*(gas_wl(g)+gas_wl(g+1))
        !
        if ((gas_sig(z)+gas_cap(g,z))*gas_drarr(z) &
             *help >= prt_tauddmc*gas_curvcent(z)) then
           hyparam = 2
        else
           hyparam = 1
!-- wavelength thermal resample
!            x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(z))
!            x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(z))
!            if (x2<pc_plkpk) then
!               bmax = x2**3/(exp(x2)-1d0)
!            elseif (x1>pc_plkpk) then
!               bmax = x1**3/(exp(x1)-1d0)
!            else
!               bmax = pc_plkpk
!            endif
!            r1 = rand()
!            r2 = rand()
!            xx0 = (1d0-r1)*x1+r1*x2
!            do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!               r1 = rand()
!               r2 = rand()
!               xx0 = (1d0-r1)*x1+r1*x2
!            enddo
!            wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
!
!-- direction sampled isotropically           
           r1 = rand()
           mu = 1.0-2.0*r1
!-- position sampled uniformly
            r1 = rand()
            r = (r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
!-- position sampled from source tilt
!            r1 = 0d0
!            r2 = 1d0
!            uul = gas_tempb(z)**4
!            uur = gas_tempb(z+1)**4
!            uumax = max(uul,uur)
!            do while (r2 > r1)
!               r3 = rand()
!               r0 = (r3*gas_rarr(z+1)**3+(1.0-r3)*gas_rarr(z)**3)**(1.0/3.0)
!               r3 = (r0-gas_rarr(z))/gas_drarr(z)
!               r1 = (r3*uur+(1d0-r3)*uul)/uumax
!               r2 = rand()
!            enddo
!            r = r0
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
              E = E/(1.0-mu*r/pc_c)
              E0 = E0/(1.0-mu*r/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif
     !else
     !   stop 'diffusion1: invalid histogram sample'
     endif
  else
     !gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt
     prt_done = .true.
     gas_numcensus(z)=gas_numcensus(z)+1
     gas_erad = gas_erad+E
  endif


end subroutine diffusion1
