!Pure diffusion routine

subroutine diffusion1(z,wl,r,mu,t,E,E0,hyparam,vacnt,partnum)

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
  integer, intent(in) :: partnum !ipart
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
!-- lumped quantities -----------------------------------------
  integer :: gminlump, gmaxlump
  real*8 :: emitlump, speclump
  real*8 :: opacleakllump, opacleakrlump, caplump
  real*8 :: specig

  integer :: glump, gunlump
  integer :: glumps(gas_ng)

!--------------------------------------------------------------
  !
!  alpeff= gas_fcoef(z)**(deleff/(1-deleff))
  alpeff = 0d0
  !
  if(gas_isvelocity) then
     help = tsp_t
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
!-- lump testing ---------------------------------------------
!
  glump = 0
  gunlump = gas_ng
  glumps = 0
!-- initializing lumped quantities
  speclump = 0d0
  emitlump = 0d0
  opacleakllump = 0d0
  opacleakrlump = 0d0
  caplump = 0d0
!
!-- find lumpable groups
  do ig = 1, g-1
     if(gas_cap(ig,z)*gas_drarr(z) &
          *help >= prt_taulump*gas_curvcent(z) &
          .and. g-ig <= gas_epslump) then
        glump=glump+1
        glumps(glump)=ig
     else
        glumps(gunlump)=ig
        gunlump=gunlump-1
     endif
  enddo
  do ig = g, gas_ng
     if(gas_cap(ig,z)*gas_drarr(z) &
          *help >= prt_taulump*gas_curvcent(z) &
          .and. ig-g <= gas_epslump) then
        glump=glump+1
        glumps(glump)=ig
     else
        glumps(gunlump)=ig
        gunlump=gunlump-1
     endif
  enddo
  if(glump==0) then
     glump=1
     glumps(1)=g
!
     forall(ig=2:g) glumps(ig)=ig-1
     forall(ig=g+1:gas_ng) glumps(ig)=ig
!
  endif
!
!-- lumping
  do ig = 1, glump
     iig = glumps(ig)
     specig = gas_siggrey(z)*gas_emitprob(iig,z)/&
          gas_cap(iig,z)
     speclump = speclump+specig
  enddo
!
  if(speclump>0d0) then
     do ig = 1, glump
        iig = glumps(ig)
        specig = gas_siggrey(z)*gas_emitprob(iig,z)/&
             gas_cap(iig,z)
!-- emission lump
        emitlump = emitlump+gas_emitprob(iig,z)
!-- Planck x-section lump
        caplump = caplump+specig*gas_cap(iig,z)/speclump
!-- inward leakage lump
        opacleakllump = opacleakllump+specig*&
             gas_opacleakl(iig,z)/speclump
!-- outward leakage lump
        opacleakrlump = opacleakrlump+specig*&
             gas_opacleakr(iig,z)/speclump
     enddo
  else
     emitlump = gas_emitprob(g,z)
     opacleakllump = gas_opacleakl(g,z)
     opacleakrlump = gas_opacleakr(g,z)
     caplump = gas_cap(g,z)
  endif
!
!-------------------------------------------------------------
!
  !
  denom = opacleakllump+opacleakrlump
  denom = denom+(1d0-alpeff)*(1d0-emitlump)*&
       (1d0-gas_fcoef(z))*caplump
!--add analog term
  if(prt_isddmcanlog) then
     denom = denom+gas_fcoef(z)*caplump
  endif
  !write(*,*) gas_emitprob(g,z),g
  r1 = rand()
  prt_tlyrand = prt_tlyrand+1
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-t
  ddmct = min(tau,tcensus)
!
!-- calculating energy depostion and density
  !
  if(.not.prt_isddmcanlog) then
     gas_edep(z)=gas_edep(z)+E*(1d0-exp(-gas_fcoef(z) &!{{{
          *caplump*pc_c*ddmct))
!--
!-- must use speclump...
     if(gas_fcoef(z)*caplump*gas_drarr(z)*help>1d-6) then
        do ig = 1, glump
           iig=glumps(ig)
           gas_eraddens(iig,z)= &
                gas_eraddens(iig,z)+E* &
                (1d0-exp(-gas_fcoef(z)*caplump*pc_c*ddmct))/ &
                (gas_fcoef(z)*caplump*pc_c*tsp_dt) &
                /real(glump)
        enddo
     else
        do ig = 1, glump
           iig=glumps(ig)
           gas_eraddens(iig,z)= &
                gas_eraddens(iig,z)+E*ddmct/tsp_dt &
                /real(glump)
        enddo
     endif
     E=E*exp(-gas_fcoef(z)*caplump*pc_c*ddmct)

!!}}}
  else
     !
     do ig = 1, glump
        iig=glumps(ig)
        gas_eraddens(iig,z)= &
             gas_eraddens(iig,z)+E*ddmct/tsp_dt &
             /real(glump)
     enddo
!
  endif
  !
  t = t+ddmct
!
!
!-- Recalculating histogram sum (rev. 120)
  denom = opacleakllump+opacleakrlump !+gas_fcoef(z)*gas_cap(g,z)
  denom = denom+(1d0-alpeff)*(1d0-emitlump)*&
       (1d0-gas_fcoef(z))*caplump
!--add analog term
  if(prt_isddmcanlog) then
     denom=denom+gas_fcoef(z)*caplump
  endif

  if (ddmct == tau) then
     r1 = rand()!{{{
     prt_tlyrand = prt_tlyrand+1
!-- right leak probability
     PR = opacleakrlump/denom
!-- left leak probability
     PL = opacleakllump/denom
!-- absorption probability
     if(prt_isddmcanlog) then
        PA = gas_fcoef(z)*caplump/denom
     else
        PA = 0d0
     endif
!
!-- left leakage sample
     if (0.0d0<=r1 .and. r1<PL) then
!{{{
!--
        if (z == 1) then
           if(gas_isshell) then
              vacnt = .true.
              prt_done = .true.
              gas_eleft = gas_eleft+E
           else
              write(6,*) 'Non-physical left leakage', g, gas_wl(g+1), wl
           endif
!
!
        elseif(speclump>0d0) then
!
!-- sample adjacent group (assumes aligned g bounds)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           denom2 = 0d0
           do ig= 1, glump
              iig = glumps(ig)
!              iig = ig
              specig = gas_siggrey(z)*gas_emitprob(iig,z)/&
                   gas_cap(iig,z)
              if((r1>=denom2).and. &
                   (r1<denom2+specig* &
                   gas_opacleakl(iig,z) &
                   /(speclump*opacleakllump))) exit

              denom2 = denom2+specig* &
                   gas_opacleakl(iig,z) &
                   /(speclump*opacleakllump)
!
           enddo
           if((gas_sig(z-1)+gas_cap(iig,z-1))*gas_drarr(z-1) &
                *help >= prt_tauddmc*gas_curvcent(z-1)) then
              z = z-1
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
              g = iig
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
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
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = -max(r1,r2)
!
!-- doppler and aberration corrections
              if(gas_isvelocity) then
                 mu = (mu+r/pc_c)/(1d0+r*mu/pc_c)
!-- velocity effects accounting
                 gas_evelo=gas_evelo+E*(1d0-1d0/(1.0-r*mu/pc_c))
!
                 E = E/(1.0-r*mu/pc_c)
                 E0 = E0/(1.0-r*mu/pc_c)
                 wl = wl*(1.0-r*mu/pc_c)
              endif
!
!-- group reset
              g = iig
!
           endif
!
!
        elseif ((gas_sig(z-1)+gas_cap(g,z-1))*gas_drarr(z-1) &
             *help >= prt_tauddmc*gas_curvcent(z-1)) then
           z = z-1
!--
!
        else

           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
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
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = -max(r1,r2)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1d0+r*mu/pc_c)
!-- velocity effects accounting
              gas_evelo=gas_evelo+E*(1d0-1d0/(1.0-r*mu/pc_c))
!
              E = E/(1.0-r*mu/pc_c)
              E0 = E0/(1.0-r*mu/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif!}}}
!-- right leakage sample
     elseif (PL<=r1 .and. r1<PL+PR) then
!!{{{
        if (z == gas_nr) then
           vacnt = .true.
           prt_done = .true.
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = max(r1,r2)
!-- outbound luminosity tally
           if(gas_isvelocity) then
              gas_eright = gas_eright+E
!               gas_luminos(gminlump:gmaxlump)=gas_luminos(gminlump:gmaxlump)+(E/tsp_dt)* &
!                    (mu+gas_rarr(gas_nr+1)/pc_c)/ &
!                    (1.0+gas_rarr(gas_nr+1)*mu/pc_c)/real(gmaxlump-gminlump+1)
              do ig = 1, glump
                 iig=glumps(ig)
                 gas_luminos(iig)= &
                      gas_luminos(iig)+&
                      (E/tsp_dt)* &
                      (1.0+gas_rarr(gas_nr+1)*mu/pc_c) &
                      /real(glump)
              enddo
           else
              gas_eright = gas_eright+E
              do ig = 1, glump
                 iig=glumps(ig)
                 gas_luminos(iig)=gas_luminos(iig)+&
                      (E/tsp_dt)*mu/real(glump)
              enddo
           endif
!
!
        elseif(speclump>0d0) then
!
!-- sample adjacent group (assumes aligned g bounds)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           denom2 = 0d0
           do ig= 1, glump
              iig = glumps(ig)
              specig = gas_siggrey(z)*gas_emitprob(iig,z)/&
                   gas_cap(iig,z)
              if((r1>=denom2).and. &
                   (r1<denom2+specig* &
                   gas_opacleakr(iig,z) &
                   /(speclump*opacleakrlump))) exit

              denom2 = denom2+specig* &
                   gas_opacleakr(iig,z) &
                   /(speclump*opacleakrlump)
           enddo
           if((gas_sig(z+1)+gas_cap(iig,z+1))*gas_drarr(z+1) &
                *help >= prt_tauddmc*gas_curvcent(z+1)) then
!
              z = z+1
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
!--
!
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
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
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = max(r1,r2)
!
!-- doppler and aberration corrections
              if(gas_isvelocity) then
                 mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
!-- velocity effects accounting
                 gas_evelo=gas_evelo+E*(1d0-1d0/(1.0-r*mu/pc_c))
!
                 E = E/(1.0-r*mu/pc_c)
                 E0 = E0/(1.0-r*mu/pc_c)
                 wl = wl*(1.0-r*mu/pc_c)
              endif
!
!-- group reset
              g = iig
!
           endif
!
!
        elseif ((gas_sig(z+1)+gas_cap(g,z+1))*gas_drarr(z+1) &
             *help >= prt_tauddmc*gas_curvcent(z+1)) then
!
           z = z+1
!--
!
        else
            r1 = rand()
            prt_tlyrand = prt_tlyrand+1
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
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
!-- velocity effects accounting
              gas_evelo=gas_evelo+E*(1d0-1d0/(1.0-r*mu/pc_c))
!
              E = E/(1.0-r*mu/pc_c)
              E0 = E0/(1.0-r*mu/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif!}}}
!-- absorption sample
     elseif (PL+PR<=r1 .and. r1<PL+PR+PA) then
        vacnt = .true.!{{{
        prt_done = .true.
        gas_edep(z) = gas_edep(z)+E
!
!!}}}
!
!-- fictitious scattering sample (physical scattering currently elastic)
!
     else
        denom2 = 0d0!{{{
!-- this may need to be changed for g95 compiler
!         do ig = 1, gminlump-1
! !           if(ig.ne.g) then
!            denom2 = denom2+gas_emitprob(ig,z)
! !           endif
!         enddo
! !
!         do ig = gmaxlump+1, gas_ng
!            denom2 = denom2+gas_emitprob(ig,z)
!         enddo
        do ig = gas_ng,glump+1,-1
           iig=glumps(ig)
           denom2 = denom2+gas_emitprob(iig,z)
        enddo
!
        denom3 = 0d0
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
!         do ig = 1, gas_ng
!            if(ig<gminlump.or.ig>gmaxlump) then
!               iig = ig
!               if((r1>=denom3).and.(r1<denom3+gas_emitprob(ig,z)/denom2)) exit
!               denom3 = denom3+gas_emitprob(ig,z)/denom2
!            endif
!         enddo
        do ig = gas_ng,glump+1,-1
           iig=glumps(ig)
           if((r1>=denom3).and.(r1<denom3+gas_emitprob(iig,z)/denom2)) exit
           denom3 = denom3+gas_emitprob(iig,z)/denom2
        enddo
!
        !write(*,*) 'Scatter: ',g,'to ',iig
        g = iig
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        wl = 1d0/((1d0-r1)/gas_wl(g)+r1/gas_wl(g+1))
        ! during DDMC phase, wavelength is only a placeholder for group 
        !wl = 0.5d0*(gas_wl(g)+gas_wl(g+1))
        !
        if ((gas_sig(z)+gas_cap(g,z))*gas_drarr(z) &
             *help >= prt_tauddmc*gas_curvcent(z)) then
           hyparam = 2
!
!--
!
        else
           hyparam = 1
!
!-- direction sampled isotropically           
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = 1.0-2.0*r1
!-- position sampled uniformly
            r1 = rand()
            prt_tlyrand = prt_tlyrand+1
            r = (r1*gas_rarr(z+1)**3+(1.0-r1)*gas_rarr(z)**3)**(1.0/3.0)
!
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r/pc_c)/(1.0+r*mu/pc_c)
!-- velocity effects accounting
              gas_evelo=gas_evelo+E*(1d0-1d0/(1.0-r*mu/pc_c))
!
              E = E/(1.0-mu*r/pc_c)
              E0 = E0/(1.0-mu*r/pc_c)
              wl = wl*(1.0-r*mu/pc_c)
           endif
        endif
     !else
     !   stop 'diffusion1: invalid histogram sample'!}}}
     endif!}}}
  else
     !gas_eraddens(g,z)=gas_eraddens(g,z)+E*ddmct/tsp_dt
     prt_done = .true.
     gas_numcensus(z)=gas_numcensus(z)+1
!     gas_erad = gas_erad+E
  endif

!--
  if(E/E0<1d-6.and..not.vacnt) then
     r1 = rand()
     prt_tlyrand=prt_tlyrand+1
     if(r1<0.5d0) then
        vacnt=.true.
        prt_done=.true.
        gas_edep(z)=gas_edep(z)+E
     else
        E=2d0*E
        E0=2d0*E0
     endif
  endif


!-- temporarily deprecated ---------------------------------------------
!
!-- wl thermal resample example
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
!                 prt_tlyrand = prt_tlyrand+1
!            r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            xx0 = (1d0-r1)*x1+r1*x2
!            do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!               r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!               r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!               xx0 = (1d0-r1)*x1+r1*x2
!            enddo
!            wl = pc_h*pc_c/(xx0*pc_kb*gas_temp(z))
!
!
!-- position sampled from source tilt
!            r1 = 0d0
!            r2 = 1d0
!            uul = gas_tempb(z)**4
!            uur = gas_tempb(z+1)**4
!            uumax = max(uul,uur)
!            do while (r2 > r1)
!               r3 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!               r0 = (r3*gas_rarr(z+1)**3+(1.0-r3)*gas_rarr(z)**3)**(1.0/3.0)
!               r3 = (r0-gas_rarr(z))/gas_drarr(z)
!               r1 = (r3*uur+(1d0-r3)*uul)/uumax
!               r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            enddo
!            r = r0

end subroutine diffusion1
