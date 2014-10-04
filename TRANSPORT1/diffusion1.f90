subroutine diffusion1(z,wl,r,mu,t,E,E0,hyparam,vacnt,partnum)

  use gasgridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  implicit none
!
  integer, intent(in) :: partnum !ipart
  integer, intent(inout) :: z, hyparam !,g
  real*8, intent(inout) :: r, mu, t, E, E0, wl
  logical, intent(inout) :: vacnt
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: deleff=0.38d0
  real*8,parameter :: alpeff=0d0
!
  integer :: ig, iig, g
  logical :: lhelp
  integer,external :: binsrch
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus, PR, PL, PA
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleakllump, opacleakrlump, mfphelp, ppl, ppr
  real*8 :: resopacleakl, resopacleakr
  integer :: glump, gunlump
  integer :: glumps(gas_ng)
  real*8 :: glumpinv,dtinv,capinv(gas_ng)
  real*8 :: help

!--------------------------------------------------------------
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  capinv = 1d0/gas_cap(:,z,1,1)

!
!-- set expansion helper
  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif


!-- find group
  g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
  if(g>gas_ng.or.g<1) then
     !particle out of wlgrid bound
     if(g>gas_ng) then
        g=gas_ng
        wl=gas_wl(gas_ng+1)
     elseif(g<1) then
        g=1
        wl=gas_wl(1)
     else
        stop 'diffusion1: missing particle group'
     endif
  endif


!-- lump testing ---------------------------------------------
!
  glump = 0
  gunlump = gas_ng
  glumps = 0
!
!-- find lumpable groups
  if(gas_cap(g,z,1,1)*gas_dxarr(z)*thelp>=prt_taulump) then
     do ig = 1, g-1
        if(gas_cap(ig,z,1,1)*gas_dxarr(z) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=ig
        else
           glumps(gunlump)=ig
           gunlump=gunlump-1
        endif
     enddo
     do ig = g, gas_ng
        if(gas_cap(ig,z,1,1)*gas_dxarr(z) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=ig
        else
           glumps(gunlump)=ig
           gunlump=gunlump-1
        endif
     enddo
  endif
!
  if(glump==0) then
     glump=1
     glumps(1)=g
!
     forall(ig=2:g) glumps(ig)=ig-1
     forall(ig=g+1:gas_ng) glumps(ig)=ig
!
  endif
  glumpinv = 1d0/glump
!
!-- lumping
  speclump = 0d0
  do ig = 1, glump
     iig = glumps(ig)
     specig = gas_siggrey(z,1,1)*gas_emitprob(iig,z,1,1)*capinv(iig)
     speclump = speclump+specig
  enddo
  if(speclump>0d0.and.glump>1) then
     speclump = 1d0/speclump
  else
     speclump = 0d0
  endif
!
  emitlump = 0d0
  caplump = 0d0
  if(speclump>0d0) then
!
!-- calculating lumped values
     do ig = 1, glump
        iig = glumps(ig)
        specig = gas_siggrey(z,1,1)*gas_emitprob(iig,z,1,1)*capinv(iig)
!-- emission lump
        emitlump = emitlump+gas_emitprob(iig,z,1,1)
!-- Planck x-section lump
        caplump = caplump+specig*gas_cap(iig,z,1,1)*speclump
     enddo
!-- leakage opacities
     opacleakllump = gas_opacleak(1,z,1,1)
     opacleakrlump = gas_opacleak(2,z,1,1)
  else
!
!-- calculating unlumped values
     emitlump = gas_emitprob(g,z,1,1)
     caplump = gas_cap(g,z,1,1)
!-- inward
     if(z==1) then
        opacleakllump = 0d0
     elseif((gas_cap(g,z-1,1,1)+ &
          gas_sig(z-1,1,1))*gas_dxarr(z-1)*thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (gas_cap(g,z,1,1)+gas_sig(z,1,1))*gas_dxarr(z)*thelp
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakllump= 1.5d0*ppl*(thelp*gas_xarr(z))**2/ &
             (thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
     else
!-- DDMC interior
        mfphelp = ((gas_sig(z,1,1)+gas_cap(g,z,1,1))*gas_dxarr(z)+&
             (gas_sig(z-1,1,1)+gas_cap(g,z-1,1,1))*gas_dxarr(z-1))*thelp
        opacleakllump=2.0d0*(thelp*gas_xarr(z))**2/ &
             (mfphelp*thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
     endif
!
!-- outward
     if(z==gas_nx) then
        lhelp = .true.
     else
        lhelp = (gas_cap(g,z+1,1,1)+ &
           gas_sig(z+1,1,1))*gas_dxarr(z+1)*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (gas_cap(g,z,1,1)+gas_sig(z,1,1))*gas_dxarr(z)*thelp
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakrlump=1.5d0*ppr*(thelp*gas_xarr(z+1))**2/ &
             (thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
     else
!-- DDMC interior
        mfphelp = ((gas_sig(z,1,1)+gas_cap(g,z,1,1))*gas_dxarr(z)+&
             (gas_sig(z+1,1,1)+gas_cap(g,z+1,1,1))*gas_dxarr(z+1))*thelp
        opacleakrlump=2.0d0*(thelp*gas_xarr(z+1))**2/ &
             (mfphelp*thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
     endif
  endif
!
!-------------------------------------------------------------
!

!-- calculate time to census or event
  denom = opacleakllump+opacleakrlump+&
       (1d0-alpeff)*(1d0-emitlump)*(1d0-gas_fcoef(z,1,1))*caplump
  if(prt_isddmcanlog) then
     denom = denom+gas_fcoef(z,1,1)*caplump
  endif

  r1 = rand()
  prt_tlyrand = prt_tlyrand+1
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-t
  ddmct = min(tau,tcensus)
!
!-- calculating energy depostion and density
  !
  if(prt_isddmcanlog) then
     gas_eraddens(z,1,1)= gas_eraddens(z,1,1)+E*ddmct*dtinv
  else
     gas_edep(z,1,1) = gas_edep(z,1,1)+E*(1d0-exp(-gas_fcoef(z,1,1) &!{{{
          *caplump*pc_c*ddmct))
!--
     if(gas_fcoef(z,1,1)*caplump*gas_dxarr(z)*thelp>1d-6) then
        help = 1d0/(gas_fcoef(z,1,1)*caplump)
        gas_eraddens(z,1,1)= &
             gas_eraddens(z,1,1)+E* &
             (1d0-exp(-gas_fcoef(z,1,1)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        gas_eraddens(z,1,1) = gas_eraddens(z,1,1)+E*ddmct*dtinv
     endif
     E=E*exp(-gas_fcoef(z,1,1)*caplump*pc_c*ddmct)
!!}}}
  endif


!-- updating particle time
  t = t+ddmct


!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     prt_done = .true.
     gas_numcensus(z,1,1)=gas_numcensus(z,1,1)+1


!-- otherwise, perform event
  else
     r1 = rand()!{{{
     prt_tlyrand = prt_tlyrand+1
     help = 1d0/denom
!-- right leak probability
     PR = opacleakrlump*help

!-- left leak probability
     PL = opacleakllump*help

!-- absorption probability
     if(prt_isddmcanlog) then
        PA = gas_fcoef(z,1,1)*caplump*help
     else
        PA = 0d0
     endif

!-- absorption sample
     if(r1>=0d0 .and. r1<PA) then
        vacnt = .true.
        prt_done = .true.
        gas_edep(z,1,1) = gas_edep(z,1,1)+E


!-- left leakage sample
     elseif (r1>=PA .and. r1<PA+PL) then
!{{{
!-- checking if at inner bound
        if (z == 1) then
           stop 'diffusion1: non-physical inward leakage'

!-- sample adjacent group (assumes aligned g bounds)
        else

           if(speclump>0d0) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              denom2 = 0d0
              help = 1d0/opacleakllump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(z,1,1)*gas_emitprob(iig,z,1,1)*capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,z-1,1,1)+ &
                      gas_sig(z-1,1,1))*gas_dxarr(z-1)*thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,z,1,1)+gas_sig(z,1,1))*gas_dxarr(z)*thelp
                    ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakl = 1.5d0*ppl*(thelp*gas_xarr(z))**2/ &
                         (thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(z,1,1)+gas_cap(iig,z,1,1))*gas_dxarr(z)+&
                         (gas_sig(z-1,1,1)+gas_cap(iig,z-1,1,1))*gas_dxarr(z-1))*thelp
                    resopacleakl = 2.0d0*(thelp*gas_xarr(z))**2/ &
                         (mfphelp*thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakl*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakl*speclump*help
              enddo
           else
              iig = g
           endif


           if((gas_sig(z-1,1,1)+gas_cap(iig,z-1,1,1))*gas_dxarr(z-1) &
                *thelp >= prt_tauddmc) then
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
              gas_methodswap(z,1,1)=gas_methodswap(z,1,1)+1
!
!-- location set right bound of left cell
              r = gas_xarr(z)
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
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
                 help = 1d0/(1.0-r*mu*cinv)
                 gas_evelo=gas_evelo+E*(1d0 - help)
!
                 E = E*help
                 E0 = E0*help
                 wl = wl*(1.0-r*mu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif

        endif!}}}


!-- right leakage sample
     elseif (r1>=PA+PL .and. r1<PA+PL+PR) then
!!{{{
!-- checking if at outer bound
        if (z == gas_nx) then
           vacnt = .true.
           prt_done = .true.
           gas_eright = gas_eright+E
!-- outbound luminosity tally
           if(gas_isvelocity) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = max(r1,r2)
              if(speclump>0d0) then
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 denom2 = 0d0
                 help = 1d0/opacleakrlump
                 do ig = 1, glump
                    iig=glumps(ig)
                    specig = gas_siggrey(z,1,1)*gas_emitprob(iig,z,1,1)*capinv(iig)
!-- calculating resolved leakage opacities
                    mfphelp = (gas_cap(iig,z,1,1)+gas_sig(z,1,1))*gas_dxarr(z)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakr = 1.5d0*ppr*(thelp*gas_xarr(z+1))**2/ &
                         (thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
                    if((r1>=denom2).and. &
                         (r1<denom2+specig*resopacleakr*speclump*help)) exit
                    denom2 = denom2+specig*resopacleakr*speclump*help
!                    gas_luminos(iig) = gas_luminos(iig)+&
!                         E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv) * glumpinv
                 enddo
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 wl=1d0/(r1/gas_wl(iig+1) + (1d0-r1)/gas_wl(iig))
!-- changing from comoving frame to observer frame
                 wl = wl/(1d0+mu*gas_xarr(gas_nx+1)*cinv)
!-- obtaining spectrum (lab) group
                 iig = binsrch(wl,gas_wl,gas_ng+1,in_ng)
                 if(iig>gas_ng.or.iig<1) then
                    if(iig>gas_ng) then
                       iig=gas_ng
                       wl=gas_wl(gas_ng+1)
                    else
                       iig=1
                       wl=gas_wl(1)
                    endif
                 endif
                 gas_luminos(iig) = gas_luminos(iig)+&
                      E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv)
                 gas_lumdev(iig) = gas_lumdev(iig)+&
                      (E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv))**2
                 gas_lumnum(iig) = gas_lumnum(iig)+1
              else
                 iig = g
                 gas_luminos(iig) = gas_luminos(iig)+&
                      E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv)
                 gas_lumdev(iig) = gas_lumdev(iig)+&
                      (E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv))**2
                 gas_lumnum(iig) = gas_lumnum(iig)+1
              endif
!               gas_luminos(iig) = gas_luminos(iig)+&
!                    E*dtinv*(1.0+gas_xarr(gas_nx+1)*mu*cinv)
           else
              gas_eright = gas_eright+E
              do ig = 1, glump
                 iig=glumps(ig)
                 gas_luminos(iig)=gas_luminos(iig)+&
                      (E*dtinv) * glumpinv
                 gas_lumdev(iig) = gas_lumdev(iig)+&
                      (E*dtinv*glumpinv)**2
                 gas_lumnum(iig) = gas_lumnum(iig)+1
              enddo
           endif
!
!
        else
!
!-- sample adjacent group (assumes aligned g bounds)
           if(speclump>0d0) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              denom2 = 0d0
              help = 1d0/opacleakrlump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(z,1,1)*gas_emitprob(iig,z,1,1)*capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,z+1,1,1)+ &
                      gas_sig(z+1,1,1))*gas_dxarr(z+1)*thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,z,1,1)+gas_sig(z,1,1))*gas_dxarr(z)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakr = 1.5d0*ppr*(thelp*gas_xarr(z+1))**2/ &
                         (thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(z,1,1)+gas_cap(iig,z,1,1))*gas_dxarr(z)+&
                         (gas_sig(z+1,1,1)+gas_cap(iig,z+1,1,1))*gas_dxarr(z+1))*thelp
                    resopacleakr = 2.0d0*(thelp*gas_xarr(z+1))**2/ &
                         (mfphelp*thelp**3*(gas_xarr(z+1)**3-gas_xarr(z)**3))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakr*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakr*speclump*help
              enddo
           else
              iig = g
           endif


           if((gas_sig(z+1,1,1)+gas_cap(iig,z+1,1,1))*gas_dxarr(z+1) &
                *thelp >= prt_tauddmc) then
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
              gas_methodswap(z,1,1)=gas_methodswap(z,1,1)+1
!
!-- location set left bound of right cell
              r = gas_xarr(z+1)
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
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
                 help = 1d0/(1.0-r*mu*cinv)
                 gas_evelo=gas_evelo+E*(1d0 - help)
!
                 E = E*help
                 E0 = E0*help
                 wl = wl*(1.0-r*mu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif
        endif!!}}}


!-- effective scattering sample
     else
!         denom2 = 0d0!{{{
!         do ig = gas_ng,glump+1,-1
!            iig=glumps(ig)
!            denom2 = denom2+gas_emitprob(iig,z,1,1)
!         enddo

        denom2 = 1d0-emitlump
        help = 1d0/denom2
!
        denom3 = 0d0
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1

        do ig = gas_ng,glump+1,-1
           iig=glumps(ig)
           if((r1>=denom3).and.(r1<denom3+gas_emitprob(iig,z,1,1)*help)) exit
           denom3 = denom3+gas_emitprob(iig,z,1,1)*help
        enddo
!
        g = iig
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        wl = 1d0/((1d0-r1)/gas_wl(g) + r1/gas_wl(g+1))

        if ((gas_sig(z,1,1)+gas_cap(g,z,1,1))*gas_dxarr(z) &
             *thelp >= prt_tauddmc) then
           hyparam = 2
        else
           hyparam = 1
           gas_methodswap(z,1,1)=gas_methodswap(z,1,1)+1
!-- direction sampled isotropically           
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = 1.0-2.0*r1
!-- position sampled uniformly
            r1 = rand()
            prt_tlyrand = prt_tlyrand+1
            r = (r1*gas_xarr(z+1)**3+(1.0-r1)*gas_xarr(z)**3)**(1.0/3.0)
!-- doppler and aberration corrections
           if(gas_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
              help = 1d0/(1.0-r*mu*cinv)
              gas_evelo=gas_evelo+E*(1d0 - help)
!
              E = E*help
              E0 = E0*help
              wl = wl*(1.0-r*mu*cinv)
           endif
        endif

     endif!}}}
  endif


!-- check for particle termination
  if(E<1d-6*E0.and..not.vacnt) then
     r1 = rand()
     prt_tlyrand=prt_tlyrand+1
     if(r1<0.5d0) then
        vacnt=.true.
        prt_done=.true.
        gas_edep(z,1,1)=gas_edep(z,1,1)+E
     else
!-- weight addition accounted for in external source
        gas_eext=gas_eext+E
!
        E=2d0*E
        E0=2d0*E0
     endif
  endif


end subroutine diffusion1


