subroutine diffusion2(ptcl,isvacant)

  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  use totalsmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch
!
  integer :: ig, iig, iiig, imu
  logical :: lhelp
  real*8 :: r1, r2, thelp, mu0
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: elabfact, dirdotu
  real*8 :: PU, PD, PR, PL, PA
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleakllump, opacleakrlump
  real*8 :: opacleakdlump, opacleakulump
  real*8 :: mfphelp, ppl, ppr
  real*8 :: resopacleakl, resopacleakr
  real*8 :: resopacleakd, resopacleaku
  integer :: glump, gunlump
  integer :: glumps(grd_ng)
  real*8 :: dtinv,capinv(grd_ng)
  real*8 :: help

  integer,pointer :: ix,iy
  real*8,pointer :: x,y,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2-grd_xarr(l)**2
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

  ix => ptcl%ix
  iy => ptcl%iy
  x => ptcl%x
  y => ptcl%y
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!--------------------------------------------------------------
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  capinv = 1d0/grd_cap(:,ix,iy,1)

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- looking up initial group
  ig = binsrch(wl,grd_wl,grd_ng+1,in_ng)
!-- checking group bounds
  if(ig>grd_ng.or.ig<1) then
     if(ig>grd_ng) then
        ig = grd_ng
     elseif(ig<1) then
        ig = 1
     else
        stop 'diffusion2 (1): particle group invalid'
     endif
  endif

!-- lump testing ---------------------------------------------
!
  glump = 0
  gunlump = grd_ng
  glumps = 0
!
!-- find lumpable groups
  if(grd_cap(ig,ix,iy,1)*min(dx(ix),dy(iy)) * &
       thelp>=prt_taulump) then
     do iig = 1, ig-1
        if(grd_cap(iig,ix,iy,1)*min(dx(ix),dy(iy)) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
     do iig = ig, grd_ng
        if(grd_cap(iig,ix,iy,1)*min(dx(ix),dy(iy)) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
  endif
!
  if(glump==0) then
     glump=1
     glumps(1)=ig
!
     forall(iig=2:ig) glumps(iig)=iig-1
     forall(iig=ig+1:grd_ng) glumps(iig)=iig
!
  endif

!
!-- lumping
  speclump = 0d0
  do iig = 1, glump
     iiig = glumps(iig)
     specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1)*capinv(iiig)
     speclump = speclump+specig
  enddo
  if(speclump>0d0.and.glump>1) then
     speclump = 1d0/speclump
  else
     speclump = 0d0
  endif

  emitlump = 0d0
  caplump = 0d0
  if(speclump>0d0) then
!
!-- calculating lumped values
     do iig = 1, glump
        iiig = glumps(iig)
        specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1)*capinv(iiig)
!-- emission lump
        emitlump = emitlump+grd_emitprob(iiig,ix,iy,1)
!-- Planck x-section lump
        caplump = caplump+specig*grd_cap(iiig,ix,iy,1)*speclump
     enddo
!-- leakage opacities
     opacleakllump = grd_opacleak(1,ix,iy,1)
     opacleakrlump = grd_opacleak(2,ix,iy,1)
     opacleakdlump = grd_opacleak(3,ix,iy,1)
     opacleakulump = grd_opacleak(4,ix,iy,1)
  else
!
!-- calculating unlumped values
     emitlump = grd_emitprob(ig,ix,iy,1)
     caplump = grd_cap(ig,ix,iy,1)

!-- inward
     if(ix==1) then
        opacleakllump = 0d0
     elseif((grd_cap(ig,ix-1,iy,1)+ &
          grd_sig(ix-1,iy,1))*min(dx(ix-1),dy(iy))* &
          thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1))*dx(ix)*thelp
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakllump= ppl*(thelp*grd_xarr(ix))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dx(ix)+&
             (grd_sig(ix-1,iy,1)+grd_cap(ig,ix-1,iy,1))*dx(ix-1))*thelp
        opacleakllump=(4d0/3d0)*(thelp*grd_xarr(ix))/ &
             (mfphelp*thelp**2*dx2(ix))
     endif

!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix+1,iy,1)+ &
           grd_sig(ix+1,iy,1))*min(dx(ix+1),dy(iy))* &
           thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1))*dx(ix)*thelp
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakrlump= ppr*(thelp*grd_xarr(ix+1))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dx(ix)+&
             (grd_sig(ix+1,iy,1)+grd_cap(ig,ix+1,iy,1))*dx(ix+1))*thelp
        opacleakrlump=(4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
             (mfphelp*thelp**2*dx2(ix))
     endif

!-- downward
     if(iy==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix,iy-1,1)+ &
             grd_sig(ix,iy-1,1))*min(dx(ix),dy(iy-1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1))*dy(iy)*thelp
        ppl = 4d0/(3d0*help+6d0*pc_dext)
        opacleakdlump=0.5d0*ppl/(thelp*dy(iy))
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dy(iy)+&
             (grd_sig(ix,iy-1,1)+grd_cap(ig,ix,iy-1,1))*dy(iy-1))*thelp
        opacleakdlump=(2d0/3d0)/(help*dy(iy)*thelp)
     endif

!-- upward
     if(iy==grd_ny) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix,iy+1,1)+ &
             grd_sig(ix,iy+1,1))*min(dx(ix),dy(iy+1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1))*dy(iy)*thelp
        ppr = 4d0/(3d0*help+6d0*pc_dext)
        opacleakulump=0.5d0*ppr/(thelp*dy(iy))
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dy(iy)+&
             (grd_sig(ix,iy+1,1)+grd_cap(ig,ix,iy+1,1))*dy(iy+1))*thelp
        opacleakulump=(2d0/3d0)/(help*dy(iy)*thelp)
     endif
  endif
!
!-------------------------------------------------------------
!

!-- calculating time to census or event
  denom = opacleakllump+opacleakrlump + &
       opacleakdlump+opacleakulump + &
       (1d0-emitlump)*(1d0-grd_fcoef(ix,iy,1))*caplump
  if(prt_isddmcanlog) then
     denom = denom+grd_fcoef(ix,iy,1)*caplump
  endif

  r1 = rand()
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  !
  if(prt_isddmcanlog) then
     grd_eraddens(ix,iy,1)= grd_eraddens(ix,iy,1)+e*ddmct*dtinv
  else
     grd_edep(ix,iy,1) = grd_edep(ix,iy,1)+e*(1d0-exp(-grd_fcoef(ix,iy,1) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ix,iy,1)*caplump*min(dx(ix),dy(iy)) * &
          thelp>1d-6) then
        help = 1d0/(grd_fcoef(ix,iy,1)*caplump)
        grd_eraddens(ix,iy,1)= &
             grd_eraddens(ix,iy,1)+e* &
             (1d0-exp(-grd_fcoef(ix,iy,1)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        grd_eraddens(ix,iy,1) = grd_eraddens(ix,iy,1)+e*ddmct*dtinv
     endif
     e=e*exp(-grd_fcoef(ix,iy,1)*caplump*pc_c*ddmct)
!!}}}
  endif

!-- updating particle time
  ptcl%t = ptcl%t+ddmct

!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     prt_done = .true.
     grd_numcensus(ix,iy,1)=grd_numcensus(ix,iy,1)+1

!-- otherwise, perform event
  else
     r1 = rand()!{{{
     help = 1d0/denom

!-- outward leak probability
     PR = opacleakrlump*help
!-- inward leak probability
     PL = opacleakllump*help
!-- upward leak probability
     PU = opacleakulump*help
!-- downward leak probability
     PD = opacleakdlump*help
!-- absorption probability
     if(prt_isddmcanlog) then
        PA = grd_fcoef(ix,iy,1)*caplump*help
     else
        PA = 0d0
     endif

!-- absorption sample
     if(r1>=0d0 .and. r1<PA) then
        isvacant = .true.
        prt_done = .true.
        grd_edep(ix,iy,1) = grd_edep(ix,iy,1)+e

!-- inward leakage sample
     elseif (r1>=PA .and. r1<PA+PL) then
!{{{
!-- checking if at inner bound
        if (ix == 1) then
           stop 'diffusion1: non-physical inward leakage'

!-- sample adjacent group (assumes aligned ig bounds)
        else
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakllump
              do iig= 1, glump
                 iiig = glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 if((grd_cap(iiig,ix-1,iy,1) + &
                      grd_sig(ix-1,iy,1))*min(dx(ix-1),dy(iy)) * &
                      thelp<prt_tauddmc) then
!-- IMC interface
                    mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                         dx(ix)*thelp
                    ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakl = ppl*(thelp*grd_xarr(ix))/ &
                         (thelp**2*dx2(ix))
                 else
!-- DDMC interface
                    mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                         dx(ix)+(grd_sig(ix-1,iy,1) + &
                         grd_cap(iiig,ix-1,iy,1))*dx(ix-1))*thelp
                    resopacleakl = (4d0/3d0)*(thelp*grd_xarr(ix))/ &
                         (mfphelp*thelp**2*dx2(ix))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakl*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakl*speclump*help
              enddo
           else
              iiig = ig
           endif

!-- updating properties
           if((grd_sig(ix-1,iy,1)+grd_cap(iiig,ix-1,iy,1)) * &
                min(dx(ix-1),dy(iy))*thelp>=prt_tauddmc) then
              ix = ix-1
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
              ig = iiig
           else
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
!
!-- method changed to IMC
              ptcl%itype = 1
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
              x = grd_xarr(ix)

!-- current particle cell set to 1 left
              ix = ix-1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu0 = -max(r1,r2)
              r1 = rand()
              mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
              if(om<0d0) om = om+pc_pi2
!-- doppler and aberration corrections
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 e = e*(1d0+dirdotu*cinv)
                 e0 = e0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              ig = iiig
!
           endif

        endif!}}}

!
!-- outward leakage sample
     elseif (r1>=PA+PL .and. r1<PA+PL+PR) then
!!{{{
!-- checking if at outer bound
        if (ix == grd_nx) then
           isvacant = .true.
           prt_done = .true.
           tot_eright = tot_eright+e
!-- outbound luminosity tally
!-- sampling x, y
           r1 = rand()
           y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
           x = grd_xarr(grd_nx+1)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           mu0 = max(r1,r2)
           r1 = rand()
           mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
           if(om<0d0) om = om+pc_pi2
           if(grd_isvelocity) then
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              elabfact = 1d0+dirdotu*cinv
           else
              elabfact = 1d0
           endif
!-- determining outbound group
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakrlump
              do iig = 1, glump
                 iiig=glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                      dx(ix)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleakr = ppr*(thelp*grd_xarr(ix+1))/ &
                      (thelp**2*dx2(ix))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakr*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakr*speclump*help
              enddo
           else
              iiig = ig
           endif
!-- sampling outbound wavlength
           r1 = rand()
           wl=1d0/(r1/grd_wl(iiig+1) + (1d0-r1)/grd_wl(iiig))
!-- changing from comoving frame to observer frame
           if(grd_isvelocity) then
              mu = (mu+y*cinv)/elabfact
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
              wl = wl/elabfact
           endif
!-- obtaining spectrum (lab) group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           iiig = binsrch(wl,flx_wl,flx_ng+1,0)
           if(iiig>flx_ng.or.iiig<1) then
              if(iiig>flx_ng) then
                 iiig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 iiig=1
                 wl=flx_wl(1)
              endif
           endif
           flx_luminos(iiig,imu,1)=flx_luminos(iiig,imu,1)+&
                e*dtinv*elabfact
           flx_lumdev(iiig,imu,1)=flx_lumdev(iiig,imu,1)+&
                (e*dtinv*elabfact)**2
           flx_lumnum(iiig,imu,1)=flx_lumnum(iiig,imu,1)+1

        else
!-- outward leaking in domain
!
!-- sample adjacent group (assumes aligned ig bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakrlump
              do iig= 1, glump
                 iiig = glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 if((grd_cap(iiig,ix+1,iy,1)+grd_sig(ix+1,iy,1)) * &
                      min(dx(ix+1),dy(iy))*thelp<prt_tauddmc) then
!-- IMC interface
                    mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                         dx(ix)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakr = ppr*(thelp*grd_xarr(ix+1))/ &
                         (thelp**2*dx2(ix))
                 else
!-- DDMC interface
                    mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                         dx(ix)+&
                         (grd_sig(ix+1,iy,1)+grd_cap(iiig,ix+1,iy,1)) * &
                         dx(ix+1))*thelp
                    resopacleakr = (4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
                         (mfphelp*thelp**2*dx2(ix))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakr*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakr*speclump*help
              enddo
           else
              iiig = ig
           endif

!-- updating particle properties
           if((grd_sig(ix+1,iy,1)+grd_cap(iiig,ix+1,iy,1)) * &
                min(dx(ix+1),dy(iy)) &
                *thelp >= prt_tauddmc) then
              ix = ix+1
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
              ig = iiig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
!
!-- method changed to IMC
              ptcl%itype = 1
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
              x = grd_xarr(ix+1)

!-- current particle cell set to 1 left
              ix = ix+1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu0 = max(r1,r2)
              r1 = rand()
              mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
              if(om<0d0) om = om+pc_pi2
!-- doppler and aberration corrections
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 e = e*(1d0+dirdotu*cinv)
                 e0 = e0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              ig = iiig
!
           endif
        endif

!
!-- downward leakage sample
     elseif(r1>=PA+PL+PR.and.r1<PA+PL+PR+PD) then

!-- checking if at outer bound
        if (iy == 1) then
           isvacant = .true.
           prt_done = .true.
           tot_eright = tot_eright+e
!-- outbound luminosity tally
!-- sampling x, y
           r1 = rand()
           x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
           y = grd_yarr(iy)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           mu = -max(r1,r2)
           r1 = rand()
           om = pc_pi2*r1
           if(grd_isvelocity) then
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              elabfact = 1d0+dirdotu*cinv
           else
              elabfact = 1d0
           endif
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakdlump
              do iig = 1, glump
                 iiig=glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                      dy(iy)*thelp
                 ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleakd = 0.5d0*ppl/(thelp*dy(iy))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakd*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakd*speclump*help
              enddo
           else
              iiig = ig
           endif
!-- sampling outbound wavlength
           r1 = rand()
           wl=1d0/(r1/grd_wl(iiig+1) + (1d0-r1)/grd_wl(iiig))
!-- changing from comoving frame to observer frame
           if(grd_isvelocity) then
              mu = (mu+y*cinv)/elabfact
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
              wl = wl/elabfact
           endif
!-- obtaining spectrum (lab) group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           iiig = binsrch(wl,flx_wl,flx_ng+1,0)
           if(iiig>flx_ng.or.iiig<1) then
              if(iiig>flx_ng) then
                 iiig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 iiig=1
                 wl=flx_wl(1)
              endif
           endif
           flx_luminos(iiig,imu,1)=flx_luminos(iiig,imu,1)+&
                e*dtinv*elabfact
           flx_lumdev(iiig,imu,1)=flx_lumdev(iiig,imu,1)+&
                (e*dtinv*elabfact)**2
           flx_lumnum(iiig,imu,1)=flx_lumnum(iiig,imu,1)+1

        else

!-- sample adjacent group (assumes aligned ig bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakdlump
              do iig= 1, glump
                 iiig = glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 if((grd_cap(iiig,ix,iy-1,1)+grd_sig(ix,iy-1,1)) * &
                      min(dx(ix),dy(iy-1))*thelp<prt_tauddmc) then
!-- IMC interface
                    mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                         dy(iy)*thelp
                    ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakd = 0.5d0*ppl/(thelp*dy(iy))
                 else
!-- DDMC interface
                    mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                         dy(iy)+&
                         (grd_sig(ix,iy-1,1)+grd_cap(iiig,ix,iy-1,1)) * &
                         dy(iy-1))*thelp
                    resopacleakd = (2d0/3d0)/(mfphelp*thelp*dy(iy))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakd*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakd*speclump*help
              enddo
           else
              iiig = ig
           endif

!-- updating particle properties
           if((grd_sig(ix,iy-1,1)+grd_cap(iiig,ix,iy-1,1)) * &
                min(dx(ix),dy(iy-1)) &
                *thelp >= prt_tauddmc) then
              iy = iy-1
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
              ig = iiig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
!
!-- method changed to IMC
              ptcl%itype = 1
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
              y = grd_yarr(iy)

!-- current particle cell set to 1 left
              iy = iy-1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu = -max(r1,r2)
              r1 = rand()
              om = pc_pi2*r1

!-- doppler and aberration corrections
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 e = e*(1d0+dirdotu*cinv)
                 e0 = e0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              ig = iiig
!
           endif
        endif

!
!-- upward leakage sample
     elseif(r1>=PA+PL+PR+PD.and.r1<PA+PL+PR+PD+PU) then

!-- checking if at outer bound
        if (iy == grd_ny) then
           isvacant = .true.
           prt_done = .true.
           tot_eright = tot_eright+e
!-- outbound luminosity tally
!-- sampling x, y
           r1 = rand()
           x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
           y = grd_yarr(iy+1)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           mu = max(r1,r2)
           r1 = rand()
           om = pc_pi2*r1
           if(grd_isvelocity) then
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              elabfact = 1d0+dirdotu*cinv
           else
              elabfact = 1d0
           endif
!-- determining outbound group
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakulump
              do iig = 1, glump
                 iiig=glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                      dy(iy)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleaku = 0.5d0*ppr/(thelp*dy(iy))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleaku*speclump*help)) exit
                 denom2 = denom2+specig*resopacleaku*speclump*help
              enddo
           else
              iiig = ig
           endif
           r1 = rand()
           wl=1d0/(r1/grd_wl(iiig+1) + (1d0-r1)/grd_wl(iiig))
!-- changing from comoving frame to observer frame
           if(grd_isvelocity) then
              mu = (mu+y*cinv)/elabfact
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
              wl = wl/elabfact
           endif
!-- obtaining spectrum (lab) group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           iiig = binsrch(wl,flx_wl,flx_ng+1,0)
           if(iiig>flx_ng.or.iiig<1) then
              if(iiig>flx_ng) then
                 iiig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 iiig=1
                 wl=flx_wl(1)
              endif
           endif
           flx_luminos(iiig,imu,1)=flx_luminos(iiig,imu,1)+&
                e*dtinv*elabfact
           flx_lumdev(iiig,imu,1)=flx_lumdev(iiig,imu,1)+&
                (e*dtinv*elabfact)**2
           flx_lumnum(iiig,imu,1)=flx_lumnum(iiig,imu,1)+1

        else

!-- sample adjacent group (assumes aligned ig bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakulump
              do iig= 1, glump
                 iiig = glumps(iig)
                 specig = grd_capgrey(ix,iy,1)*grd_emitprob(iiig,ix,iy,1) * &
                      capinv(iiig)
!-- calculating resolved leakage opacities
                 if((grd_cap(iiig,ix,iy+1,1)+grd_sig(ix,iy+1,1)) * &
                      min(dx(ix),dy(iy+1))*thelp<prt_tauddmc) then
!-- IMC interface
                    mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                         dy(iy)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleaku = 0.5d0*ppr/(thelp*dy(iy))
                 else
!-- DDMC interface
                    mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                         dy(iy)+&
                         (grd_sig(ix,iy+1,1)+grd_cap(iiig,ix,iy+1,1)) * &
                         dy(iy+1))*thelp
                    resopacleaku = (2d0/3d0)/(mfphelp*thelp*dy(iy))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleaku*speclump*help)) exit
                 denom2 = denom2+specig*resopacleaku*speclump*help
              enddo
           else
              iiig = ig
           endif

!-- updating particle properties
           if((grd_sig(ix,iy+1,1)+grd_cap(iiig,ix,iy+1,1)) * &
                min(dx(ix),dy(iy+1)) &
                *thelp >= prt_tauddmc) then
              iy = iy+1
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
              ig = iiig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/grd_wl(iiig+1)+(1d0-r1)/grd_wl(iiig))
!
!-- method changed to IMC
              ptcl%itype = 1
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
              y = grd_yarr(iy+1)

!-- current particle cell set to 1 left
              iy = iy+1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu = max(r1,r2)
              r1 = rand()
              om = pc_pi2*r1

!-- doppler and aberration corrections
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 e = e*(1d0+dirdotu*cinv)
                 e0 = e0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              ig = iiig
!
           endif
        endif

!-- effective scattering
     else
        denom2 = 1d0-emitlump
        help = 1d0/denom2
!
        denom3 = 0d0
        r1 = rand()

        do iig = grd_ng,glump+1,-1
           iiig=glumps(iig)
           if((r1>=denom3).and.(r1<denom3+grd_emitprob(iiig,ix,iy,1)*help)) exit
           denom3 = denom3+grd_emitprob(iiig,ix,iy,1)*help
        enddo
!
        ig = iiig
        r1 = rand()
        wl = 1d0/((1d0-r1)/grd_wl(ig) + r1/grd_wl(ig+1))

        if ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1)) * &
             min(dx(ix),dy(iy)) &
             *thelp >= prt_tauddmc) then
           ptcl%itype = 2
        else
           ptcl%itype = 1
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- direction sampled isotropically           
           r1 = rand()
           mu = 1d0 - 2d0*r1
           r1 = rand()
           om = pc_pi2*r1
!-- position sampled uniformly
           r1 = rand()
           x = sqrt(r1*grd_xarr(ix+1)**2+(1d0-r1)*grd_xarr(ix)**2)
           r1 = rand()
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)

!-- doppler and aberration corrections
           if(grd_isvelocity) then
!-- calculating transformation factors
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              om = atan2(sqrt(1d0-mu**2)*sin(om), &
                   sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
              mu = (mu+y/pc_c)/(1d0+dirdotu/pc_c)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming azimuthal angle to lab
              if(om<0d0) om=om+pc_pi2
!-- transforming wavelength to lab
              wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
              e = e*(1d0+dirdotu*cinv)
              e0 = e0*(1d0+dirdotu*cinv)
           endif
        endif
     endif
  endif!}}}


end subroutine diffusion2
