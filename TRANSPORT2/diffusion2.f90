subroutine diffusion2(ptcl,isvacant,icell,specarr)

  use gridmod
  use groupmod
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
  integer,intent(inout) :: icell(3)
  real*8,intent(inout) :: specarr(grp_ng)
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer,external :: binsrch, emitgroup
!
  integer :: ig, iig, iiig, imu
  logical :: lhelp
  real*8 :: r1, r2, thelp, mu0
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: dirdotu, gm
!-- lumped quantities
  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleak(4)
  real*8 :: probleak(4) !leakage probabilities
  real*8 :: pa !absorption probability
  real*8 :: mfphelp, pp
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer :: glumps(grp_ng)
  real*8 :: dtinv, tempinv, capgreyinv
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
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  tempinv = 1d0/grd_temp(ix,iy,1)
  capgreyinv = max(1d0/grd_capgrey(ix,iy,1),0d0) !catch nans

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- looking up initial group
  ig = binsrch(wl,grp_wl,grp_ng+1,in_ng)
!-- checking group bounds
  if(ig>grp_ng.or.ig<1) then
     if(ig==grp_ng+1) then
        ig = grp_ng
     elseif(ig==0) then
        ig = 1
     else
        stop 'diffusion2: particle group invalid'
     endif
  endif

!
!-- opacity regrouping --------------------------
  glump = 0
  gunlump = grp_ng
  glumps = 0
!
!-- find lumpable groups
  if(grd_cap(ig,ix,iy,1)*min(dx(ix),dy(iy)) * thelp>=prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ix,iy,1)*min(dx(ix),dy(iy))*thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
  endif
! write(0,*) ipart,istep,glump,g,ix,iy

!
!-- only do this if needed
  if(glump>0 .and. .not.all(icell==[ix,iy,1])) then
     icell = [ix,iy,1]
     specarr = specintv(tempinv) !this is slow!
  endif

!
  if(glump==0) then
     forall(iig=1:grp_ng) glumps(iig)=iig
  endif

!
!-- lumping
  speclump = 0d0
  do iig=1,glump
     iiig = glumps(iig)
     specig = specarr(iiig)
     speclump = speclump + specig
  enddo
  if(speclump>0d0) then
     speclump = 1d0/speclump
  else
     speclump = 0d0
  endif

!write(0,*) impi,glump,speclump
!
  emitlump = 0d0
  caplump = 0d0
!-- calculate lumped values
  if(speclump>0d0) then
     if(glump==grp_ng) then!{{{
        emitlump = 1d0
        caplump = grd_capgrey(ix,iy,1)
     else
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- emission lump
           emitlump = emitlump + specig*capgreyinv*grd_cap(iiig,ix,iy,1)
!-- Planck x-section lump
           caplump = caplump + specig*grd_cap(iiig,ix,iy,1)*speclump
        enddo
        emitlump = min(emitlump,1d0)
     endif
!-- leakage opacities
     opacleak = grd_opacleak(:4,ix,iy,1)
!!}}}
  else
!
!-- calculating unlumped values
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ix,iy,1)!{{{
     caplump = grd_cap(ig,ix,iy,1)

!-- inward
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,ix-1,iy,1)+ &
          grd_sig(ix-1,iy,1))*min(dx(ix-1),dy(iy))* &
          thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1))*dx(ix)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= pp*(thelp*grd_xarr(ix))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dx(ix)+&
             (grd_sig(ix-1,iy,1)+grd_cap(ig,ix-1,iy,1))*dx(ix-1))*thelp
        opacleak(1)=(4d0/3d0)*(thelp*grd_xarr(ix))/ &
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
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)= pp*(thelp*grd_xarr(ix+1))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dx(ix)+&
             (grd_sig(ix+1,iy,1)+grd_cap(ig,ix+1,iy,1))*dx(ix+1))*thelp
        opacleak(2)=(4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
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
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(3)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dy(iy)+&
             (grd_sig(ix,iy-1,1)+grd_cap(ig,ix,iy-1,1))*dy(iy-1))*thelp
        opacleak(3)=(2d0/3d0)/(help*dy(iy)*thelp)
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
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(4)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1))*dy(iy)+&
             (grd_sig(ix,iy+1,1)+grd_cap(ig,ix,iy+1,1))*dy(iy+1))*thelp
        opacleak(4)=(2d0/3d0)/(help*dy(iy)*thelp)
     endif
!}}}
  endif
!
!--------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
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
  if(prt_isddmcanlog) then
     grd_eraddens(ix,iy,1) = grd_eraddens(ix,iy,1)+e*ddmct*dtinv
  else
     grd_edep(ix,iy,1) = grd_edep(ix,iy,1)+e*(1d0-exp(-grd_fcoef(ix,iy,1) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ix,iy,1)*caplump*min(dx(ix),dy(iy))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ix,iy,1)*caplump)
        grd_eraddens(ix,iy,1)= &
             grd_eraddens(ix,iy,1)+e* &
             (1d0-exp(-grd_fcoef(ix,iy,1)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        grd_eraddens(ix,iy,1) = grd_eraddens(ix,iy,1)+e*ddmct*dtinv
     endif
     e = e*exp(-grd_fcoef(ix,iy,1)*caplump*pc_c*ddmct)
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
     return
  endif


!-- otherwise, perform event
  r1 = rand()
  help = 1d0/denom

!-- leakage probabilities
  probleak = opacleak*help

!-- absorption probability
  if(prt_isddmcanlog) then
     pa = grd_fcoef(ix,iy,1)*caplump*help
  else
     pa = 0d0
  endif

!-- absorption
  if(r1<pa) then
     isvacant = .true.
     prt_done = .true.
     grd_edep(ix,iy,1) = grd_edep(ix,iy,1)+e

!-- ix->ix-1 leakage
  elseif(r1>=pa.and.r1<pa+probleak(1)) then
!{{{
!-- sanity check
     if (ix == 1) stop 'diffusion1: non-physical inward leakage'

!-- sample next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rand()
        denom2 = 0d0
        help = 1d0/opacleak(1)
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if((grd_cap(iiig,ix-1,iy,1) + &
                grd_sig(ix-1,iy,1))*min(dx(ix-1),dy(iy)) * &
                thelp<prt_tauddmc) then
!-- IMC interface
              mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = pp*(thelp*grd_xarr(ix))/ &
                   (thelp**2*dx2(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                   dx(ix)+(grd_sig(ix-1,iy,1) + &
                   grd_cap(iiig,ix-1,iy,1))*dx(ix-1))*thelp
              resopacleak = (4d0/3d0)*(thelp*grd_xarr(ix))/ &
                   (mfphelp*thelp**2*dx2(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     r1 = rand()
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,ix-1,iy,1)+grd_sig(ix-1,iy,1)) * &
                min(dx(ix-1),dy(iy))*thelp<prt_tauddmc

     if(.not.lhelp) then
!-- ix->ix-1
        ix = ix-1
     else
!-- sampling x,y
        x = grd_xarr(ix)
        r1 = rand()
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
!-- sampling direction
        r1 = rand()
        r2 = rand()
        mu0 = -max(r1,r2)
        r1 = rand()
        mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
        om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
        if(om<0d0) om = om+pc_pi2
        if(grd_isvelocity) then
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)+gm*x*cinv * &
                (1d0+gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu+gm*y*cinv*(1d0+gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0+dirdotu*cinv))
!-- DIRDOTU LAB RESET
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           help = 1d0/(1d0-dirdotu*cinv)
!-- transforming wavelength to lab
           wl = wl*(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl%itype = 1
        grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- ix->ix-1
        ix = ix-1
     endif
!}}}

!-- ix->ix+1 leakage
  elseif (r1>=pa+probleak(1).and.r1<pa+sum(probleak(1:2))) then
!{{{
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rand()
        denom2 = 0d0
        help = 1d0/opacleak(2)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(ix==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,ix+1,iy,1)+grd_sig(ix+1,iy,1)) * &
                   min(dx(ix+1),dy(iy))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = pp*(thelp*grd_xarr(ix+1))/ &
                   (thelp**2*dx2(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                   dx(ix)+&
                   (grd_sig(ix+1,iy,1)+grd_cap(iiig,ix+1,iy,1)) * &
                   dx(ix+1))*thelp
              resopacleak = (4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
                   (mfphelp*thelp**2*dx2(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavlength
     r1 = rand()
     wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,ix+1,iy,1)+grd_sig(ix+1,iy,1)) * &
             min(dx(ix+1),dy(iy))*thelp<prt_tauddmc
     endif

     if(.not.lhelp) then
!-- ix->ix+1
        ix = ix+1
     else
!-- sampling x,y
        x = grd_xarr(ix+1)
        r1 = rand()
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
!-- sampling direction
        r1 = rand()
        r2 = rand()
        mu0 = max(r1,r2)
        r1 = rand()
        mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
        om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)+gm*x*cinv * &
                (1d0+gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu+gm*y*cinv*(1d0+gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0+dirdotu*cinv))
!-- DIRDOTU LAB RESET
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           help = 1d0/(1d0-dirdotu*cinv)
!-- transforming wavelength to lab
           wl = wl*(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (ix==grd_nx) then
!-- escaping at ix=nx
           isvacant = .true.
           prt_done = .true.
           tot_eout = tot_eout+e
!-- luminosity tally
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
           flx_luminos(iiig,imu,1) = flx_luminos(iiig,imu,1) + e*dtinv
           flx_lumdev(iiig,imu,1) = flx_lumdev(iiig,imu,1) + (e*dtinv)**2
           flx_lumnum(iiig,imu,1) = flx_lumnum(iiig,imu,1) + 1
        else
!-- converting to IMC
           ptcl%itype = 1
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- ix->ix+1
           ix = ix+1
        endif
     endif
!}}}

!-- iy->iy-1 leakage
  elseif(r1>=pa+sum(probleak(1:2)).and.r1<pa+sum(probleak(1:3))) then
!{{{
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rand()
        denom2 = 0d0
        help = 1d0/opacleak(3)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,ix,iy-1,1)+grd_sig(ix,iy-1,1)) * &
                   min(dx(ix),dy(iy-1))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                   dy(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dy(iy))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                   dy(iy)+&
                   (grd_sig(ix,iy-1,1)+grd_cap(iiig,ix,iy-1,1)) * &
                   dy(iy-1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dy(iy))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavlength
     r1 = rand()
     wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iy==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,ix,iy-1,1)+grd_sig(ix,iy-1,1)) * &
             min(dx(ix),dy(iy-1))*thelp<prt_tauddmc
     endif

     if(.not.lhelp) then
!-- iy->iy-1
        iy = iy-1
     else
!-- sampling x,y
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
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)+gm*x*cinv * &
                (1d0+gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu+gm*y*cinv*(1d0+gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0+dirdotu*cinv))
!-- DIRDOTU LAB RESET
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           help = 1d0/(1d0-dirdotu*cinv)
!-- transforming wl to lab
           wl = wl*(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (iy==1) then
!-- escaping at iy=1
           isvacant = .true.
           prt_done = .true.
           tot_eout = tot_eout+e
!-- luminosity tally
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
           flx_luminos(iiig,imu,1) = flx_luminos(iiig,imu,1) + e*dtinv
           flx_lumdev(iiig,imu,1) = flx_lumdev(iiig,imu,1) + (e*dtinv)**2
           flx_lumnum(iiig,imu,1) = flx_lumnum(iiig,imu,1) + 1
        else
!-- converting to IMC
           ptcl%itype = 1
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- iy->iy-1
           iy = iy-1
        endif
     endif
!}}}

!-- iy->iy+1 leakage
  elseif(r1>=pa+sum(probleak(1:3)).and.r1<pa+sum(probleak(1:4))) then
!{{{
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rand()
        denom2 = 0d0
        help = 1d0/opacleak(4)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==grd_ny) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,ix,iy+1,1)+grd_sig(ix,iy+1,1)) * &
                   min(dx(ix),dy(iy+1))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ix,iy,1)+grd_sig(ix,iy,1)) * &
                   dy(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dy(iy))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ix,iy,1)+grd_cap(iiig,ix,iy,1)) * &
                   dy(iy)+&
                   (grd_sig(ix,iy+1,1)+grd_cap(iiig,ix,iy+1,1)) * &
                   dy(iy+1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dy(iy))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     r1 = rand()
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iy==grd_ny) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,ix,iy+1,1)+grd_sig(ix,iy+1,1)) * &
             min(dx(ix),dy(iy+1))*thelp<prt_tauddmc
     endif

     if(.not.lhelp) then
!-- iy->iy+1
        iy = iy+1
     else
!-- sampling x,y
        r1 = rand()
        x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
        y = grd_yarr(iy+1)
!-- sampling direction
        r1 = rand()
        r2 = rand()
        mu = max(r1,r2)
        r1 = rand()
        om = pc_pi2*r1
!-- doppler and aberration corrections
        if(grd_isvelocity) then
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)+gm*x*cinv * &
                (1d0+gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu+gm*y*cinv*(1d0+gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0+dirdotu*cinv))
!-- DIRDOTU LAB RESET
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           help = 1d0/(1d0-dirdotu*cinv)
!-- transforming wl to lab
           wl = wl*(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (iy == grd_ny) then
!-- escaping at iy=ny
           isvacant = .true.
           prt_done = .true.
           tot_eout = tot_eout+e
!-- luminosity tally
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
           flx_luminos(iiig,imu,1) = flx_luminos(iiig,imu,1) + e*dtinv
           flx_lumdev(iiig,imu,1) = flx_lumdev(iiig,imu,1) + (e*dtinv)**2
           flx_lumnum(iiig,imu,1) = flx_lumnum(iiig,imu,1) + 1
        else
!-- converting to IMC
           ptcl%itype = 1
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- iy->iy+1
           iy = iy+1
        endif
     endif
!}}}

!-- effective scattering
  else
!
     if(glump==grp_ng) stop 'diffusion2: effective scattering with glump==ng'

     r1 = rand()

     if(glump==0) then
        iiig = emitgroup(r1,ix,iy,1)
     else
        denom3 = 0d0
        denom2 = 1d0-emitlump
        denom2 = 1d0/denom2
        do iig = glump+1,grp_ng
           iiig=glumps(iig)
           if(all(icell==[ix,iy,1])) then
              help = specarr(iiig)*grd_cap(iiig,ix,iy,1)*capgreyinv
           else
              help = specint0(tempinv,iiig)*grd_cap(iiig,ix,iy,1)*capgreyinv
           endif
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
!
     ig = iiig
     r1 = rand()
     wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))

     if((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1)) * &
          min(dx(ix),dy(iy)) &
          *thelp < prt_tauddmc) then
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
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)+gm*x*cinv * &
                (1d0+gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu+gm*y*cinv*(1d0+gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0+dirdotu*cinv))
!-- DIRDOTU LAB RESET
           dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
           help = 1d0/(1d0-dirdotu*cinv)
!-- transforming wl to lab
           wl = wl*(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
     endif

  endif

end subroutine diffusion2
