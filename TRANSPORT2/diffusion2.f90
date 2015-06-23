pure subroutine diffusion2(ptcl,ptcl2,rndstate,edep,eraddens,totevelo,icspec,specarr,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens
  real*8,intent(inout) :: totevelo
  integer,intent(inout) :: icspec
  real*8,intent(inout) :: specarr(grp_ng)
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
!
  integer :: iig, iiig
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
  real*8 :: dist, help

  integer,pointer :: ix, iy, ic, ig
  integer,parameter :: iz=1
  real*8,pointer :: x,y,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2-grd_xarr(l)**2
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

  ix => ptcl2%ix
  iy => ptcl2%iy
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x
  y => ptcl%y
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  tempinv = 1d0/grd_temp(ic)
  capgreyinv = max(1d0/grd_capgrey(ic),0d0) !catch nans

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!
!-- opacity regrouping --------------------------
  glump = 0
  gunlump = grp_ng
  glumps = 0
!
!-- find lumpable groups
  dist = min(dx(ix),dy(iy)) * thelp
  if(grd_cap(ig,ic)*dist >= prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ic)*dist >= prt_taulump .and. &
             (grd_sig(ic) + grd_cap(iig,ic))*dist >= prt_tauddmc) then
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
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic))*dist < prt_tauddmc) then
     ierr = 100
     return
  endif

!
!-- only do this if needed
  if(glump>0 .and. icspec/=ic) then
     icspec = ic
     specarr = specintv(tempinv,0) !this is slow!
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
        caplump = grd_capgrey(ic)
     else
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- emission lump
           emitlump = emitlump + specig*capgreyinv*grd_cap(iiig,ic)
!-- Planck x-section lump
           caplump = caplump + specig*grd_cap(iiig,ic)*speclump
        enddo
        emitlump = min(emitlump,1d0)
     endif
!-- leakage opacities
     opacleak = grd_opacleak(:4,ic)
!!}}}
  else
!
!-- calculating unlumped values
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ic)!{{{
     caplump = grd_cap(ig,ic)

!-- inward
     if(ix/=1) l = grd_icell(ix-1,iy,iz)
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,l)+ &
        grd_sig(l))*min(dx(ix-1),dy(iy))* &
        thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= pp*(thelp*grd_xarr(ix))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))*thelp
        opacleak(1)=(4d0/3d0)*(thelp*grd_xarr(ix))/ &
             (mfphelp*thelp**2*dx2(ix))
     endif

!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix+1),dy(iy))* &
           thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)= pp*(thelp*grd_xarr(ix+1))/ &
             (thelp**2*dx2(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))*thelp
        opacleak(2)=(4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
             (mfphelp*thelp**2*dx2(ix))
     endif

!-- downward
     if(iy==1) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy-1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),dy(iy-1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        dist = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        pp = 4d0/(3d0*dist+6d0*pc_dext)
        opacleak(3)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        dist = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy-1))*thelp
        opacleak(3)=(2d0/3d0)/(dist*dy(iy)*thelp)
     endif

!-- upward
     if(iy==grd_ny) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy+1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),dy(iy+1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        dist = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        pp = 4d0/(3d0*dist+6d0*pc_dext)
        opacleak(4)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        dist = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy+1))*thelp
        opacleak(4)=(2d0/3d0)/(dist*dy(iy)*thelp)
     endif
!}}}
  endif
!
!--------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
       (1d0-emitlump)*(1d0-grd_fcoef(ic))*caplump
  if(prt_isddmcanlog) then
     denom = denom+grd_fcoef(ic)*caplump
  endif

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     eraddens = e*ddmct*dtinv
  else
     edep = e*(1d0-exp(-grd_fcoef(ic) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*min(dx(ix),dy(iy))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        eraddens = e*ddmct*dtinv
     endif
     e = e*exp(-grd_fcoef(ic)*caplump*pc_c*ddmct)
!!}}}
  endif


!-- updating particle time
  ptcl%t = ptcl%t+ddmct


!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
     return
  endif


!-- otherwise, perform event
  call rnd_r(r1,rndstate)
  help = 1d0/denom

!-- leakage probabilities
  probleak = opacleak*help

!-- absorption probability
  if(prt_isddmcanlog) then
     pa = grd_fcoef(ic)*caplump*help
  else
     pa = 0d0
  endif

!-- absorption
  if(r1<pa) then
     ptcl2%isvacant = .true.
     ptcl2%done = .true.
     edep = e

!-- ix->ix-1 leakage
  elseif(r1>=pa.and.r1<pa+probleak(1)) then
!{{{
!-- sanity check
     if(ix==1) then
!       stop 'diffusion1: non-physical inward leakage'
        ierr = 101
        return
     endif

     l = grd_icell(ix-1,iy,iz)
!-- sample next group
     if(speclump<=0d0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(1)
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if((grd_cap(iiig,l) + &
                grd_sig(l))*min(dx(ix-1),dy(iy)) * &
                thelp<prt_tauddmc) then
!-- IMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = pp*(thelp*grd_xarr(ix))/ &
                   (thelp**2*dx2(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dx(ix)+(grd_sig(l) + &
                   grd_cap(iiig,l))*dx(ix-1))*thelp
              resopacleak = (4d0/3d0)*(thelp*grd_xarr(ix))/ &
                   (mfphelp*thelp**2*dx2(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix-1),dy(iy))*thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y
        x = grd_xarr(ix)
        call rnd_r(r1,rndstate)
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
!-- must be inside cell
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu0 = -max(r1,r2)
        call rnd_r(r1,rndstate)
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
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
     endif
!
!-- update particle: ix->ix-1
     ix = ix-1
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- ix->ix+1 leakage
  elseif (r1>=pa+probleak(1).and.r1<pa+sum(probleak(1:2))) then
!{{{
     if(ix/=grd_nx) l = grd_icell(ix+1,iy,iz)
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(2)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(ix==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix+1),dy(iy))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = pp*(thelp*grd_xarr(ix+1))/ &
                   (thelp**2*dx2(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dx(ix)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dx(ix+1))*thelp
              resopacleak = (4d0/3d0)*(thelp*grd_xarr(ix+1))/ &
                   (mfphelp*thelp**2*dx2(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavlength
     call rnd_r(r1,rndstate)
     wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix+1),dy(iy))*thelp<prt_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y
        x = grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu0 = max(r1,r2)
        call rnd_r(r1,rndstate)
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
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (ix==grd_nx) then
!-- escaping at ix=nx
           ptcl2%isvacant = .true.
           ptcl2%done = .true.
           ptcl2%lflux = .true.
           return
        else
!-- converting to IMC
           ptcl2%itype = 1
        endif
     endif
!
!-- update particle: ix->ix+1
     ix = ix+1
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- iy->iy-1 leakage
  elseif(r1>=pa+sum(probleak(1:2)).and.r1<pa+sum(probleak(1:3))) then
!{{{
     if(iy/=1) l = grd_icell(ix,iy-1,iz)
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(3)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy-1))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dy(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dy(iy))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dy(iy)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dy(iy-1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dy(iy))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavlength
     call rnd_r(r1,rndstate)
     wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iy==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix),dy(iy-1))*thelp<prt_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y
        call rnd_r(r1,rndstate)
        x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
        y = grd_yarr(iy)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = -max(r1,r2)
        call rnd_r(r1,rndstate)
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
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (iy==1) then
!-- escaping at iy=1
           ptcl2%isvacant = .true.
           ptcl2%done = .true.
           ptcl2%lflux = .true.
           return
        else
!-- converting to IMC
           ptcl2%itype = 1
        endif
     endif
!
!-- update particle: iy->iy-1
     iy = iy-1
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- iy->iy+1 leakage
  elseif(r1>=pa+sum(probleak(1:3)).and.r1<pa+sum(probleak(1:4))) then
!{{{
     if(iy/=grd_ny) l = grd_icell(ix,iy+1,iz)
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(4)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==grd_ny) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy+1))*thelp<prt_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dy(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dy(iy))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dy(iy)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dy(iy+1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dy(iy))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iy==grd_ny) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix),dy(iy+1))*thelp<prt_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y
        call rnd_r(r1,rndstate)
        x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
        y = grd_yarr(iy+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        call rnd_r(r1,rndstate)
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
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if (iy == grd_ny) then
!-- escaping at iy=ny
           ptcl2%isvacant = .true.
           ptcl2%done = .true.
           ptcl2%lflux = .true.
           return
        else
!-- converting to IMC
           ptcl2%itype = 1
        endif
     endif
!
!-- update particle: iy->iy+1
     iy = iy+1
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- effective scattering
  else
!
     if(glump==grp_ng) then
!       stop 'diffusion2: effective scattering with glump==ng'
        ierr = 102
        return
     endif

     call rnd_r(r1,rndstate)

     if(glump==0) then
        iiig = emitgroup(r1,ic)
        if(iiig>grp_ng) then
!          stop 'diffusion2: emitgroup ig>ng'
           ierr = 103
           return
        endif
     else
        denom3 = 0d0
        denom2 = 1d0-emitlump
        denom2 = 1d0/denom2
        do iig = glump+1,grp_ng
           iiig=glumps(iig)
           if(icspec==ic) then
              help = specarr(iiig)*grd_cap(iiig,ic)*capgreyinv
           else
              help = specint0(tempinv,iiig)*grd_cap(iiig,ic)*capgreyinv
           endif
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
!
     ig = iiig
     call rnd_r(r1,rndstate)
     wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))

     if((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),dy(iy)) &
          *thelp < prt_tauddmc) then
        ptcl2%itype = 1
!-- direction sampled isotropically           
        call rnd_r(r1,rndstate)
        mu = 1d0 - 2d0*r1
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = sqrt(r1*grd_xarr(ix+1)**2+(1d0-r1)*grd_xarr(ix)**2)
        call rnd_r(r1,rndstate)
        y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
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
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
     endif

  endif

end subroutine diffusion2
