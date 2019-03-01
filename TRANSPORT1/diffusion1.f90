!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine diffusion1(ptcl,ptcl2,cache,rndstate,edep,eraddens,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(grp_t_cache),target,intent(inout) :: cache
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens
  real*8,intent(inout) :: totevelo
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
  logical :: lredir
!
  integer :: iig, iiig, iznext
  logical :: lhelp
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: elabfact, xi, eta
!-- lumped quantities
  real*8 :: emitlump, caplump, doplump
  real*8 :: specig
  real*8 :: opacleak(6)
  real*8 :: probleak(6) !leakage probabilities
  real*8 :: pa, pdop !absorption, doppler probability
  real*8 :: mfphelp, pp
  real*8 :: resopacleak, resdopleak
  integer :: glump, gunlump
  integer*2,pointer :: glumps(:)
  logical*2,pointer :: llumps(:)
  real*8,pointer :: capgreyinv
  real*8,pointer :: speclump
  real*8 :: dist, help, mux, muy, muz
!
  integer,pointer :: ix,iy,iz,ic,ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym,dx2,dx3
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l) = grd_xarr(l+1)**2 - grd_xarr(l)**2
  dx3(l) = grd_xarr(l+1)**3 - grd_xarr(l)**3
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

  capgreyinv => cache%capgreyinv
  speclump => cache%speclump
  glumps => cache%glumps
  llumps => cache%llumps

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0

!-- direction resample flag
  lredir = .false.

!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
  dist = min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * thelp

!
!-- update cache
  if(ic/=cache%ic) then
     cache%ic = ic!{{{
     cache%istat = 0 !specarr is not cached yet
     capgreyinv = max(1d0/grd_capgrey(ic),0d0) !catch nans

!
!-- opacity regrouping --------------------------
     glump = 0
     gunlump = grp_ng
     glumps = 0
!
!-- find lumpable groups
     speclump = grd_opaclump(7,ic)
     if(speclump==0d0) then
        glump=0
     else
        do iig=1,grp_ng
           if(grd_cap(iig,ic)*dist >= trn_taulump .and. &
                (grd_sig(ic) + grd_cap(iig,ic))*dist >= trn_tauddmc) then
              llumps(iig) = .true.
              glump=glump+1
              glumps(glump) = int(iig,2)
           else
              llumps(iig) = .false.
              glumps(gunlump) = int(iig,2)
              gunlump=gunlump-1
           endif
        enddo
     endif
!
!-- calculate lumped values
     if(glump==grp_ng) then
        emitlump = 1d0
        caplump = grd_capgrey(ic)
        doplump = 0d0
     else
!-- Planck x-section lump
        caplump = grd_opaclump(8,ic)*speclump
        emitlump = grd_opaclump(8,ic)*capgreyinv
        emitlump = min(emitlump,1d0)
        doplump = grd_opaclump(10,ic)*speclump
     endif
!
!-- save
     cache%nlump = glump
     cache%emitlump = emitlump
     cache%caplump = caplump
     cache%doplump = doplump
!}}}
  endif !cache%ic /= ic

!
!-- in lump?
  if(grd_cap(ig,ic)*dist >= trn_taulump) then
     glump = cache%nlump
  else
     glump = 0
  endif
!
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic))*dist < trn_tauddmc) then
     ierr = 100
     return
  endif
!
!-- retrieve from cache
  if(glump>0) then
     emitlump = cache%emitlump
     caplump = cache%caplump
     doplump = cache%doplump
  else
!-- outside the lump
     emitlump = specint0(grd_tempinv(ic),ig)*capgreyinv*grd_cap(ig,ic)
     caplump = grd_cap(ig,ic)
     if(grd_isvelocity .and. ig/=grp_ng) then
        doplump = dopspeccalc(grd_tempinv(ic),ig)/(cache%specarr(ig)*pc_c*tsp_t)
     else
        doplump = 0d0
     endif
  endif
!
!-- calculate lumped values
  if(glump>0) then
!-- leakage opacities
     opacleak = grd_opaclump(1:6,ic)
  else
!!{{{
!-- calculating unlumped values
!-- ix->ix-1 in (opacleak(1))
     if(ix==1) then
        lhelp = .true.
     else
        l = grd_icell(ix-1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix-1),xm(ix-1)*dyac(iy) , &
           xm(ix-1)*ym(iy)*dz(iz))*thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)=1.5d0*pp*grd_xarr(ix)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))*thelp
        opacleak(1)=2d0*grd_xarr(ix)**2/(dx3(ix)*thelp*mfphelp)
     endif

!
!-- ix->ix+1 (opacleak(2))
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+grd_sig(l)) * &
             min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy)*dz(iz)) * &
             thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=1.5d0*pp*grd_xarr(ix+1)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))*thelp
        opacleak(2)=2d0*grd_xarr(ix+1)**2/(dx3(ix)*thelp*mfphelp)
     endif

!
!-- iy->iy-1 (opacleak(3))
     if(iy==1) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy-1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy-1), &
             xm(ix)*ym(iy-1)*dz(iz))*thelp<trn_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(3)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dyac(iy) + &
             (grd_sig(l)+grd_cap(ig,l))*dyac(iy-1))
        opacleak(3)=2d0*sqrt(1d0-grd_yarr(iy)**2)*dx(ix) / &
             (dy(iy)*dx3(ix)*mfphelp*thelp**2)
     endif

!
!-- iy->iy+1 (opacleak(4))
     if(iy==grd_ny) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy+1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy+1), &
             xm(ix)*ym(iy+1)*dz(iz))*thelp<trn_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(4)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy+1)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dyac(iy) + &
             (grd_sig(l)+grd_cap(ig,l))*dyac(iy+1))
        opacleak(4)=2d0*sqrt(1d0-grd_yarr(iy+1)**2)*dx(ix) / &
             (dy(iy)*dx3(ix)*mfphelp*thelp**2)
     endif

!
!-- azimuthal leakage opacities
     if(grd_nz==1) then
        opacleak(5:)=0d0
        iznext = iz
     else
!-- iz->iz-1 (opacleak(5))
        if(iz==1) then
           l = grd_icell(ix,iy,grd_nz)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(grd_nz))*thelp<trn_tauddmc
           iznext = grd_nz
        else
           l = grd_icell(ix,iy,iz-1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz-1))*thelp<trn_tauddmc
           iznext = iz-1
        endif
!
        if(lhelp) then
!-- DDMC interface
           mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
           opacleak(5)=0.75d0*pp*dx2(ix)*dyac(iy) / &
                (dy(iy)*dx3(ix)*dz(iz)*thelp)
        else
!-- DDMC interior
           mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(5)=2d0*dyac(iy)*dx(ix) / &
                (ym(iy)*dy(iy)*dz(iz)*dx3(ix)*mfphelp*thelp**2)
        endif
!-- iz->iz+1 (opacleak(6))
        if(iz==grd_nz) then
           l = grd_icell(ix,iy,1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(1))*thelp<trn_tauddmc
           iznext = 1
        else
           l = grd_icell(ix,iy,iz+1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz+1))*thelp<trn_tauddmc
           iznext = iz+1
        endif
!
        if(lhelp) then
!-- DDMC interface
           mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
           opacleak(6)=0.75d0*pp*dx2(ix)*dyac(iy) / &
                (dy(iy)*dx3(ix)*dz(iz)*thelp)
        else
!-- DDMC interior
           mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(6)=2d0*dyac(iy)*dx(ix) / &
                (ym(iy)*dy(iy)*dz(iz)*dx3(ix)*mfphelp*thelp**2)
        endif
     endif!}}}
  endif
!
!--------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
       (1d0-emitlump)*(1d0-grd_fcoef(ic))*caplump + &
       doplump
  if(trn_isddmcanlog) then
     denom = denom+grd_fcoef(ic)*caplump
  endif
  denom = 1d0/denom

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)*denom*cinv)
  tcensus = tsp_t1-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(trn_isddmcanlog) then
     eraddens = e*ddmct*tsp_dtinv
  else
     edep = e * &
          (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*min(dx(ix),xm(ix)*dyac(iy) , &
          xm(ix)*ym(iy)*dz(iz))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = &
             e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*tsp_dtinv
     else
        eraddens = e*ddmct*tsp_dtinv
     endif
!
     if(edep/=edep) then
!       write(0,*) e,grd_fcoef(ic),caplump,ddmct,glump,speclump,ig,grd_tempinv(ic)
!       stop 'diffusion1: invalid energy deposition'
        ierr = 101
        return
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
  if (tcensus < tau) then
     ptcl2%stat = 'cens'
     return
  endif

!-- otherwise, perform event
  call rnd_r(r1,rndstate)

!-- leakage probabilities
  probleak = opacleak*denom
!  write(*,*) probleak
!-- absorption probability
  if(trn_isddmcanlog) then
     pa = grd_fcoef(ic)*caplump*denom
  else
     pa = 0d0
  endif

!-- redshift
  if(grd_isvelocity) then
     pdop = doplump*denom
  else
     pdop = 0d0
  endif

!-- update specarr cache only when necessary. this is slow
  if(r1>=pa .and. r1<pa+sum(probleak(1:6)) .and. speclump>0d0 .and. &
        iand(cache%istat,2)==0) then
     cache%istat = cache%istat + 2
     call specintv(grd_tempinv(ic),grp_ng,cache%specarr,&
        mask=llumps,maskval=.true.)
  endif

!-- absorption
  if(r1<pa) then
     ptcl2%idist = -1
     ptcl2%stat = 'dead'
     edep = e

!-- doppler shift
  elseif (r1>=pa .and. r1<pa+pdop) then

     if(glump==0) then
        iiig = ig
     else
!-- sample group
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/doplump
        do iig=1,glump
           iiig = glumps(iig)
           if(iiig == grp_ng) cycle
           if(grd_cap(iiig+1,ic)*dist >= trn_taulump) cycle
           specig = cache%specarr(iiig)
           resdopleak = dopspeccalc(grd_tempinv(ic),iiig)/(pc_c*tsp_t)
           denom2 = denom2+resdopleak*speclump*help
           if(denom2>r1) exit
        enddo
     endif

!-- reshift particle in this group
     ig = iiig+1
     wl = grp_wl(ig)
     ig = min(ig,grp_ng)

!-- method changes to IMC
     if((grd_sig(ic)+grd_cap(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- direction sampled isotropically
        lredir = .true.
        call rnd_r(r1,rndstate)
        mu = 1d0 - 2d0*r1
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(ix+1)**3+(1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
        call rnd_r(r1,rndstate)
        y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- velocity effects accounting
        mu = (mu+x*cinv)/(1.0+x*mu*cinv)
        wl = wl*(1.0-x*mu*cinv)
        help = 1d0/(1.0-x*mu*cinv)
        totevelo = totevelo+e*(1d0 - help)
        e = e*help
        e0 = e0*help
     endif

!-- ix->ix-1 leakage
  elseif(r1>=pa+pdop.and.r1<pa+pdop+probleak(1)) then
     ptcl2%idist = -3
!-- sanity check!{{{
     if(ix==1) then
!       write(0,*) stop 'diffusion1: invalid probleak(1)'
        ierr = 102
        return
     endif
!
     l = grd_icell(ix-1,iy,iz)
!-- sample next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(1)
        do iig=1,glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix-1),xm(ix-1)*dyac(iy),xm(ix-1)*ym(iy) * &
                   dz(iz))*thelp<trn_tauddmc
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*pp*grd_xarr(ix)**2/(dx3(ix)*thelp)
           else
!-- IMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dx(ix)+&
                   (grd_sig(l)+grd_cap(iiig,l))*dx(ix-1))*thelp
              resopacleak = 2d0*grd_xarr(ix)**2/(mfphelp*thelp*dx3(ix))
           endif
           denom2 = denom2+specig*resopacleak*speclump*help
           if(denom2>r1) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix-1),xm(ix-1)*dyac(iy),xm(ix-1)*ym(iy)*dz(iz)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        x = grd_xarr(ix)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = -max(r1,r2)
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
        if(grd_isvelocity) then
           elabfact = 1d0+x*mu*cinv
        else
           elabfact = 1d0
        endif
!-- transforming mu to lab
        if(grd_isvelocity) then
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-x*mu*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
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
     ig = iiig!}}}

!-- ix->ix+1 leakage
  elseif(r1>=pa+pdop+probleak(1).and.r1<pa+pdop+sum(probleak(1:2))) then
     ptcl2%idist = -4
!{{{
!-- sampling next group
     if(ix/=grd_nx) l = grd_icell(ix+1,iy,iz)
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(2)
        do iig=1,glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(ix==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy) * &
                   dz(iz))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*pp*grd_xarr(ix+1)**2/(dx3(ix)*thelp)
           else
!-- IMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dx(ix)+&
                   (grd_sig(l)+grd_cap(iiig,l))*dx(ix+1))*thelp
              resopacleak = 2d0*grd_xarr(ix+1)**2/(mfphelp*thelp*dx3(ix))
           endif
           denom2 = denom2+specig*resopacleak*speclump*help
           if(denom2>r1) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy)*dz(iz)) * &
             thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        x = grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
        if(grd_isvelocity) then
           elabfact = 1d0+x*mu*cinv
        else
           elabfact = 1d0
        endif
!-- transforming mu to lab
        if(grd_isvelocity) then
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-x*mu*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if(ix==grd_nx) then
!-- escaping at ix=nx
           ptcl2%stat = 'flux'
!-- observer time correction
           ptcl%t=ptcl%t-mu*x*thelp*cinv
!-- luminosity tally
           eta = sqrt(1d0-mu**2)*cos(om)
           xi = sqrt(1d0-mu**2)*sin(om)
           mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
           muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
           muz = mu*y-eta*sqrt(1d0-y**2)
!-- redifine for flux tally
           mu = muz
           om = atan2(muy,mux)
           if(om<0d0) om=om+pc_pi2
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
     ig = iiig!}}}

!-- iy->iy-1 leakage
  elseif(r1>=pa+pdop+sum(probleak(1:2)).and.r1<pa+pdop+sum(probleak(1:3))) then
     ptcl2%idist = -5
!-- sanity check!{{{
     if(grd_ny==1) then
!       stop 'diffusion1: probleak(3) and grd_ny=1'
        ierr = 103
        return
     endif
     if(iy==1) then
!       stop 'diffusion1: invalid probleak(3)'
        ierr = 104
        return
     endif

!-- sampling next group
     l = grd_icell(ix,iy-1,iz)
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(3)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy-1),xm(ix)*ym(iy-1)*dz(iz)) * &
                thelp<trn_tauddmc
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*dyac(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy)**2) / &
                   (dy(iy)*dx3(ix)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dyac(iy) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dyac(iy-1))
              resopacleak=2d0*sqrt(1d0-grd_yarr(iy)**2)*dx(ix) / &
                   (dy(iy)*dx3(ix)*mfphelp*thelp**2)
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif
!
!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!
!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy-1),xm(ix)*ym(iy-1)*dz(iz)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        y = grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        eta = max(r1,r2)
        call rnd_r(r1,rndstate)
        xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
        om = atan2(xi,eta)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact=1d0+mu*x*cinv
        else
           elabfact=1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-x*mu*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
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
!-- update particle: iy->iy-1
     iy = iy-1
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- iy->iy+1 leakage
  elseif(r1>=pa+pdop+sum(probleak(1:3)).and.r1<pa+pdop+sum(probleak(1:4))) then
     ptcl2%idist = -6
!-- sanity check!{{{
     if(grd_ny==1) then
!       stop 'diffusion1: probleak(4) and grd_ny=1'
        ierr = 105
        return
     endif
     if(iy==grd_ny) then
!       stop 'diffusion1: invalid probleak(4)'
        ierr = 106
        return
     endif

!-- sampling next group
     l = grd_icell(ix,iy+1,iz)
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(4)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = cache%specarr(iiig)
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy+1),xm(ix)*ym(iy+1)*dz(iz)) * &
                thelp<trn_tauddmc
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*dyac(iy)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy+1)**2) / &
                   (dy(iy)*dx3(ix)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dyac(iy) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dyac(iy+1))
              resopacleak=2d0*sqrt(1d0-grd_yarr(iy+1)**2)*dx(ix) / &
                   (dy(iy)*dx3(ix)*mfphelp*thelp**2)
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
          min(dx(ix),xm(ix)*dyac(iy+1),xm(ix)*ym(iy+1)*dz(iz)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        y = grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        eta = -max(r1,r2)
        call rnd_r(r1,rndstate)
        xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
        om = atan2(xi,eta)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact = 1d0+mu*x*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-mu*x*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
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
!-- update particle: iy->iy+1
     iy = iy+1
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- iz->iz-1 leakage
  elseif(r1>=pa+pdop+sum(probleak(1:4)).and.r1<pa+pdop+sum(probleak(1:5))) then
     ptcl2%idist = -7
!-- sanity check!{{{
     if(grd_nz==1) then
!       stop 'diffusion1: probleak(5) and grd_nz=1'
        ierr = 107
        return
     endif
!-- setting helper index
     if(iz==1) then
        iznext=grd_nz
     else
        iznext=iz-1
     endif
     l = grd_icell(ix,iy,iznext)
!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(5)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
                thelp<trn_tauddmc
!
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*ym(iy)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*dyac(iy) / &
                   (dy(iy)*dx3(ix)*dz(iz)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dz(iz) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dz(iznext))
              resopacleak=2d0*dyac(iy)*dx(ix)/(ym(iy)*dy(iy)*dz(iz)*dx3(ix) * &
                   mfphelp*thelp**2)
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        if(iznext==iz-1) then
           z = grd_zarr(iz)
        else
!-- iz = 1, iznext = grd_nz
           z = pc_pi2
        endif
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        xi = -max(r1,r2)
        call rnd_r(r1,rndstate)
        eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
        om = atan2(xi,eta)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact=1d0+mu*x*cinv
        else
           elabfact=1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-x*mu*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
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
!-- update particle: iz->iz-1
     iz = iznext
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- iz->iz+1 leakage
  elseif(r1>=pa+pdop+sum(probleak(1:5)).and.r1<pa+pdop+sum(probleak(1:6))) then
     ptcl2%idist = -8
!-- sanity check!{{{
     if(grd_nz==1) then
!       stop 'diffusion1: probleak(6) and grd_nz=1'
        ierr = 108
        return
     endif
!-- setting helper index
     if(iz==grd_nz) then
        iznext=1
     else
        iznext=iz+1
     endif
     l = grd_icell(ix,iy,iznext)
!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(6)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
                thelp<trn_tauddmc
!
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*ym(iy)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*dyac(iy) / &
                   (dy(iy)*dx3(ix)*dz(iz)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dz(iz) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dz(iznext))
              resopacleak=2d0*dyac(iy)*dx(ix)/(ym(iy)*dy(iy)*dz(iz)*dx3(ix) * &
                   mfphelp*thelp**2)
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        if(iznext==iz+1) then
           z = grd_zarr(iz+1)
        else
!-- iz = grd_nz, iznext = 1
           z = 0d0
        endif
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- sampling direction
        lredir = .true.
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        xi = max(r1,r2)
        call rnd_r(r1,rndstate)
        eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
        om = atan2(xi,eta)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact=1d0+mu*x*cinv
        else
           elabfact=1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-x*mu*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
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
!-- update particle: iz->iz+1
     iz = iznext
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- effective scattering
  else
     ptcl2%idist = -2
!!{{{
     if(glump==grp_ng) then
!       stop 'diffusion1: effective scattering with glump==ng'
        ierr = 109
        return
     endif

     if(glump==0) then
!-- sample group
        if(cache%emitlump<.99d0 .or. trn_nolumpshortcut) then
           call rnd_r(r1,rndstate)
           iiig = emitgroup(r1,ic)
!-- don't sample, it will end up in the lump anyway
        else
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(9,ic))
           return
        endif
     else
!
!-- update specarr cache. this is slow
        if(iand(cache%istat,1)==0) then
           cache%istat = cache%istat + 1
           call specintv(grd_tempinv(ic),grp_ng,cache%specarr,&
              mask=llumps,maskval=.false.)
        endif
!
        call rnd_r(r1,rndstate)
        denom2 = 1d0/(1d0-emitlump)
        denom3 = 0d0
        do iig=grp_ng,glump+1,-1
           iiig = glumps(iig)
           help = cache%specarr(iiig)*grd_cap(iiig,ic)*capgreyinv
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
     ig = iiig

     if((grd_sig(ic)+grd_cap(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
!-- direction sampled isotropically
        lredir = .true.
        call rnd_r(r1,rndstate)
        mu = 1d0 - 2d0*r1
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(ix+1)**3+(1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
        call rnd_r(r1,rndstate)
        y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!-- doppler and aberration corrections
        if(grd_isvelocity) then
!-- calculating transformation factors
           elabfact = 1d0+mu*x*cinv
           mu = (mu+x*cinv)/elabfact
!-- ELABFACT LAB RESET
           elabfact=1d0-mu*x*cinv
           help = 1d0/elabfact
!-- transforming wl to lab
           wl = wl*elabfact
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
     endif
!}}}
  endif

  if(lredir) then
!-- spherical projections
     eta = sqrt(1d0-mu**2)*cos(om)
     xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
     ptcl2%mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
     ptcl2%muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
     ptcl2%muz = mu*y-eta*sqrt(1d0-y**2)
  endif

end subroutine diffusion1
! vim: fdm=marker
