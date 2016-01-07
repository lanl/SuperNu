pure subroutine diffusion2(ptcl,ptcl2,cache,rndstate,edep,eraddens,totevelo,ierr)

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
  logical :: lredir
!
  integer :: iig, iiig, iznext
  logical :: lhelp
  real*8 :: r1, r2, thelp, mu0
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: dirdotu, gm, xi
!-- lumped quantities
  real*8 :: emitlump, caplump
  real*8 :: specig
  real*8 :: opacleak(6)
  real*8 :: probleak(6) !leakage probabilities
  real*8 :: pa !absorption probability
  real*8 :: mfphelp, pp
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer*2,pointer :: glumps(:)
  logical*2,pointer :: llumps(:)
  real*8,pointer :: capgreyinv
  real*8,pointer :: speclump
  real*8 :: dist, help

  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy,xm,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2-grd_xarr(l)**2
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  xm(l) = .5d0*(grd_xarr(l+1) + grd_xarr(l))
  dz(l) = grd_zarr(l+1) - grd_zarr(l)

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

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
  dist = min(dx(ix),dy(iy),xm(ix)*dz(iz)) * thelp

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
     else
!-- Planck x-section lump
        caplump = grd_opaclump(8,ic)*speclump
        emitlump = grd_opaclump(8,ic)*capgreyinv
        emitlump = min(emitlump,1d0)
     endif
!
!-- save
     cache%nlump = glump
     cache%emitlump = emitlump
     cache%caplump = caplump
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
  else
!-- outside the lump
     emitlump = specint0(grd_tempinv(ic),ig)*capgreyinv*grd_cap(ig,ic)
     caplump = grd_cap(ig,ic)
  endif

!-- calculate lumped values
  if(glump>0) then
!-- leakage opacities
     opacleak = grd_opaclump(1:6,ic)
  else
!!{{{
!-- calculating unlumped values
!-- inward
     if(ix/=1) l = grd_icell(ix-1,iy,iz)
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,l)+ &
        grd_sig(l))*min(dx(ix-1),dy(iy))* &
        thelp<trn_tauddmc) then
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
           thelp<trn_tauddmc
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
             thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(3)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy-1))*thelp
        opacleak(3)=(2d0/3d0)/(mfphelp*dy(iy)*thelp)
     endif

!-- upward
     if(iy==grd_ny) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy+1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),dy(iy+1)) * &
             thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(4)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy+1))*thelp
        opacleak(4)=(2d0/3d0)/(mfphelp*dy(iy)*thelp)
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
              grd_sig(l))*min(dx(ix),dy(iy), &
              xm(ix)*dz(grd_nz))*thelp<trn_tauddmc
           iznext = grd_nz
        else
           l = grd_icell(ix,iy,iz-1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),dy(iy), &
              xm(ix)*dz(iz-1))*thelp<trn_tauddmc
           iznext = iz-1
        endif

        if(lhelp) then
!-- DDMC interface
           mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
           opacleak(5)=.5d0*pp/(xm(ix)*dz(iz)*thelp)
        else
!-- DDMC interior
           mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(5)=2d0/(3d0*xm(ix)**2*dz(iz)*mfphelp*thelp**2)
        endif

!-- iz->iz+1 (opacleak(6))
        if(iz==grd_nz) then
           l = grd_icell(ix,iy,1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),dy(iy), &
              xm(ix)*dz(1))*thelp<trn_tauddmc
           iznext = 1
        else
           l = grd_icell(ix,iy,iz+1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),dy(iy), &
              xm(ix)*dz(iz+1))*thelp<trn_tauddmc
           iznext = iz+1
        endif

        if(lhelp) then
!-- DDMC interface
           mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
           opacleak(6)=.5d0*pp/(xm(ix)*dz(iz)*thelp)
        else
!-- DDMC interior
           mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(6)=2d0/(3d0*xm(ix)**2*dz(iz)*mfphelp*thelp**2)
        endif
     endif!}}}
  endif
!
!--------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
       (1d0-emitlump)*(1d0-grd_fcoef(ic))*caplump
  if(trn_isddmcanlog) then
     denom = denom+grd_fcoef(ic)*caplump
  endif
  denom = 1d0/denom

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)*denom/pc_c)
  tcensus = tsp_t1-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(trn_isddmcanlog) then
     eraddens = e*ddmct*tsp_dtinv
  else
     edep = e*(1d0-exp(-grd_fcoef(ic) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*min(dx(ix),dy(iy) , &
          xm(ix)*dz(iz))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*tsp_dtinv
     else
        eraddens = e*ddmct*tsp_dtinv
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
     ptcl2%stat = 'cens'
     return
  endif


!-- otherwise, perform event
  call rnd_r(r1,rndstate)

!-- leakage probabilities
  probleak = opacleak*denom

!-- absorption probability
  if(trn_isddmcanlog) then
     pa = grd_fcoef(ic)*caplump*denom
  else
     pa = 0d0
  endif

!-- update specarr cache only when necessary. this is slow
  if(r1>=pa .and. r1<pa+sum(probleak(1:6)) .and. glump>0 .and. &
        iand(cache%istat,2)==0) then
     cache%istat = cache%istat + 2
     call specintv(grd_tempinv(ic),grp_ng,cache%specarr,&
        mask=llumps,maskval=.true.)
  endif

!-- absorption
  if(r1<pa) then
     ptcl2%stat = 'dead'
     edep = e
     ptcl2%idist = -1

!-- ix->ix-1 leakage
  elseif(r1>=pa.and.r1<pa+probleak(1)) then
     ptcl2%idist = -3
!{{{
!-- sanity check
     if(ix==1) then
!       stop 'diffusion1: non-physical inward leakage'
        ierr = 101
        return
     endif

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
           if((grd_cap(iiig,l) + &
                grd_sig(l))*min(dx(ix-1),dy(iy),xm(ix-1)*dz(iz)) * &
                thelp<trn_tauddmc) then
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
                min(dx(ix-1),dy(iy),xm(ix-1)*dz(iz))*thelp<trn_tauddmc

     if(lhelp) then
!-- direction resample flag
        lredir = .true.
!-- sampling x,y,z
        x = grd_xarr(ix)
        call rnd_r(r1,rndstate)
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
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
     ptcl2%idist = -4
!{{{
     if(ix/=grd_nx) l = grd_icell(ix+1,iy,iz)
!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(2)
        do iig = 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(ix==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix+1),dy(iy),xm(ix+1)*dz(iz)) * &
                   thelp<trn_tauddmc
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
             min(dx(ix+1),dy(iy),xm(ix+1)*dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- direction resample flag
        lredir = .true.
!-- sampling x,y
        x = grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = grd_yarr(iy)*(1d0-r1)+grd_yarr(iy+1)*r1
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
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
           ptcl2%stat = 'flux'
!-- observer time correction
           ptcl%t=ptcl%t-(mu*y+sqrt(1d0-mu**2)*cos(om)*x)*thelp*cinv
!-- redefine for flux tally
           om = pc_pi-(z+om)  !-- direction angle
           if(om<0d0) om = om+pc_pi2
           if(om<0d0) om = om+pc_pi2
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
     ptcl2%idist = -5
!{{{
     if(iy/=1) l = grd_icell(ix,iy-1,iz)
!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(3)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy-1),xm(ix)*dz(iz)) * &
                   thelp<trn_tauddmc
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
             min(dx(ix),dy(iy-1),xm(ix)*dz(iz)) * &
             thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- direction resample flag
        lredir = .true.
!-- sampling x,y
        call rnd_r(r1,rndstate)
        x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
        y = grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
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
           ptcl2%stat = 'flux'
!-- observer time correction
           ptcl%t = ptcl%t-(mu*y+sqrt(1d0-mu**2)*cos(om)*x)*thelp*cinv
!-- redefine for flux tally
           om = pc_pi-(z+om)  !-- direction angle
           if(om<0d0) om = om+pc_pi2
           if(om<0d0) om = om+pc_pi2
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
     ptcl2%idist = -6
!{{{
     if(iy/=grd_ny) l = grd_icell(ix,iy+1,iz)
!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(4)
        do iig= 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==grd_ny) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy+1),xm(ix)*dz(iz)) * &
                   thelp<trn_tauddmc
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
             min(dx(ix),dy(iy+1),xm(ix)*dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- direction resample flag
        lredir = .true.
!-- sampling x,y
        call rnd_r(r1,rndstate)
        x = sqrt(grd_xarr(ix)**2*(1d0-r1)+grd_xarr(ix+1)**2*r1)
        y = grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
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
           ptcl2%stat = 'flux'
!-- observer time correction
           ptcl%t = ptcl%t-(mu*y+sqrt(1d0-mu**2)*cos(om)*x)*thelp*cinv
!-- redefine for flux tally
           om = pc_pi-(z+om)  !-- direction angle
           if(om<0d0) om = om+pc_pi2
           if(om<0d0) om = om+pc_pi2
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
     
!-- iz->iz-1 leakage
  elseif(r1>=pa+sum(probleak(1:4)).and.r1<pa+sum(probleak(1:5))) then
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
                min(dx(ix),dy(iy),xm(ix)*dz(iznext)) * &
                thelp<trn_tauddmc
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=.5d0*pp/(xm(ix)*dz(iz)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dz(iz) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dz(iznext))
              resopacleak=2d0/(3d0*xm(ix)**2*dz(iz)*mfphelp*thelp**2)
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),xm(ix)*dz(iznext))*thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = sqrt((1d0-r1)*grd_xarr(ix)**2+r1*grd_xarr(ix+1)**2)
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
!-- resampling y-cosing
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
        om = atan2(xi,sqrt(1d0-xi**2)*cos(pc_pi2*r1))
        if(om<0d0) om=om+pc_pi2
!-- transforming to lab
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
!-- converting to IMC
        ptcl2%itype = 1
     endif
!
!-- update particle: iz->iz-1
     iz = iznext
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- iz->iz+1 leakage
  elseif(r1>=pa+sum(probleak(1:5)).and.r1<pa+sum(probleak(1:6))) then
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
                min(dx(ix),dy(iy),xm(ix)*dz(iznext)) * &
                thelp<trn_tauddmc
!
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=.5d0*pp/(xm(ix)*dz(iz)*thelp)
           else
!-- DDMC interior
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dz(iz) + &
                   (grd_sig(l)+grd_cap(iiig,l))*dz(iznext))
              resopacleak=2d0/(3d0*xm(ix)**2*dz(iz)*mfphelp*thelp**2)
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),xm(ix)*dz(iznext)) * &
          thelp<trn_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = sqrt((1d0-r1)*grd_xarr(ix)**2+r1*grd_xarr(ix+1)**2)
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
!-- resampling y-cosing
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
        om = atan2(xi,sqrt(1d0-xi**2)*cos(pc_pi2*r1))
        if(om<0d0) om=om+pc_pi2
!-- transforming to lab
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
!-- converting to IMC
        ptcl2%itype = 1
     endif
!
!-- update particle: iz->iz+1
     iz = iznext
     ic = grd_icell(ix,iy,iz)
     ig = iiig
!}}}

!-- effective scattering
  else
     ptcl2%idist = -2
!{{{
     if(glump==grp_ng) then
!       stop 'diffusion2: effective scattering with glump==ng'
        ierr = 102
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
           iiig = nint(grd_opaclump(9,ic))
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
        denom3 = 0d0
        denom2 = 1d0/(1d0-emitlump)
        do iig = grp_ng,glump+1,-1
           iiig = glumps(iig)
           help = cache%specarr(iiig)*grd_cap(iiig,ic)*capgreyinv
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
!
     ig = iiig

     if((grd_sig(ic)+grd_cap(ig,ic))*dist < trn_tauddmc) then
        ptcl2%itype = 1
!-- direction resample flag
        lredir = .true.
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
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
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
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
!}}}
  endif

!-- update planar projections
  if(lredir) then
!-- planar projections (invariant until collision)
     ptcl2%mux = x*sin(om)/sin(z+om)  !-- intercept
     ptcl2%muy = x*sin(z)/sin(z+om)  !-- distance to intercept
     ptcl2%muz = pc_pi-(z+om)  !-- direction angle
     if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
     if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
  endif

end subroutine diffusion2
! vim: fdm=marker
