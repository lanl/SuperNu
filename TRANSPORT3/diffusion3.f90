pure subroutine diffusion3(ptcl,ptcl2,cache,rndstate,edep,eraddens,totevelo,ierr)

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
  integer :: iig, iiig
  logical :: lhelp
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: elabfact, xi, eta
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
  real*8 :: dist, help, alb, eps, beta
!
  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
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

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
  dist = min(dx(ix),dy(iy),dz(iz))*thelp

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
!{{{
!-- calculating unlumped values
!-- x left (opacleak(1))
     if(ix==1) then
        lhelp = .true.
     else
        l = grd_icell(ix-1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix-1),dy(iy),dz(iz))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)=0.5d0*pp/(thelp*dx(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))*thelp
        opacleak(1)=(2d0/3d0)/(mfphelp*dx(ix)*thelp)
     endif

!-- x right (opacleak(2))
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix+1),dy(iy),dz(iz))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=0.5d0*pp/(thelp*dx(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))*thelp
        opacleak(2)=(2d0/3d0)/(mfphelp*dx(ix)*thelp)
     endif

!-- y down (opacleak(3))
     if(iy==1) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy-1,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix),dy(iy-1),dz(iz))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(3)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy-1))*thelp
        opacleak(3)=(2d0/3d0)/(mfphelp*dy(iy)*thelp)
     endif

!-- y up (opacleak(4))
     if(iy==grd_ny) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy+1,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix),dy(iy+1),dz(iz))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dy(iy)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(4)=0.5d0*pp/(thelp*dy(iy))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dy(iy)+&
             (grd_sig(l)+grd_cap(ig,l))*dy(iy+1))*thelp
        opacleak(4)=(2d0/3d0)/(mfphelp*dy(iy)*thelp)
     endif

!-- z bottom (opacleak(5))
     if(iz==1) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy,iz-1)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix),dy(iy),dz(iz-1))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dz(iz)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(5)=0.5d0*pp/(thelp*dz(iz))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz)+&
             (grd_sig(l)+grd_cap(ig,l))*dz(iz-1))*thelp
        opacleak(5)=(2d0/3d0)/(mfphelp*dz(iz)*thelp)
     endif

!-- z top (opacleak(6))
     if(iz==grd_nz) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy,iz+1)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix),dy(iy),dz(iz+1))* &
           thelp<trn_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dz(iz)*thelp
        alb = grd_fcoef(ic)*grd_cap(ig,ic)/ &
             (grd_cap(ig,ic)+grd_sig(ic))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!        pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(6)=0.5d0*pp/(thelp*dz(iz))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz)+&
             (grd_sig(l)+grd_cap(ig,l))*dz(iz+1))*thelp
        opacleak(6)=(2d0/3d0)/(mfphelp*dz(iz)*thelp)
     endif
!}}}
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
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(trn_isddmcanlog) then
     eraddens = e*ddmct*tsp_dtinv
  else
     edep = e * &
          (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = &
             e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*tsp_dtinv
     else
        eraddens = eraddens + e*ddmct*tsp_dtinv
     endif
!
     if(edep/=edep) then
!       write(0,*) e,grd_fcoef(ic),caplump,ddmct,glump,speclump,ig,grd_tempinv(ic)
!       stop 'diffusion3: invalid energy deposition'
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
  if (ddmct /= tau) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
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
  if(r1>=pa .and. r1<pa+sum(probleak) .and. glump>0 .and. &
        iand(cache%istat,2)==0) then
     cache%istat = cache%istat + 2
     call specintv(grd_tempinv(ic),grp_ng,cache%specarr,&
        mask=llumps,maskval=.true.)
  endif

!-- absorption
  if(r1<pa) then
     ptcl2%idist = -1
     ptcl2%isvacant = .true.
     ptcl2%done = .true.
     edep = e

!-- ix->ix-1 leakage
  elseif(r1>=pa.and.r1<pa+probleak(1)) then
     ptcl2%idist = -3
!{{{
     if(ix/=1) l = grd_icell(ix-1,iy,iz)

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
           if(ix==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix-1),dy(iy),dz(iz))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dx(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dx(ix)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dx(ix-1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dx(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(ix==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix-1),dy(iy),dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        x = grd_xarr(ix)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        xi = -max(r1,r2)
        call rnd_r(r1,rndstate)
        eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
        om = atan2(eta,xi)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
        if(ix==1) then
!-- escaping at ix=1
           ptcl2%isvacant = .true.
           ptcl2%done = .true.
           ptcl2%lflux = .true.
           return
        else
!-- converting to IMC
           ptcl2%itype = 1
!-- ix->ix-cell(ix,iy,iz)
        endif
     endif
!
!-- update particle: ix->ix-1
     ix = ix-1
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- ix->ix+1 leakage
  elseif(r1>=pa+probleak(1).and.r1<pa+sum(probleak(1:2))) then
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
                   min(dx(ix+1),dy(iy),dz(iz))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dx(ix)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dx(ix))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dx(ix)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dx(ix+1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dx(ix))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
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
             min(dx(ix+1),dy(iy),dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        x = grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        xi = max(r1,r2)
        call rnd_r(r1,rndstate)
        eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
        om = atan2(eta,xi)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
     ig = iiig!}}}

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
        do iig = 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy-1),dz(iz))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dy(iy)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
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

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iy==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix),dy(iy-1),dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = (1d0-r1)*grd_xarr(ix)+r1*grd_xarr(ix+1)
        y = grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        eta = -max(r1,r2)
        call rnd_r(r1,rndstate)
        xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
        om = atan2(eta,xi)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
        if(iy==1) then
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
     ig = iiig!}}}

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
        do iig = 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iy==grd_ny) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy+1),dz(iz))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dy(iy)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
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
             min(dx(ix),dy(iy+1),dz(iz))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = (1d0-r1)*grd_xarr(ix)+r1*grd_xarr(ix+1)
        y = grd_yarr(iy+1)
        call rnd_r(r1,rndstate)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        eta = max(r1,r2)
        call rnd_r(r1,rndstate)
        xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
        mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
        om = atan2(eta,xi)
        if(om<0d0) om=om+pc_pi2
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
        if(iy==grd_ny) then
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
     ig = iiig!}}}

!-- iz->iz-1 leakage
  elseif(r1>=pa+sum(probleak(1:4)).and.r1<pa+sum(probleak(1:5))) then
     ptcl2%idist = -7
!{{{
     if(iz/=1) l = grd_icell(ix,iy,iz-1)

!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(5)
        do iig = 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iz==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy),dz(iz-1))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dz(iz)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dz(iz))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dz(iz)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dz(iz-1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dz(iz))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iz==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix),dy(iy),dz(iz-1))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = (1d0-r1)*grd_xarr(ix)+r1*grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        z = grd_zarr(iz)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = -max(r1,r2)
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
        xi = sqrt(1d0-mu**2)*cos(om)
        eta = sqrt(1d0-mu**2)*sin(om)
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
        if(iz==1) then
!-- escaping at iz=1
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
!-- update particle: iz->iz-1
     iz = iz-1
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- iz->iz+1 leakage
     ptcl2%idist = -8
  elseif(r1>=pa+sum(probleak(1:5)).and.r1<pa+sum(probleak(1:6))) then
!{{{
     if(iz/=grd_nz) l = grd_icell(ix,iy,iz+1)

!-- sampling next group
     if(glump==0) then
        iiig = ig
     else
        call rnd_r(r1,rndstate)
        denom2 = 0d0
        help = 1d0/opacleak(6)
        do iig = 1, glump
           iiig = glumps(iig)
           specig = cache%specarr(iiig)
!-- calculating resolved leakage opacities
           if(iz==grd_nz) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix),dy(iy),dz(iz+1))*thelp<trn_tauddmc
           endif
           if(lhelp) then
!-- IMC interface or boundary
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   dz(iz)*thelp
              alb = grd_fcoef(ic)*grd_cap(iiig,ic)/ &
                   (grd_cap(iiig,ic)+grd_sig(ic))
              eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
              beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
                   2.25d0*alb**2*mfphelp**4)
              pp = 0.5d0*eps*beta/(beta-0.75*eps*mfphelp)
!              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 0.5d0*pp/(thelp*dz(iz))
           else
!-- DDMC interface
              mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic)) * &
                   dz(iz)+&
                   (grd_sig(l)+grd_cap(iiig,l)) * &
                   dz(iz+1))*thelp
              resopacleak = (2d0/3d0)/(mfphelp*thelp*dz(iz))
           endif
           denom2 = denom2 + specig*resopacleak*speclump*help
           if(r1<denom2) exit
        enddo
     endif

!-- sampling wavelength
     call rnd_r(r1,rndstate)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(iz==grd_nz) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix),dy(iy),dz(iz+1))*thelp<trn_tauddmc
     endif

     if(lhelp) then
!-- sampling x,y,z
        call rnd_r(r1,rndstate)
        x = (1d0-r1)*grd_xarr(ix)+r1*grd_xarr(ix+1)
        call rnd_r(r1,rndstate)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        z = grd_zarr(iz+1)
!-- sampling direction
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
        xi = sqrt(1d0-mu**2)*cos(om)
        eta = sqrt(1d0-mu**2)*sin(om)
        if(grd_isvelocity) then
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
        else
           elabfact = 1d0
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
!-- transforming mu to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- transforming om to lab
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
        if(iz==grd_nz) then
!-- escaping at iz=nz
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
!-- update particle: iz->iz+1
     iz = iz+1
     ic = grd_icell(ix,iy,iz)
     ig = iiig!}}}

!-- effective scattering
  else
     ptcl2%idist = -2
!!{{{
     if(glump==grp_ng) then
!       stop 'diffusion3: effective scattering with glump==ng'
        ierr = 102
        return
     endif
!
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
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
!-- direction sampled isotropically           
        call rnd_r(r1,rndstate)
        mu = 1d0 - 2d0*r1
        call rnd_r(r1,rndstate)
        om = pc_pi2*r1
        xi = sqrt(1d0-mu**2)*cos(om)
        eta = sqrt(1d0-mu**2)*sin(om)
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        x = r1*grd_xarr(ix+1)+(1d0-r1)*grd_xarr(ix)
        call rnd_r(r1,rndstate)
        y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
        call rnd_r(r1,rndstate)
        z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- doppler and aberration corrections
        if(grd_isvelocity) then
!-- calculating transformation factors
           elabfact = 1d0+(x*xi+y*eta+z*mu)*cinv
!-- transforming z-axis direction cosine to lab
           mu = (mu+z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
           om = atan2(eta+y*cinv,xi+x*cinv)
           if(om<0d0) om=om+pc_pi2
!-- ELABFACT LAB RESET
           xi = sqrt(1d0-mu**2)*cos(om)
           eta= sqrt(1d0-mu**2)*sin(om)
           elabfact=1d0-(x*xi+y*eta+z*mu)*cinv
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
     endif!}}}

  endif

end subroutine diffusion3
