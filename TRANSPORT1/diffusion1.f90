subroutine diffusion1(ptcl,ptcl2,icspec,specarr)

  use randommod
  use miscmod
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
  type(packet2),target,intent(inout) :: ptcl2
  integer,intent(inout) :: icspec
  real*8,intent(inout) :: specarr(grp_ng)
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer,external :: emitgroup
!
  integer :: iig, iiig, imu, iom, iznext
  logical :: lhelp
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: elabfact, xi, eta
!-- lumped quantities
  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleak(6)
  real*8 :: probleak(6) !leakage probabilities
  real*8 :: pa !absorption probability
  real*8 :: mfphelp, pp
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer :: glumps(grp_ng)
  real*8 :: dtinv, tempinv, capgreyinv
  real*8 :: help, mux, muy, muz
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
!
!-- shortcuts
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
  if(grd_cap(ig,ic)*min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
       thelp>=prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ic)*min(dx(ix),xm(ix)*dyac(iy) , &
             xm(ix)*ym(iy)*dz(iz))*thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
  endif
! write(0,*) ipart,istep,glump,g,ix,iy,iz

!
!-- only do this if needed
  if(glump>0 .and. icspec/=ic) then
     icspec = ic
     specarr = specintv(tempinv,0) !this is slow!
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
     opacleak = grd_opacleak(:,ic)
!!}}}
  else
!
!-- calculating unlumped values
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ic)!{{{
     caplump = grd_cap(ig,ic)

!-- ix->ix-1 in (opacleak(1))
     if(ix==1) then
        lhelp = .true.
     else
        l = grd_icell(ix-1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*min(dx(ix-1),xm(ix-1)*dyac(iy) , &
           xm(ix-1)*ym(iy)*dz(iz))*thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(1)=1.5d0*pp*grd_xarr(ix)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))*thelp
        opacleak(1)=2d0*grd_xarr(ix)**2/(dx3(ix)*thelp*help)
     endif

!
!-- ix->ix+1 (opacleak(2))
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+grd_sig(l)) * &
             min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy)*dz(iz)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(2)=1.5d0*pp*grd_xarr(ix+1)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))*thelp
        opacleak(2)=2d0*grd_xarr(ix+1)**2/(dx3(ix)*thelp*help)
     endif

!
!-- iy->iy-1 (opacleak(3))
     if(iy==1) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy-1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy-1), &
             xm(ix)*ym(iy-1)*dz(iz))*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(3)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ic)+grd_cap(ig,ic))*dyac(iy) + &
             (grd_sig(l)+grd_cap(ig,l))*dyac(iy-1))
        opacleak(3)=2d0*sqrt(1d0-grd_yarr(iy)**2)*dx(ix) / &
             (dy(iy)*dx3(ix)*help*thelp**2)
     endif

!
!-- iy->iy+1 (opacleak(4))
     if(iy==grd_ny) then
        lhelp = .true.
     else
        l = grd_icell(ix,iy+1,iz)
        lhelp = (grd_cap(ig,l)+ &
             grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy+1), &
             xm(ix)*ym(iy+1)*dz(iz))*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(4)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy+1)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ic)+grd_cap(ig,ic))*dyac(iy) + &
             (grd_sig(l)+grd_cap(ig,l))*dyac(iy+1))
        opacleak(4)=2d0*sqrt(1d0-grd_yarr(iy+1)**2)*dx(ix) / &
             (dy(iy)*dx3(ix)*help*thelp**2)
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
              xm(ix)*ym(iy)*dz(grd_nz))*thelp<prt_tauddmc
           iznext = grd_nz
        else
           l = grd_icell(ix,iy,iz-1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz-1))*thelp<prt_tauddmc
           iznext = iz-1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           opacleak(5)=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
        else
!-- DDMC interior
           help = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(5)=2d0*dyac(iy)*dx(ix) / &
                (ym(iy)*dy(iy)*dz(iz)*dx3(ix)*help*thelp**2)
        endif
!-- iz->iz+1 (opacleak(6))
        if(iz==grd_nz) then
           l = grd_icell(ix,iy,1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(1))*thelp<prt_tauddmc
           iznext = 1
        else
           l = grd_icell(ix,iy,iz+1)
           lhelp = (grd_cap(ig,l)+ &
              grd_sig(l))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz+1))*thelp<prt_tauddmc
           iznext = iz+1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,ic)+grd_sig(ic))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           opacleak(6)=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
        else
!-- DDMC interior
           help = ((grd_sig(ic)+grd_cap(ig,ic))*dz(iz) + &
                (grd_sig(l)+grd_cap(ig,l))*dz(iznext))
           opacleak(6)=2d0*dyac(iy)*dx(ix) / &
                (ym(iy)*dy(iy)*dz(iz)*dx3(ix)*help*thelp**2)
        endif
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

  r1 = rnd_r(rnd_state)
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     grd_eraddens(ic)= grd_eraddens(ic)+e*ddmct*dtinv
  else
     grd_edep(ic) = grd_edep(ic)+e * &
          (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*min(dx(ix),xm(ix)*dyac(iy) , &
          xm(ix)*ym(iy)*dz(iz))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        grd_eraddens(ic)= &
             grd_eraddens(ic)+e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        grd_eraddens(ic) = grd_eraddens(ic)+e*ddmct*dtinv
     endif
!
     if(grd_edep(ic)/=grd_edep(ic)) then
!       write(0,*) e,grd_fcoef(ic),caplump,ddmct,glump,speclump,ig,tempinv
        stop 'diffusion1: invalid energy deposition'
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
     grd_numcensus(ic) = grd_numcensus(ic)+1
     return
  endif

!-- otherwise, perform event
  r1 = rnd_r(rnd_state)
  help = 1d0/denom

!-- leakage probabilities
  probleak = opacleak*help
!  write(*,*) probleak
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
     grd_edep(ic) = grd_edep(ic)+e

!-- ix->ix-1 leakage
  elseif(r1>=pa.and.r1<pa+probleak(1)) then
!-- sanity check
     if(ix==1) stop 'diffusion1: invalid probleak(1)'
!
     l = grd_icell(ix-1,iy,iz)
!-- sample next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(1)
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix-1),xm(ix-1)*dyac(iy),xm(ix-1)*ym(iy) * &
                   dz(iz))*thelp<prt_tauddmc
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix-1),xm(ix-1)*dyac(iy),xm(ix-1)*ym(iy)*dz(iz)) * &
          thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        x = grd_xarr(ix)
        r1 = rnd_r(rnd_state)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        r1 = rnd_r(rnd_state)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
        z = min(z,grd_zarr(iz+1))
        z = max(z,grd_zarr(iz))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        mu = -max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
     endif
!-- ix->ix-1
     ix = ix-1
     ic = grd_icell(ix,iy,iz)

!-- ix->ix+1 leakage
  elseif(r1>=pa+probleak(1).and.r1<pa+sum(probleak(1:2))) then

!-- sampling next group
     if(ix/=grd_nx) l = grd_icell(ix+1,iy,iz)
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(2)
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           if(ix==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                   min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy) * &
                   dz(iz))*thelp<prt_tauddmc
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
             min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy)*dz(iz)) * &
             thelp<prt_tauddmc
     endif

     if(.not.lhelp) then
!-- ix->ix+1
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
     else
!-- sampling x,y,z
        x = grd_xarr(ix+1)
        r1 = rnd_r(rnd_state)
        y = (1d0-r1)*grd_yarr(iy)+r1*grd_yarr(iy+1)
        r1 = rnd_r(rnd_state)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
        z = min(z,grd_zarr(iz+1))
        z = max(z,grd_zarr(iz))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        mu = max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
        if(ix==grd_nx) then
!-- escaping at ix=nx
           ptcl2%isvacant = .true.
           ptcl2%done = .true.
!-- luminosity tally
           eta = sqrt(1d0-mu**2)*cos(om)
           xi = sqrt(1d0-mu**2)*sin(om)
           mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
           muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
           muz = mu*y-eta*sqrt(1d0-y**2)
           help = atan2(muy,mux)
           if(help<0d0) help=help+pc_pi2
!-- retrieving lab frame flux group, polar, azimuthal bin
           iom = binsrch(help,flx_om,flx_nom+1)
           imu = binsrch(muz,flx_mu,flx_nmu+1)
           ig = binsrch(wl,flx_wl,flx_ng+1)
!-- checking group bounds
           if(ig>flx_ng.or.ig<1) then
              if(ig>flx_ng) then
                 ig=flx_ng
              else
                 ig=1
              endif
           endif
!-- tallying outbound energy
           tot_eout = tot_eout+e
!-- tallying outbound luminosity
           flx_luminos(ig,imu,iom) = flx_luminos(ig,imu,iom)+e*dtinv
           flx_lumdev(ig,imu,iom) = flx_lumdev(ig,imu,iom)+(e*dtinv)**2
           flx_lumnum(ig,imu,iom) = flx_lumnum(ig,imu,iom)+1
           return
        else
!-- converting to IMC
           ptcl2%itype = 1
           grd_methodswap(ic) = grd_methodswap(ic)+1
!-- ix->ix+1
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        endif
     endif

!-- iy->iy-1 leakage
  elseif(r1>=pa+sum(probleak(1:2)).and.r1<pa+sum(probleak(1:3))) then
!-- sanity check
     if(grd_ny==1) stop 'diffusion1: probleak(3) and grd_ny=1'
     if(iy==1) stop 'diffusion1: invalid probleak(3)'

!-- sampling next group
     l = grd_icell(ix,iy-1,iz)
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(3)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy-1),xm(ix)*ym(iy-1)*dz(iz)) * &
                thelp<prt_tauddmc
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

!-- sampling wavelength
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy-1),xm(ix)*ym(iy-1)*dz(iz)) * &
          thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        r1 = rnd_r(rnd_state)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        y = grd_yarr(iy)
        r1 = rnd_r(rnd_state)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
        z = min(z,grd_zarr(iz+1))
        z = max(z,grd_zarr(iz))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        eta = max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
     endif
!-- iy->iy+1
     iy = iy-1
     ic = grd_icell(ix,iy,iz)

!-- iy->iy+1 leakage
  elseif(r1>=pa+sum(probleak(1:3)).and.r1<pa+sum(probleak(1:4))) then
!-- sanity check
     if(grd_ny==1) stop 'diffusion1: probleak(4) and grd_ny=1'
     if(iy==grd_ny) stop 'diffusion1: invalid probleak(4)'

!-- sampling next group
     l = grd_icell(ix,iy+1,iz)
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(4)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy+1),xm(ix)*ym(iy+1)*dz(iz)) * &
                thelp<prt_tauddmc
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- checking adjacent
     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy+1),xm(ix)*ym(iy+1)*dz(iz)) * &
          thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        r1 = rnd_r(rnd_state)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        y = grd_yarr(iy+1)
        r1 = rnd_r(rnd_state)
        z = (1d0-r1)*grd_zarr(iz)+r1*grd_zarr(iz+1)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
        z = min(z,grd_zarr(iz+1))
        z = max(z,grd_zarr(iz))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        eta = -max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
     endif
!-- iy->iy+1
     iy = iy+1
     ic = grd_icell(ix,iy,iz)

!-- iz->iz-1 leakage
  elseif(r1>=pa+sum(probleak(1:4)).and.r1<pa+sum(probleak(1:5))) then
!-- sanity check
     if(grd_nz==1) stop 'diffusion1: probleak(5) and grd_nz=1'
!-- setting helper index
     if(iz==1) then
        iznext=grd_nz
     else
        iznext=iz-1
     endif
     l = grd_icell(ix,iy,iznext)
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(5)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
                thelp<prt_tauddmc
!
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*ym(iy)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),dz(iznext))*thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        r1 = rnd_r(rnd_state)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        r1 = rnd_r(rnd_state)
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
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        xi = -max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
     endif
!-- iz->iz-1
     iz = iznext
     ic = grd_icell(ix,iy,iz)

!-- iz->iz+1 leakage
  elseif(r1>=pa+sum(probleak(1:5)).and.r1<pa+sum(probleak(1:6))) then
!-- sanity check
     if(grd_nz==1) stop 'diffusion1: probleak(6) and grd_nz=1'
!-- setting helper index
     if(iz==grd_nz) then
        iznext=1
     else
        iznext=iz+1
     endif
     l = grd_icell(ix,iy,iznext)
!-- sampling next group
     if(speclump<=0d0) then
        iiig = ig
     else
        r1 = rnd_r(rnd_state)
        denom2 = 0d0
        help = 1d0/opacleak(6)
        do iig = 1, glump
           iiig=glumps(iig)
           specig = specarr(iiig)
!-- calculating resolved leakage opacities
           lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
                min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iznext)) * &
                thelp<prt_tauddmc
!
           if(lhelp) then
!-- DDMC interface
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic)) * &
                   xm(ix)*ym(iy)*dz(iz)*thelp
              pp = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

     lhelp = (grd_cap(iiig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),dz(iznext))*thelp<prt_tauddmc

     if(lhelp) then
!-- sampling x,y,z
        r1 = rnd_r(rnd_state)
        x = ((1d0-r1)*grd_xarr(ix)**3+r1*grd_xarr(ix+1)**3)**(1d0/3d0)
        r1 = rnd_r(rnd_state)
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
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
!-- sampling direction
        r1 = rnd_r(rnd_state)
        r2 = rnd_r(rnd_state)
        xi = max(r1,r2)
        r1 = rnd_r(rnd_state)
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
!-- converting to IMC
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
     endif
!-- iz->iz+1
     iz = iznext
     ic = grd_icell(ix,iy,iz)

!-- effective scattering
  else
!
     if(glump==grp_ng) stop 'diffusion1: effective scattering with glump==ng'
!
     r1 = rnd_r(rnd_state)

     if(glump==0) then
        iiig = emitgroup(r1,ic)
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
     r1 = rnd_r(rnd_state)
     wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))

     if((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) &
          *thelp < prt_tauddmc) then
        ptcl2%itype = 1
        grd_methodswap(ic) = grd_methodswap(ic)+1
!-- direction sampled isotropically           
        r1 = rnd_r(rnd_state)
        mu = 1d0 - 2d0*r1
        r1 = rnd_r(rnd_state)
        om = pc_pi2*r1
!-- position sampled uniformly
        r1 = rnd_r(rnd_state)
        x = (r1*grd_xarr(ix+1)**3+(1d0-r1)*grd_xarr(ix)**3)**(1d0/3d0)
        r1 = rnd_r(rnd_state)
        y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
        r1 = rnd_r(rnd_state)
        z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
        y = min(y,grd_yarr(iy+1))
        y = max(y,grd_yarr(iy))
        z = min(z,grd_zarr(iz+1))
        z = max(z,grd_zarr(iz))
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
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- transforming energy weights to lab
           e = e*help
           e0 = e0*help
        endif
     endif

  endif

end subroutine diffusion1
