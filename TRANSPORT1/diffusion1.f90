subroutine diffusion1(ptcl,ig,isvacant,icell,specarr)

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
  integer,intent(inout) :: ig
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
  real*8 :: help
!
  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l) = grd_xarr(l+1)**2 - grd_xarr(l)**2
  dx3(l) = grd_xarr(l+1)**3 - grd_xarr(l)**3
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)

  ix => ptcl%ix
  iy => ptcl%iy
  iz => ptcl%iz
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
  tempinv = 1d0/grd_temp(ix,iy,iz)
  capgreyinv = max(1d0/grd_capgrey(ix,iy,iz),0d0) !catch nans

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
  if(grd_cap(ig,ix,iy,iz)*min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
       thelp>=prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ix,iy,iz)*min(dx(ix),xm(ix)*dyac(iy) , &
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
  if(glump>0 .and. .not.all(icell==[ix,iy,iz])) then
     icell = [ix,iy,iz]
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
        caplump = grd_capgrey(ix,iy,iz)
     else
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- emission lump
           emitlump = emitlump + specig*capgreyinv*grd_cap(iiig,ix,iy,iz)
!-- Planck x-section lump
           caplump = caplump + specig*grd_cap(iiig,ix,iy,iz)*speclump
        enddo
        emitlump = min(emitlump,1d0)
     endif
!-- leakage opacities
     opacleak = grd_opacleak(:,ix,iy,iz)
!!}}}
  else
!
!-- calculating unlumped values
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ix,iy,iz)!{{{
     caplump = grd_cap(ig,ix,iy,iz)

!-- ix->ix-1 in (opacleak(1))
     if(ix==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix-1,iy,iz)+ &
           grd_sig(ix-1,iy,iz))*min(dx(ix-1),xm(ix-1)*dyac(iy) , &
           xm(ix-1)*ym(iy)*dz(iz))*thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*dx(ix)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(1)=1.5d0*pp*grd_xarr(ix)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dx(ix)+&
             (grd_sig(ix-1,iy,iz)+grd_cap(ig,ix-1,iy,iz))*dx(ix-1))*thelp
        opacleak(1)=2d0*grd_xarr(ix)**2/(dx3(ix)*thelp*help)
     endif

!
!-- ix->ix+1 (opacleak(2))
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix+1,iy,iz)+grd_sig(ix+1,iy,iz)) * &
             min(dx(ix+1),xm(ix+1)*dyac(iy),xm(ix+1)*ym(iy)*dz(iz)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*dx(ix)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(2)=1.5d0*pp*grd_xarr(ix+1)**2/(dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dx(ix)+&
             (grd_sig(ix+1,iy,iz)+grd_cap(ig,ix+1,iy,iz))*dx(ix+1))*thelp
        opacleak(2)=2d0*grd_xarr(ix+1)**2/(dx3(ix)*thelp*help)
     endif

!
!-- iy->iy-1 (opacleak(3))
     if(iy==1) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix,iy-1,iz)+ &
             grd_sig(ix,iy-1,iz))*min(dx(ix),xm(ix)*dyac(iy-1), &
             xm(ix)*ym(iy-1)*dz(iz))*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(3)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dyac(iy) + &
             (grd_sig(ix,iy-1,iz)+grd_cap(ig,ix,iy-1,iz))*dyac(iy-1))
        opacleak(3)=2d0*sqrt(1d0-grd_yarr(iy)**2)*dx(ix) / &
             (dy(iy)*dx3(ix)*help*thelp**2)
     endif

!
!-- iy->iy+1 (opacleak(4))
     if(iy==grd_ny) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix,iy+1,iz)+ &
             grd_sig(ix,iy+1,iz))*min(dx(ix),xm(ix)*dyac(iy+1), &
             xm(ix)*ym(iy+1)*dz(iz))*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*xm(ix)*dyac(iy)*thelp
        pp = 4d0/(3d0*help+6d0*pc_dext)
        opacleak(4)=0.75d0*pp*dx2(ix)*sqrt(1d0-grd_yarr(iy+1)**2) / &
             (dy(iy)*dx3(ix)*thelp)
     else
!-- DDMC interior
        help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dyac(iy) + &
             (grd_sig(ix,iy+1,iz)+grd_cap(ig,ix,iy+1,iz))*dyac(iy+1))
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
           lhelp = (grd_cap(ig,ix,iy,grd_nz)+ &
              grd_sig(ix,iy,grd_nz))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(grd_nz))*thelp<prt_tauddmc
           iznext = grd_nz
        else
           lhelp = (grd_cap(ig,ix,iy,iz-1)+ &
              grd_sig(ix,iy,iz-1))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz-1))*thelp<prt_tauddmc
           iznext = iz-1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           opacleak(5)=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
        else
!-- DDMC interior
           help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dz(iz) + &
                (grd_sig(ix,iy,iznext)+grd_cap(ig,ix,iy,iznext))*dz(iznext))
           grd_opacleak(5)=2d0*dyac(iy)*dx(ix) / &
                (ym(iy)*dy(iy)*dz(iz)*dx3(ix)*help*thelp**2)
        endif
!-- iz->iz+1 (opacleak(6))
        if(iz==grd_nz) then
           lhelp = (grd_cap(ig,ix,iy,1)+ &
              grd_sig(ix,iy,1))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(1))*thelp<prt_tauddmc
           iznext = 1
        else
           lhelp = (grd_cap(ig,ix,iy,iz+1)+ &
              grd_sig(ix,iy,iz+1))*min(dx(ix),xm(ix)*dyac(iy), &
              xm(ix)*ym(iy)*dz(iz+1))*thelp<prt_tauddmc
           iznext = iz+1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz))*xm(ix)*ym(iy) * &
              dz(iz)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           opacleak(6)=0.75d0*pp*dx2(ix)*dyac(iy)/(dy(iy)*dx3(ix)*dz(iz))
        else
!-- DDMC interior
           help = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz))*dz(iz) + &
                (grd_sig(ix,iy,iznext)+grd_cap(ig,ix,iy,iznext))*dz(iznext))
           grd_opacleak(6)=2d0*dyac(iy)*dx(ix) / &
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
       (1d0-emitlump)*(1d0-grd_fcoef(ix,iy,iz))*caplump
  if(prt_isddmcanlog) then
     denom = denom+grd_fcoef(ix,iy,iz)*caplump
  endif

  r1 = rnd_r(rnd_state)
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     grd_eraddens(ix,iy,iz)= grd_eraddens(ix,iy,iz)+e*ddmct*dtinv
  else
     grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz)+e * &
          (1d0-exp(-grd_fcoef(ix,iy,iz)*caplump*pc_c*ddmct))
     if(grd_fcoef(ix,iy,iz)*caplump*min(dx(ix),xm(ix)*dyac(iy) , &
          xm(ix)*ym(iy)*dz(iz))*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ix,iy,iz)*caplump)
        grd_eraddens(ix,iy,iz)= &
             grd_eraddens(ix,iy,iz)+e* &
             (1d0-exp(-grd_fcoef(ix,iy,iz)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        grd_eraddens(ix,iy,iz) = grd_eraddens(ix,iy,iz)+e*ddmct*dtinv
     endif
!
     if(grd_edep(ix,iy,iz)/=grd_edep(ix,iy,iz)) then
!       write(0,*) e,grd_fcoef(ix,iy,iz),caplump,ddmct,glump,speclump,ig,tempinv
        stop 'diffusion3: invalid energy deposition'
     endif
     e = e*exp(-grd_fcoef(ix,iy,iz)*caplump*pc_c*ddmct)
!!}}}
  endif

end subroutine diffusion1
