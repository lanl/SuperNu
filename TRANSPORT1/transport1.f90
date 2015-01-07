subroutine transport1(ptcl,ic,ig,isvacant)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use totalsmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  integer,intent(inout) :: ic, ig
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer,external :: emitgroup

  logical :: lout
  integer :: imu, iom, ihelp
  real*8 :: elabfact, eta, xi, mux,muy,muz
  real*8 :: dtinv, thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop,d
  real*8 :: r1,r2

  integer :: iynext,iynext1,iynext2,iznext,iyold,idby1,idby2
  integer :: idby1old,idby2old
  real*8 :: yhelp1,yhelp2,yhelp3,yhelp4,dby1,dby2
  real*8 :: dby1old,dby2old
  real*8 :: zhelp
  real*8 :: xold,yold,muold,dbyold,etaold,muzold

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
! dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5d0*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25d0*(grd_yarr(l+1)+grd_yarr(l))**2)

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
!-- shortcut
  dtinv = 1d0/tsp_dt
!-- spherical projections
  eta = sqrt(1d0-mu**2)*cos(om)
  xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
  mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
  muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
  muz = mu*y-eta*sqrt(1d0-y**2)
  idby1=0
  idby2=0
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     elabfact=1d0-mu*x*cinv
     thelp = tsp_t
  else
     elabfact = 1d0
     thelp = 1d0
  endif
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- census distance
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- radial boundary distance (x)
  if (ix == 1) then
     dbx = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  elseif (mu < -sqrt(1d0-(grd_xarr(ix)/x)**2)) then
     dbx = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
  else
     dbx = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  endif
!-- sanity check
  if(dbx/=dbx) stop 'transport1: dbx/=dbx'
!
!-- polar boundary distance (y)
!-- dby1: iy->iy-1
  yhelp1=grd_yarr(iy)**2-muz**2
  yhelp2=mu*grd_yarr(iy)**2-muz*y
  yhelp3=grd_yarr(iy)**2-y**2
  if(yhelp1==0d0.and.yhelp3==0d0) then
     idby1=1
!-- particle, direction on cone
     dby1 = 2d0*pc_c*tsp_dt*thelpinv
     iynext1=iy
  elseif(yhelp1==0d0) then
     if((muz>=0d0.and.muz==grd_yarr(iy)).or. &
          (muz>=0d0.and.muz==-grd_yarr(iy)).or. &
          yhelp2==0d0) then
        idby1=2
!-- direction parallel to lower cone or nonphysical
        dby1 = 2d0*pc_c*tsp_dt*thelpinv
        iynext1=iy
     else
        idby1=3
        dby1=-0.5*x*yhelp3/yhelp2
        iynext1=iy-1
     endif
  elseif(abs(yhelp3)<1d-15*abs(y)) then
!-- particle on lower cone
     if(abs(y-grd_yarr(iy))<1d-15*abs(y)) then

        dby1=-2d0*x*yhelp2/yhelp1
        if(dby1>0d0) then
           if(cos(om)>=0d0) then
              write(*,*) yold, y, iy, x, ix
              write(*,*) dbyold, iyold, muzold, idby1old, idby2old
              write(*,*) dby1old, dby2old
              stop 'transport1: did not cross lower y bound'
           elseif(grd_yarr(iy)<=0d0) then
!-- choose dby2
              idby1=4
              dby1 = 2d0*pc_c*tsp_dt*thelpinv
              iynext1=iy
           elseif(grd_yarr(iy)>0d0) then
!-- cone internal transfer
              idby1=5
              iynext1=iy-1
           else
              stop 'transport1: ptcl on cone and dby1 invalid'
           endif
        else
           idby1=6
           dby1 = 2d0*pc_c*tsp_dt*thelpinv
           iynext1=iy
        endif

     else
        idby1=7
        dby1=-2d0*x*yhelp2/yhelp1
        iynext1=iy-1
     endif

  else
     yhelp4=yhelp2**2-yhelp1*yhelp3
     if(abs(yhelp4)<1d-12*abs(yhelp2)) yhelp4=0d0
     if(yhelp4<0d0) then
        idby1=8
!-- not intersecting lower cone
        dby1 = 2d0*pc_c*tsp_dt*thelpinv
        iynext1=iy
     else
!-- intersecting lower cone at at least one point
        if(cos(om)<0d0.and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
!-- choose dby2
           idby1=9
           dby1 = 2d0*pc_c*tsp_dt*thelpinv
           iynext1=iy
        else
           idby1=10
           yhelp4=sqrt(yhelp4)
           yhelp1=1d0/yhelp1
           help=x*(-yhelp2+yhelp4)*yhelp1
           dby1=x*(-yhelp2-yhelp4)*yhelp1
           if(help<0d0) help=2d0*pc_c*tsp_dt*thelpinv
           if(dby1<0d0) dby1=2d0*pc_c*tsp_dt*thelpinv
           dby1=min(help,dby1)
           iynext1=iy-1
        endif
     endif
  endif

!-- dby2: iy->iy+1
  yhelp1=grd_yarr(iy+1)**2-muz**2
  yhelp2=mu*grd_yarr(iy+1)**2-muz*y
  yhelp3=grd_yarr(iy+1)**2-y**2
  if(yhelp1==0d0.and.yhelp3==0d0) then
     idby2=1
!-- particle, direction on cone
     dby2 = 2d0*pc_c*tsp_dt*thelpinv
     iynext2=iy
  elseif(yhelp1==0d0) then
     if((muz<=0d0.and.muz==-grd_yarr(iy+1)).or. &
          (muz<=0d0.and.muz==grd_yarr(iy+1)).or. &
          yhelp2==0d0) then
        idby2=2
!-- direction parallel to upper cone or nonphysical
        dby2 = 2d0*pc_c*tsp_dt*thelpinv
        iynext2=iy
     else
        idby2=3
        dby2=-0.5*x*yhelp3/yhelp2
        iynext2=iy+1
     endif
  elseif(abs(yhelp3)<1d-15*abs(y)) then
!-- particle on upper cone
     if(abs(y-grd_yarr(iy+1))<1d-15*abs(y)) then

        dby2=-2d0*x*yhelp2/yhelp1
        if(dby2>0d0) then
           if(cos(om)<0d0) then
              write(*,*) yold, y, iy, x, ix
              write(*,*) dbyold, iyold, muzold, idby1old, idby2old
              write(*,*) dby1old, dby2old
              stop 'transport1: did not cross upper y bound'
           elseif(grd_yarr(iy+1)>=0d0) then
!-- choose dby1
              idby2=4
              dby2 = 2d0*pc_c*tsp_dt*thelpinv
              iynext2=iy
           elseif(grd_yarr(iy+1)<0d0) then
!-- cone internal transfer
              idby2=5
              iynext2=iy+1
           else
              stop 'transport1: ptcl on cone and dby2 invalid'
           endif
        else
           idby2=6
           dby2 = 2d0*pc_c*tsp_dt*thelpinv
           iynext2=iy
        endif
     else
        idby2=7
        dby2=-2d0*x*yhelp2/yhelp1
        iynext2=iy+1
     endif

  else
     yhelp4=yhelp2**2-yhelp1*yhelp3
     if(abs(yhelp4)<1d-12*abs(yhelp2)) yhelp4=0d0
     if(yhelp4<0d0) then
        idby2=8
!-- not intersecting upper cone
        dby2 = 2d0*pc_c*tsp_dt*thelpinv
        iynext2=iy
     else
!-- intersecting upper cone at at least one point
        if(cos(om)>=0d0.and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
!-- choose dby1
           idby2=9
           dby2 = 2d0*pc_c*tsp_dt*thelpinv
           iynext2=iy
        else
           idby2=10
           yhelp4=sqrt(yhelp4)
           yhelp1=1d0/yhelp1
           help=x*(-yhelp2+yhelp4)*yhelp1
           dby2=x*(-yhelp2-yhelp4)*yhelp1
           if(help<0d0) help=2d0*pc_c*tsp_dt*thelpinv
           if(dby2<0d0) dby2=2d0*pc_c*tsp_dt*thelpinv
           dby2=min(help,dby2)
           iynext2=iy+1
        endif
     endif
  endif
!  write(*,*) idby1,dby1,idby2,dby2
!  if(y<grd_yarr(iy+1).and.y>grd_yarr(iy)) write(*,*) idby1,idby2
!  if(dby1==0d0.and.idby1==4) write(*,*) '1: ',idby1, y, iy, dby1
!  if(dby2==0d0.and.idby2==4) write(*,*) '2: ',idby2, y, iy
  ! if(dby1==0d0.and.dby2==0d0) stop 'transport1: invalid dby[1,2]'
  if(dby1<0d0.and.dby2<0d0) then
     write(*,*) iy, y
     write(*,*) idby1, dby1, idby2, dby2
     stop 'transport1: dby1<0 and dby2<0'
  endif
  ! if(dby1==0d0) then
  !    write(*,*) '1: ', iy, y, dby2, eta
  ! endif
  ! if(dby2==0d0) then
  !    write(*,*) '2: ', iy, y, dby1, eta
  ! endif
  if(dby1<=0d0) dby1=2d0*pc_c*tsp_dt*thelpinv
  if(dby2<=0d0) dby2=2d0*pc_c*tsp_dt*thelpinv
  dby=min(dby1,dby2)
  if(dby==dby1) then
     iynext=iynext1
  else
     iynext=iynext2
  endif

!-- azimuthal boundary distance (z)
  if(xi==0d0.or.grd_nz==1) then
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  elseif(xi>0d0) then
!-- counterclockwise
     iznext=iz+1
     if(iznext==grd_nz+1) iznext=1
     zhelp = muy*cos(grd_zarr(iz+1))-mux*sin(grd_zarr(iz+1))
     if(zhelp==0d0) then
        dbz = 2d0*pc_c*tsp_dt*thelpinv
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz+1)-z)/zhelp
        if(dbz<0d0) dbz = 2d0*pc_c*tsp_dt*thelpinv
     endif
  else
!-- clockwise
     iznext=iz-1
     if(iznext==0) iznext=grd_nz
     zhelp = muy*cos(grd_zarr(iz))-mux*sin(grd_zarr(iz))
     if(zhelp==0d0) then
        dbz = 2d0*pc_c*tsp_dt*thelpinv
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz)-z)/zhelp
        if(dbz<0d0) dbz = 2d0*pc_c*tsp_dt*thelpinv
     endif
  endif

!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     r1 = rnd_r(rnd_state)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif

!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ic))
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif

!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = pc_c*(elabfact-wl*grp_wlinv(ig+1))
     if(ddop<0d0) then
        ddop = 2d0*pc_c*tsp_dt*thelpinv
     endif
  else
     ddop = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dcen,dbx,dby,dbz,dthm,dcol,ddop)
  if(d<0d0) then
     write(*,*) dcen,dbx,dby,dbz,dthm,dcol,ddop,ix,iy,iz,x,y,z,xi,eta,mu
     stop 'transport1: negative distance'
  endif

!-- storing old position
  xold = x
  yold = y
  muold = mu
!-- updating radius
  x = sqrt((1d0-mu**2)*x**2+(d+x*mu)**2)
  if(x<1d-15*grd_xarr(2).and.muold==-1d0) then
!-- sanity check
     if(d==dbx) stop 'transport1: x<1d-15*xarr(2),d==dbx,mu=-1'
!-- excluding dbz
     dbz = 2d0*pc_c*tsp_dt*thelpinv
!-- resetting direction
     mu = 1d0
     eta = 0d0
     xi = 0d0
  else
     if(x<1d-15*grd_xarr(2)) stop 'transport1: x=0 and muold/=-1'
!-- updating radial projection of direction
     mu = (xold*mu+d)/x
!-- updating polar projection of position
     y = (xold*yold+muz*d)/x
!-- updating azimuthal angle of position
     z = atan2(xold*sqrt(1d0-yold**2)*sin(z)+muy*d , &
          xold*sqrt(1d0-yold**2)*cos(z)+mux*d)
     if(z<0d0) z=z+pc_pi2
!-- updating azimuthal angle of direction (about radius)
     eta = y*(cos(z)*mux+sin(z)*muy)-sqrt(1d0-y**2)*muz
     xi = cos(z)*muy-sin(z)*mux
     om = atan2(xi,eta)
     if(om<0d0) om=om+pc_pi2
  endif

!-- direction sanity check
  if(abs(mux**2+muy**2+muz**2-1d0)>1d-9) then
     write(*,*) mux**2+muy**2+muz**2,mux,muy,muz
     stop 'transport1: invalid mux,muy,muz'
  endif
  if(abs(mu**2+eta**2+xi**2-1d0)>1d-9) then
     write(*,*) mu**2+eta**2+xi**2,mu,eta,xi
     stop 'transport1: invalid mu,eta,xi'
  endif

!-- updating time
  ptcl%t = ptcl%t + thelp*d*cinv

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     grd_eraddens(ic)=grd_eraddens(ic)+e*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)* &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
          thelp>1d-6) then
        grd_eraddens(ic) = grd_eraddens(ic)+e* &
             (1.0d0-exp(-grd_fcoef(ic)*elabfact* &
             grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*elabfact * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ic)=grd_eraddens(ic)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ic)=grd_edep(ic)+e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ic)/=grd_edep(ic)) then
        stop 'transport1: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
          elabfact*d*thelp)
  endif
!
!-- updating transformation factors
  if(grd_isvelocity) elabfact=1d0-mu*x*cinv

!
!-- census
  if(d==dcen) then
     prt_done = .true.
     grd_numcensus(ic)=grd_numcensus(ic)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rnd_r(rnd_state)
     mu = 1d0 - 2d0*r1
     r1 = rnd_r(rnd_state)
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+mu*x*cinv)
  elseif(d==dbx) then
!-- checking if escaped domain
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- ending particle
        isvacant = .true.
        prt_done = .true.
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
     endif
  endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- lab wavelength
        wl = wl*(1d0-mu*x*cinv)/elabfact
        help = elabfact/(1d0-mu*x*cinv)
!-- velocity effects accounting
        tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
        e = e*help
        e0 = e0*help
     endif

!
!-- effective collision
  elseif(d==dcol) then
     r1 = rnd_r(rnd_state)
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ic)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ic)=grd_edep(ic)+e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        r1 = rnd_r(rnd_state)
        if(grp_ng>1) ig = emitgroup(r1,ic)
!-- uniformly in new group
        r1 = rnd_r(rnd_state)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- transforming to lab
        if(grd_isvelocity) then
           wl = wl*(1d0-mu*x*cinv)
           help = elabfact/(1d0-mu*x*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!-- checking for DDMC in new group
        if((grd_cap(ig,ic)+grd_sig(ic)) * &
             min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
             thelp>=prt_tauddmc.and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*mu*x*cinv
!
              e = e*(1d0-mu*x*cinv)
              e0 = e0*(1d0-mu*x*cinv)
              wl = wl/(1d0-mu*x*cinv)
           endif
        endif
     endif

!
!-- radial bound
  elseif(d==dbx .and. dbx<dby) then

     if(mu>=0d0) then
        ihelp=ix+1
        x = grd_xarr(ix+1)
     else
        if(ix==1) stop 'transport1: ix=1 and mu<0'
        ihelp=ix-1
        x = grd_xarr(ix)
     endif

     l = grd_icell(ihelp,iy,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ihelp),xm(ihelp)*dyac(iy),xm(ihelp) * &
          ym(iy)*dz(iz))*thelp<prt_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        ix = ihelp
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming x-cosine to cmf
           mu = (mu-x*cinv)/elabfact
!-- amplification factor
           if(mu<0d0) then
              help = 1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              tot_evelo = tot_evelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              grd_eamp(ic) = grd_eamp(ic) + &
                   e*2d0*0.55d0*max(0d0,help-2d0)*x*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
              e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if(r1<help*(1d0+1.5d0*abs(mu))) then
           if(grd_isvelocity) then
              ptcl%itype = 2
              grd_methodswap(ic)=grd_methodswap(ic)+1
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ihelp
           ic = grd_icell(ix,iy,iz)
        else
!-- resampling x-cosine
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
           mu=(ix-ihelp)*max(r1,r2)
!-- resampling azimuthal
           r1 = rnd_r(rnd_state)
           om = pc_pi2*r1
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
        endif
     endif

!
!-- polar bound
  elseif(d==dby) then

!-- iznext=iz except
     iznext=iz
     if(x<1d-15*grd_xarr(2).and.muold==-1d0) then
!-- reflecting y
        y=-y
        iynext=binsrch(y,grd_yarr,grd_ny+1)
!-- reflecting z
        z=z+pc_pi
        if(z>pc_pi2) z=z-pc_pi2
        if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1)
     elseif(iynext==iy-1) then
        if(abs(y-grd_yarr(iy))>1d-9) then
           write(*,*) iy,'y: ',y,'yarr(iy): ',grd_yarr(iy)
           stop 'transport1: y/=yarr(iy)'
        endif
        y=grd_yarr(iy)
        if(iynext==0) then
!-- reflecting z
           z=z+pc_pi
           if(z>pc_pi2) z=z-pc_pi2
           if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1)
           iynext=1
        endif
     elseif(iynext==iy+1) then
        if(abs(y-grd_yarr(iy+1))>1d-9) then
           write(*,*) iy,'y: ',y,'yarr(iy+1): ',grd_yarr(iy+1)
           stop 'transport1: y/=yarr(iy+1)'
        endif
        y=grd_yarr(iy+1)
        if(iynext==grd_ny+1) then
!-- reflecting z
           z=z+pc_pi
           if(z>pc_pi2) z=z-pc_pi2
           if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1)
           iynext=grd_ny
        endif
     else
!-- sanity check
        write(*,*) dby
        write(*,*) y,grd_yarr(iy),grd_yarr(iy+1),iy,iynext
        stop 'transport1: invalid polar bound crossing'
     endif

     l = grd_icell(ix,iynext,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iynext),xm(ix)*ym(iynext) * &
          dz(iznext))*thelp<prt_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        iyold=iy
        iy = iynext
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming y-cosine to cmf
           mu=(mu-x*cinv)/(1d0-x*mu*cinv)
           eta = sqrt(1d0-mu**2)*cos(om)
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             xm(ix)*dyac(iynext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if(r1 < help*(1d0+1.5d0*abs(eta))) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iynext
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
           eta = (iynext-iy)*max(r1,r2)
           r1 = rnd_r(rnd_state)
           xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
!-- resampling x-cosine
           mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(xi,eta)
           if(om<0d0) om=om+pc_pi2
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
        endif
     endif

!
!-- azimuthal bound
  elseif(d==dbz) then
!-- sanity check
     if(grd_nz==1) stop 'transport1: invalid z crossing'
     if(iznext==grd_nz.and.iz==1) then
        z=pc_pi2
     elseif(iznext==1.and.iz==grd_nz) then
        z=0d0
     elseif(iznext==iz-1) then
        z=grd_zarr(iz)
     else
!-- iznext==iz+1
        if(iznext/=iz+1) stop 'transport1: invalid iznext'
        z=grd_zarr(iz+1)
     endif
     l = grd_icell(ix,iy,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy) * &
          dz(iznext))*thelp<prt_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming z-cosine to cmf
           mu = (mu-x*cinv)/elabfact
           xi = sqrt(1d0-mu**2)*sin(om)
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             xm(ix)*ym(iy)*dz(iznext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(xi))) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
!-- resampling z-cosine
           if(iznext==iz+1.or.(iznext==1.and.iz==grd_nz)) then
              xi = -max(r1,r2)
           else
              xi = max(r1,r2)
           endif
           r1 = rnd_r(rnd_state)
           eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
!-- resampling x-cosine
           mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(xi,eta)
           if(om<0d0) om=om+pc_pi2
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
!-- reverting z
           if(iznext==grd_nz.and.iz==1) then
              z=0d0
           elseif(iznext==1.and.iz==grd_nz) then
              z=pc_pi2
           endif
        endif
     endif

!
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) stop 'transport1: ddop and no velocity'
     if(ig<grp_ng) then
!-- shifting group
        ig = ig+1
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        r1 = rnd_r(rnd_state)
        wl = 1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if((grd_cap(ig,ic)+grd_sig(ic)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
          thelp>=prt_tauddmc.and..not.in_puretran) then
        ptcl%itype = 2
        grd_methodswap(ic)=grd_methodswap(ic)+1
        if(grd_isvelocity) then
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  else
     stop 'transport1: invalid distance'
  endif
  idby1old=idby1
  idby2old=idby2
  dby1old=dby1
  dby2old=dby2
  muzold=muz
  etaold=eta
  dbyold=dby

  if((y>grd_yarr(iy+1) .or. y<grd_yarr(iy))) then
     write(0,*) 'theta not in cell (x): ',ix,xold,x,grd_xarr(ix),grd_xarr(ix+1)
     write(0,*) 'theta not in cell (y): ',iy,yold,y,grd_yarr(iy),grd_yarr(iy+1)
     write(0,*) 'old (y): ',iyold,dbyold
     write(0,*) 'dir: ',mux,muy,muz,mu,eta,xi
     write(0,*) 'dby: ',dby,'dbx: ',dbx,'max d: ',2d0*pc_c*tsp_dt*thelpinv
     write(0,*) 'dby1: ',dby1,'dby2: ',dby2,'etaold: ',etaold
     write(0,*) 'idby1: ',idby1,'idby2: ',idby2
     write(0,*) d
     write(0,*)
  endif

end subroutine transport1
