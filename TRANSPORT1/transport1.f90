subroutine transport1(ptcl,ig,isvacant)

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
  integer,intent(inout) :: ig
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

  integer :: iynext,iznext
  real*8 :: yhelp1,yhelp2,yhelp3,dby1,dby2
  real*8 :: zhelp1,zhelp2,dbz1,dbz2
  real*8 :: xold,yold,muold

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
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
!-- shortcut
  dtinv = 1d0/tsp_dt
!-- spherical projections
  eta = sqrt(1d0-mu**2)*sin(om)
  xi = sqrt(1d0-mu**2)*cos(om)
!-- planar projections (invariant until collision)
  mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
  muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
  muz = mu*y-eta*sqrt(1d0-y**2)
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
  yhelp1 = y**2-(1d0-mu**2)*grd_yarr(iy)**2-2d0*muz*mu*y+muz**2
  yhelp2 = y**2-(1d0-mu**2)*grd_yarr(iy+1)**2-2d0*muz*mu*y+muz**2
  if(yhelp1<0d0.and.yhelp2<0d0) then
     dby = 2d0*pc_c*tsp_dt*thelpinv
  elseif(yhelp1>=0d0) then
!-- iy->iy-1
     yhelp1 = sqrt(yhelp1)
     yhelp3 = grd_yarr(iy)**2-muz**2
     if(yhelp3==0d0) then
        dby = 2d0*pc_c*tsp_dt*thelpinv
     else
        yhelp3=1d0/yhelp3
        dby1 = x*(muz*y-mu*grd_yarr(iy)**2-grd_yarr(iy)*yhelp1)*yhelp3
        dby2 = x*(muz*y-mu*grd_yarr(iy)**2+grd_yarr(iy)*yhelp1)*yhelp3
        if(dby1<0d0) dby1=2d0*pc_c*tsp_dt*thelpinv
        if(dby2<0d0) dby2=2d0*pc_c*tsp_dt*thelpinv
        iynext = iy-1
     endif
  else
!-- iy->iy+1
     yhelp2 = sqrt(yhelp2)
     yhelp3 = grd_yarr(iy)**2-muz**2
     if(yhelp3==0d0) then
        dby = 2d0*pc_c*tsp_dt*thelpinv
     else
        yhelp3=1d0/yhelp3
        dby1 = x*(muz*y-mu*grd_yarr(iy+1)**2-grd_yarr(iy+1)*yhelp2)*yhelp3
        dby2 = x*(muz*y-mu*grd_yarr(iy+1)**2+grd_yarr(iy+1)*yhelp2)*yhelp3
        if(dby1<0d0) dby1=2d0*pc_c*tsp_dt*thelpinv
        if(dby2<0d0) dby2=2d0*pc_c*tsp_dt*thelpinv
        dby = min(dby1,dby2)
        iynext = iy+1
     endif
  endif

!-- azimuthal boundary distance (z)
  zhelp1 = muy*cos(grd_zarr(iz))-mux*sin(grd_zarr(iz))
  if(zhelp1==0d0) then
     dbz1=2d0*pc_c*tsp_dt*thelpinv
  else
     zhelp1=1d0/zhelp1
!-- dbz1: iz->iz-1
     dbz1 = x*sqrt(1d0-y**2)*sin(grd_zarr(iz)-z)*zhelp1
     if(dbz1<0d0) dbz1=2d0*pc_c*tsp_dt*thelpinv
  endif
  zhelp2 = muy*cos(grd_zarr(iz+1))-mux*sin(grd_zarr(iz+1))
  if(zhelp2==0d0) then
     dbz2=2d0*pc_c*tsp_dt*thelpinv
  else
     zhelp2=1d0/zhelp2
!-- dbz1: iz->iz+1
     dbz2 = x*sqrt(1d0-y**2)*sin(grd_zarr(iz+1)-z)*zhelp2
     if(dbz2<0d0) dbz2=2d0*pc_c*tsp_dt*thelpinv
  endif
  dbz = min(dbz1,dbz2)
  if(dbz==dbz1) then
     iznext=iz-1
     if(iznext==0) iznext=grd_nz
  else
     iznext=iz+1
     if(iznext==grd_nz+1) iznext=1
  endif

!-- Thomson scattering distance
  if(grd_sig(ix,iy,iz)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ix,iy,iz))
  else
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif

!-- effective collision distance
  if(grd_cap(ig,ix,iy,iz)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ix,iy,iz))
  elseif(grd_fcoef(ix,iy,iz)<1d0.and.grd_fcoef(ix,iy,iz)>=0d0) then
     r1 = rand()
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ix,iy,iz))*grd_cap(ig,ix,iy,iz))
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

!-- updating radius
  xold = x
  x = sqrt((1d0-mu**2)*x**2+(d+x*mu)**2)
!-- updating radial projection of direction
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (xold*mu+d)/x
  endif
!-- updating polar projection of position
  yold = y
  if(x/=0d0) y = (xold*yold+muz*d)/x
!-- updating azimuthal angle of position
  z = atan2(xold*sqrt(1d0-yold**2)*sin(z)+muy*d , &
       xold*sqrt(1d0-yold**2)*cos(z)+mux*d)
  if(z<0d0) z=z+pc_pi2
!-- updating azimuthal angle of direction (about radius)
  eta = y*(cos(z)*mux+sin(z)*muy)-sqrt(1d0-y**2)*muz
  xi = cos(z)*muy-sin(z)*mux
  om = atan2(xi,eta)
  if(om<0d0) om=om+pc_pi2
!-- updating time
  ptcl%t = ptcl%t + thelp*d*cinv

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     grd_eraddens(ix,iy,iz)=grd_eraddens(ix,iy,iz)+e*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz)* &
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        grd_eraddens(ix,iy,iz) = grd_eraddens(ix,iy,iz)+e* &
             (1.0d0-exp(-grd_fcoef(ix,iy,iz)*elabfact* &
             grd_cap(ig,ix,iy,iz)*d*thelp))* &
             elabfact/(grd_fcoef(ix,iy,iz)*elabfact * &
             grd_cap(ig,ix,iy,iz)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ix,iy,iz)=grd_eraddens(ix,iy,iz)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e* &
          (1d0-exp(-grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ix,iy,iz)/=grd_edep(ix,iy,iz)) then
        stop 'transport1: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz) * &
          elabfact*d*thelp)
  endif
!
!-- updating transformation factors
  if(grd_isvelocity) elabfact=1d0-mu*x*cinv

!
!-- census
  if(d==dcen) then
     prt_done = .true.
     grd_numcensus(ix,iy,iz)=grd_numcensus(ix,iy,iz)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rand()
     mu = 1d0 - 2d0*r1
     r1 = rand()
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
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ix,iy,iz)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        r1 = rand()
        ig = emitgroup(r1,ix,iy,iz)
!-- uniformly in new group
        r1 = rand()
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
        if((grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz)) * &
             min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
             thelp>=prt_tauddmc.and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
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
  elseif(d==dbx) then

!
!-- polar bound
  elseif(d==dby) then

!
!-- azimuthal bound
  elseif(d==dbz) then

!
!-- Doppler shift
  elseif(d==ddop) then

  else
     stop 'transport1: invalid distance'
  endif

end subroutine transport1
