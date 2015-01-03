subroutine transport1_gamgrey(ptcl,ic)

  use randommod
  use miscmod
  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  integer,intent(inout) :: ic
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time
!
  integer :: imu, iom, ihelp
  real*8 :: elabfact, eta, xi, mux,muy,muz
  real*8 :: r1, thelp,thelpinv, help
  real*8 :: dcol,dbx,dby,dbz,d

  integer :: iynext,iznext
  real*8 :: yhelp1,yhelp2,dby1,dby2
  real*8 :: zhelp1,zhelp2,dbz1,dbz2
  real*8 :: xold,yold,muold

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0

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

!-- spherical projections
  eta = sqrt(1d0-mu**2)*cos(om)
  xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
  mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
  muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
  muz = mu*y-eta*sqrt(1d0-y**2)
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     elabfact = 1d0-mu*x*cinv
     thelp = tsp_t
  else
     elabfact = 1d0
     thelp = 1d0
  endif
  thelpinv = 1d0/thelp

!-- radial boundary distance (x)
  if (ix==1) then
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
  if(muz>=0d0) then
     iynext = iy+1
  else
     iynext = iy-1
  endif
  ihelp = max(iy,iynext)
  yhelp1 = y**2-(1d0-mu**2)*grd_yarr(ihelp)**2-2d0*muz*mu*y+muz**2
  if(yhelp1<0d0) then
     dby = 2d0*pc_c*tsp_dt*thelpinv
  else
     yhelp1 = sqrt(yhelp1)
     yhelp2 = grd_yarr(ihelp)**2-muz**2
     if(yhelp2==0d0) then
        dby = 2d0*pc_c*tsp_dt*thelpinv
     else
        yhelp2 = 1d0/yhelp2
        dby1 = x*(muz*y-mu*grd_yarr(ihelp)**2-grd_yarr(ihelp)*yhelp1)*yhelp2
        dby2 = x*(muz*y-mu*grd_yarr(ihelp)**2+grd_yarr(ihelp)*yhelp1)*yhelp2
        if(dby1<0d0) dby1 = 2d0*pc_c*tsp_dt*thelpinv
        if(dby2<0d0) dby2 = 2d0*pc_c*tsp_dt*thelpinv
        dby = min(dby1,dby2)
     endif
  endif

!-- azimuthal boundary distance (z)
  if(grd_nz==1) then
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  else
     zhelp1 = muy*cos(grd_zarr(iz))-mux*sin(grd_zarr(iz))
     if(zhelp1==0d0) then
        dbz1 = 2d0*pc_c*tsp_dt*thelpinv
     else
        zhelp1 = 1d0/zhelp1
!-- dbz1: iz->iz-1
        dbz1 = x*sqrt(1d0-y**2)*sin(grd_zarr(iz)-z)*zhelp1
        if(dbz1<0d0) dbz1 = 2d0*pc_c*tsp_dt*thelpinv
     endif
     zhelp2 = muy*cos(grd_zarr(iz+1))-mux*sin(grd_zarr(iz+1))
     if(zhelp2==0d0) then
        dbz2 = 2d0*pc_c*tsp_dt*thelpinv
     else
        zhelp2 = 1d0/zhelp2
!-- dbz1: iz->iz+1
        dbz2 = x*sqrt(1d0-y**2)*sin(grd_zarr(iz+1)-z)*zhelp2
        if(dbz2<0d0) dbz2 = 2d0*pc_c*tsp_dt*thelpinv
     endif
     dbz = min(dbz1,dbz2)
     if(dbz==dbz1) then
        iznext = iz-1
        if(iznext==0) iznext = grd_nz
     else
        iznext = iz+1
        if(iznext==grd_nz+1) iznext = 1
     endif
  endif

!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(grd_capgrey(ic)>0d0) then
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        dcol = -log(r1)*thelpinv/(grd_capgrey(ic)*elabfact)
     else
        dcol = 2d0*abs(pc_c*dt*thelpinv) !> dcen
     endif
  else
     dcol = 2d0*abs(pc_c*dt*thelpinv) !> dcen
  endif

!
!-- finding minimum distance
  d = min(dbx,dby,dbz,dcol)
  if(d<0d0) then
     write(*,*) dbx,dby,dbz,dcol,ix,iy,iz,x,y,z,xi,eta,mu
     stop 'transport1_gamgrey: negative distance'
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
  if(z<0d0) z = z+pc_pi2
!-- updating azimuthal angle of direction (about radius)
  eta = y*(cos(z)*mux+sin(z)*muy)-sqrt(1d0-y**2)*muz
  xi = cos(z)*muy-sin(z)*mux
  om = atan2(xi,eta)
  if(om<0d0) om = om+pc_pi2

!-- depositing nonanalog absorbed energy
  if(.not.prt_isimcanlog) then
     grd_edep(ic) = grd_edep(ic)+e* &
          (1d0-exp(-grd_capgrey(ic)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ic)/=grd_edep(ic)) then
        stop 'transport1_gamgrey: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgrey(ic)*elabfact*d*thelp)
  endif
!
!-- updating transformation factors
  if(grd_isvelocity) elabfact = 1d0-mu*x*cinv

  if(d==dcol) then
!-- sanity check
     if(.not.prt_isimcanlog) stop &
          'transport1_gamgrey: not isimcanlog and dcol<db[xyz]'
     prt_done = .true.
     grd_edep(ic) = grd_edep(ic) + e*elabfact
  elseif(d==dbx) then
     if(mu>=0d0) then
        if(ix==grd_nx) then
!-- ending particle
           prt_done = .true.
           help = atan2(muy,mux)
           if(help<0d0) help = help+pc_pi2
!-- retrieving lab frame flux group, polar, azimuthal bin
           iom = binsrch(help,flx_om,flx_nom+1)
           imu = binsrch(muz,flx_mu,flx_nmu+1)
           flx_gamluminos(imu,iom) = flx_gamluminos(imu,iom)+e/tsp_dt
           flx_gamlumdev(imu,iom) = flx_gamlumdev(imu,iom)+(e/tsp_dt)**2
           flx_gamlumnum(imu,iom) = flx_gamlumnum(imu,iom)+1
           return
        else
           x = grd_xarr(ix+1)
           ix = ix+1
        endif
     else
        x = grd_xarr(ix)
        ix = ix-1
     endif
     ic = grd_icell(ix,iy,iz)
  elseif(d==dby) then
     if(iynext==iy-1) then
        y = grd_yarr(iy)
     elseif(iynext==iy+1) then
        y = grd_yarr(iy+1)
     else
!-- sanity check
        stop 'transport1_gamgrey: invalid polar bound crossing'
     endif
     iy = iynext
     ic = grd_icell(ix,iy,iz)
  else !d==dbz
!-- sanity check
     if(grd_nz==1) stop 'transport1_gamgrey: invalid z crossing'
     if(iznext==grd_nz.and.iz==1) then
        z = pc_pi2
     elseif(iznext==1.and.iz==grd_nz) then
        z = 0d0
     elseif(iznext==iz-1) then
        z = grd_zarr(iz)
     else
!-- iznext==iz+1
        if(iznext/=iz+1) stop 'transport1_gamgrey: invalid iznext'
        z = grd_zarr(iz+1)
     endif
     iz = iznext
     ic = grd_icell(ix,iy,iz)
  endif

end subroutine transport1_gamgrey
