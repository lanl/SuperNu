subroutine transport1_gamgrey(ptcl,ptcl2)

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
  type(packet2),target,intent(inout) :: ptcl2
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time
!
  integer :: imu, iom
  real*8 :: elabfact, eta, xi, mux,muy,muz
  real*8 :: r1, thelp,thelpinv, help
  real*8 :: dcol,dbx,dby,dbz,d
  real*8 :: darr(4)

  integer :: iynext,iynext1,iynext2,iznext,idby1,idby2
  real*8 :: yhelp1,yhelp2,yhelp3,yhelp4,dby1,dby2
  real*8 :: zhelp
  real*8 :: xold,yold,zold,muold
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix,iy,iz,ic,ig
  real*8,pointer :: x,y,z,mu,om,e,e0

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

!-- spherical projections
  eta = sqrt(1d0-mu**2)*cos(om)
  xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
  mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
  muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
  muz = mu*y-eta*sqrt(1d0-y**2)
!-- direction sanity check
  if(abs(mux**2+muy**2+muz**2-1d0)>1d-9) then
     write(*,*) 'transport1_gamgrey: invalid mux,muy,muz', &
        mux**2+muy**2+muz**2-1d0,mux,muy,muz
  endif
!-- normalize direction
  help = 1d0/sqrt(mux**2+muy**2+muz**2)
  mux = mux*help
  muy = muy*help
  muz = muz*help
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

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*dt*thelpinv) !> dcen

!-- radial boundary distance (x)
  if (ix==1) then
     dbx = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  elseif (mu < -sqrt(1d0-(grd_xarr(ix)/x)**2)) then
     dbx = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
  else
     dbx = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  endif
!-- sanity check
  if(dbx/=dbx) stop 'transport1_gamgrey: dbx/=dbx'
!
!-- polar boundary distance (y)
  yhelp1=grd_yarr(iy)**2-muz**2
  yhelp2=mu*grd_yarr(iy)**2-muz*y
  if(y==grd_yarr(iy+1).and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
     yhelp3=0d0
  else
     yhelp3=grd_yarr(iy)**2-y**2
  endif
  if(yhelp1==0d0.and.yhelp3==0d0) then
     idby1=1
!-- particle, direction on cone
     dby1 = far
     iynext1=iy
  elseif(yhelp1==0d0) then
     if((muz>=0d0.and.muz==grd_yarr(iy)).or. &
          (muz>=0d0.and.muz==-grd_yarr(iy)).or. &
          yhelp2==0d0) then
        idby1=2
!-- direction parallel to lower cone or nonphysical
        dby1 = far
        iynext1=iy
     else
        idby1=3
        dby1=-0.5*x*yhelp3/yhelp2
        iynext1=iy-1
     endif
  else
     yhelp4=yhelp2**2-yhelp1*yhelp3
     if(abs(yhelp4)<1d-12*abs(yhelp2)) yhelp4=0d0
     if(yhelp4<0d0) then
        idby1=4
!-- not intersecting lower cone
        dby1 = far
        iynext1=iy
     else
!-- intersecting lower cone at at least one point
        if(cos(om)<0d0.and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
!-- choose dby2
           idby1=5
           dby1 = far
           iynext1=iy
        else
!write(0,*) yhelp4,yhelp1,yhelp2,yhelp3
           if(yhelp3==0d0) then
              idby1=6
              yhelp4=abs(yhelp2)
           else
              idby1=7
              yhelp4=sqrt(yhelp4)
           endif
           yhelp1=1d0/yhelp1
           help=x*(-yhelp2+yhelp4)*yhelp1
           dby1=x*(-yhelp2-yhelp4)*yhelp1
           if(help<=0d0) help=far
           if(dby1<=0d0) dby1=far
           dby1=min(help,dby1)
           iynext1=iy-1
        endif
     endif
  endif

!-- dby2: iy->iy+1
  yhelp1=grd_yarr(iy+1)**2-muz**2
  yhelp2=mu*grd_yarr(iy+1)**2-muz*y
  if(y==grd_yarr(iy).and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
     yhelp3=0d0
  else
     yhelp3=grd_yarr(iy+1)**2-y**2
  endif
  if(yhelp1==0d0.and.yhelp3==0d0) then
     idby2=1
!-- particle, direction on cone
     dby2 = far
     iynext2=iy
  elseif(yhelp1==0d0) then
     if((muz<=0d0.and.muz==-grd_yarr(iy+1)).or. &
          (muz<=0d0.and.muz==grd_yarr(iy+1)).or. &
          yhelp2==0d0) then
        idby2=2
!-- direction parallel to upper cone or nonphysical
        dby2 = far
        iynext2=iy
     else
        idby2=3
        dby2=-0.5*x*yhelp3/yhelp2
        iynext2=iy+1
     endif
  else
     yhelp4=yhelp2**2-yhelp1*yhelp3
     if(abs(yhelp4)<1d-12*abs(yhelp2)) yhelp4=0d0
     if(yhelp4<0d0) then
        idby2=4
!-- not intersecting upper cone
        dby2 = far
        iynext2=iy
     else
!-- intersecting upper cone at at least one point
        if(cos(om)>=0d0.and.abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9) then
!-- choose dby1
!-- choose dby1
           idby2=5
           dby2 = far
           iynext2=iy
        else
           if(yhelp3==0d0) then
              idby2=6
              yhelp4=abs(yhelp2)
           else
              idby2=7
              yhelp4=sqrt(yhelp4)
           endif
           yhelp1=1d0/yhelp1
           help=x*(-yhelp2+yhelp4)*yhelp1
           dby2=x*(-yhelp2-yhelp4)*yhelp1
           if(help<=0d0) help=far
           if(dby2<=0d0) dby2=far
           dby2=min(help,dby2)
           iynext2=iy+1
        endif
     endif
  endif
  if(dby1<0d0.and.dby2<0d0) then
     write(*,*) iy, y
     write(*,*) idby1, dby1, idby2, dby2
     stop 'transport1_gg: dby1<0 and dby2<0'
  endif
  if(dby1<=0d0) dby1=far
  if(dby2<=0d0) dby2=far
  dby=min(dby1,dby2)
  if(dby==dby1) then
     iynext=iynext1
  else
     iynext=iynext2
  endif

!-- azimuthal boundary distance (z)
  iznext = iz
  if(xi==0d0.or.grd_nz==1) then
     dbz = far
  elseif(abs(y)==1d0) then
!-- azimuthal position undetermined
     dbz = far
  elseif(xi>0d0.and.grd_zarr(iz+1)-pc_pi<z) then
!-- counterclockwise
     iznext=iz+1
     zhelp = muy*cos(grd_zarr(iz+1))-mux*sin(grd_zarr(iz+1))
     if(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz+1)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  elseif(xi<0d0.and.z<pc_pi+grd_zarr(iz)) then
!-- clockwise
     iznext=iz-1
     zhelp = muy*cos(grd_zarr(iz))-mux*sin(grd_zarr(iz))
     if(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  else
     dbz = far
  endif

!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(grd_capgrey(ic)>0d0) then
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        dcol = -log(r1)*thelpinv/(grd_capgrey(ic)*elabfact)
     else
        dcol = far
     endif
  else
     dcol = far
  endif

!
!-- finding minimum distance
  darr = [dbx,dby,dbz,dcol]
  if(any(darr/=darr) .or. any(darr<0d0)) then
     write(0,*) darr
     write(0,*) ix,iy,iz,x,y,z,mu,eta,xi,om
     stop 'transport1_gamgrey: invalid distance'
  endif
  d = minval(darr)

!-- storing old position
  xold = x
  yold = y
  zold = z
  muold = mu
!-- updating radius
  x = sqrt((1d0-mu**2)*x**2+(d+x*mu)**2)
  if(x/=x) stop 'transport1: x/=x'
  if(x<1d-15*grd_xarr(2).and.muold==-1d0) then
!-- sanity check
     if(d==dbx) stop 'transport1_gamgrey: x<1d-15*xarr(2),d==dbx,mu=-1'
!-- excluding dbz
     dbz = far
!-- resetting direction
     mu = 1d0
     eta = 0d0
     xi = 0d0
  else
     if(x<1d-15*grd_xarr(2)) stop &
          'transport1_gamgrey: x=0 and muold/=-1'
!-- updating radial projection of direction
     mu = (xold*mu+d)/x
     mu = max(mu,-1d0)
     mu = min(mu,1d0)
!-- updating polar projection of position
     y = (xold*yold+muz*d)/x
     y = max(y,-1d0)
     y = min(y,1d0)
!-- updating azimuthal angle of position
     z = atan2(xold*sqrt(1d0-yold**2)*sin(z)+muy*d , &
          xold*sqrt(1d0-yold**2)*cos(z)+mux*d)
     if(d==dbz .and. abs(z)<1d-9.and.iz==1) then
        z = 0d0
     elseif(d==dbz .and. abs(z)<1d-9.and.iz==grd_nz) then
        z = pc_pi2
     elseif(z<0d0) then
        z=z+pc_pi2
     endif
!-- any iz possible if previous position was undetermined
     if(abs(yold)==1d0) then
        iz = binsrch(z,grd_zarr,grd_nz+1,.false.)
        ic = grd_icell(ix,iy,iz)
     endif
!-- updating azimuthal angle of direction (about radius)
     eta = y*(cos(z)*mux+sin(z)*muy)-sqrt(1d0-y**2)*muz
     xi = cos(z)*muy-sin(z)*mux
     om = atan2(xi,eta)
     if(om<0d0) om = om+pc_pi2
!-- direction sanity check
     if(abs(mu**2+eta**2+xi**2-1d0)>1d-9) then
        write(*,*) 'transport1_gamgrey: invalid mu,eta,xi',mu**2+eta**2+xi**2-1d0,mu,eta,xi
     endif
!-- normalize direction
     help = 1d0/sqrt(mu**2+eta**2+xi**2)
     mu = mu*help
     eta = eta*help
     xi = xi*help
  endif

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
     ptcl2%done = .true.
     grd_edep(ic) = grd_edep(ic) + e*elabfact
  elseif(d==dbx .and. dbx<dby) then
     if(mu>=0d0) then
        if(ix==grd_nx) then
!-- ending particle
           ptcl2%done = .true.
           help = atan2(muy,mux)
           if(help<0d0) help = help+pc_pi2
!-- retrieving lab frame flux group, polar, azimuthal bin
           iom = binsrch(help,flx_om,flx_nom+1,.false.)
           imu = binsrch(muz,flx_mu,flx_nmu+1,.false.)
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

!-- iznext=iz except
     iznext=iz
     if(x<1d-15*grd_xarr(2).and.muold==-1d0) then
!-- reflecting y
        y=-y
        iynext=binsrch(y,grd_yarr,grd_ny+1,.false.)
!-- reflecting z
        z=z+pc_pi !z is not updated with atan2 calculation
        if(z>pc_pi2) z=z-pc_pi2
        if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1,.false.)
     elseif(iynext==iy-1) then
        if(abs(y-grd_yarr(iy))>1d-9) then
           write(*,*) 'transport1_gamgrey: y/=yarr(iy)', iy,y,grd_yarr(iy)
        endif
        y=grd_yarr(iy)
        if(iynext==0) then
!-- reflecting z
           z=zold+pc_pi
           if(z>pc_pi2) z=z-pc_pi2
           if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1,.false.)
           iynext=1
        endif
     elseif(iynext==iy+1) then
        if(abs(y-grd_yarr(iy+1))>1d-9) then
           write(*,*) 'transport1_gamgrey: y/=yarr(iy+1)',iy,y,grd_yarr(iy+1)
        endif
        y=grd_yarr(iy+1)
        if(iynext==grd_ny+1) then
!-- reflecting z
           z=zold+pc_pi
           if(z>pc_pi2) z=z-pc_pi2
           if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1,.false.)
           iynext=grd_ny
        endif
     else
!-- sanity check
        write(*,*) dby,dby1,dby2,idby1,idby2,ptcl2%ipart,ptcl2%istep
        write(*,*) x,y,mu,grd_yarr(iy),grd_yarr(iy+1),iy,iynext
        stop 'transport1_gamgrey: invalid polar bound crossing'
     endif

     iz = iznext
     iy = iynext
     ic = grd_icell(ix,iy,iz)
  elseif(d==dbz) then
!-- sanity check
     if(grd_nz==1) stop 'transport1_gamgrey: invalid z crossing'
     if(iznext==iz-1) then
        z=grd_zarr(iz)
        if(iznext==0) then
           iznext = grd_nz
           z = pc_pi2
        endif
     elseif(iznext==iz+1) then
        z=grd_zarr(iz+1)
        if(iznext==grd_nz+1) then
           iznext = 1
           z = 0d0
        endif
     else
        stop 'transport1: invalid iznext'
     endif
     iz = iznext
     ic = grd_icell(ix,iy,iz)
  else
     stop 'transport1_gamgrey: invalid distance mode'
  endif

end subroutine transport1_gamgrey
