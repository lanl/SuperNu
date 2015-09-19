pure subroutine transport2_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)

  use randommod
  use miscmod
  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time

  logical :: lredir !direction resampled
  logical :: loutx,louty
  integer :: ixnext,iynext
  real*8 :: elabfact, dirdotu, gm
  real*8,pointer :: mux,muy,muz
  real*8 :: thelp, thelpinv
  real*8 :: dcol,dbx,dby
  real*8 :: darr(3)
  real*8 :: xold, omold, zold
  real*8 :: r1
!-- distance out of physical reach
  real*8 :: far

  real*8 :: mux1,mux2
  real*8 :: angrat1,angrat2,angrat3

  integer,pointer :: ix, iy, ic, ig
  integer,parameter :: iz=1
  real*8,pointer :: x,y,z,mu,om,e,d
!-- statement functions
  integer :: l
  real*8 :: muxf,zz
  muxf(zz) = ptcl%x*sin(ptcl2%muz+zz)/sin(ptcl2%muz)
  !muyf(zz) = ptcl%x*sin(zz)/sin(ptcl2%muz)

  ix => ptcl2%ix
  iy => ptcl2%iy
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e

  mux => ptcl2%mux
  muy => ptcl2%muy
  muz => ptcl2%muz

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0

!-- direction resample flag
  lredir = .false.
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
     thelp = tsp_t
  else
     dirdotu = 0d0
     elabfact = 1d0
     thelp = 1d0
  endif
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*dt*thelpinv) !> dcen
!
!-- boundary distances
!-- to x-bound
  if(abs(mu)==1d0) then
!-- making greater than dcen
     dbx = far
  else
     if(abs(sin(om))*x<grd_xarr(ix) .and. cos(om)<0d0 .and. ix/=1) then
!-- inner boundary
        dbx = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
        ixnext = ix-1
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbx = 0d0
        ixnext = ix+1
     else
!-- outer boundary
        dbx = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
        ixnext = ix+1
     endif
  endif
  if(dbx/=dbx) then
!    stop 'transport2_gamgrey: dbx nan'
     ierr = 1
     return
  endif

!-- to y-bound
  if(mu>0d0) then
     dby = (grd_yarr(iy+1)-y)/mu
     iynext = iy+1
  elseif(mu<0d0) then
     dby = (grd_yarr(iy)-y)/mu
     iynext = iy-1
  else
!-- making greater than dcen
     dby = far
  endif
!
!-- effective collision distance
  if(grd_capgam(ic)<=0d0) then
!-- making greater than dcen
     dcol = far
  elseif(trn_isimcanlog) then
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ic))
  else
!-- making greater than dcen
     dcol = far
  endif
!
!-- finding minimum distance
  darr = [dcol,dbx,dby]
  ptcl2%idist = minloc(darr,dim=1)
  d = minval(darr)
  if(any(darr/=darr)) then
     ierr = 4
     return
  endif
  if(d<0d0) then
     ierr = 5
     return
  endif
!
!-- store prev position
  xold = x
  zold = z
  omold = om
!
!-- update distance to intercept
  muy = muy + d*sqrt(1d0-mu**2)
!
!-- update position
  if(d>0d0) then
     y = y + mu*d!{{{

     if(abs(abs(cos(muz))-1d0)<1d-2) then
!-- muz calculation unreliable
        x = sqrt(xold**2 + (1d0-mu**2)*d**2 + 2d0*xold*sqrt(1d0-mu**2)*d*cos(omold))
     else
        x = sqrt(mux**2 + muy**2 - 2d0*mux*muy*cos(muz))
     endif
!
!-- special case: particle close to corners
     if(abs(x-xold)<1d-15*(x+xold)) then
        if(xold==grd_xarr(ix)) x = xold
        if(xold==grd_xarr(ix+1)) x = xold
     endif
!
!-- update direction
     if(x==0) then
!-- todo: implement exception for x==0
        ierr = 8
        return
     endif

!-- trigoniometric ratios
     angrat1 = (x**2 + mux**2 - muy**2)/(2*x*mux)
     angrat2 = sin(muz)*muy/x

     if(abs(angrat1)<1d0 .and. abs(angrat1)<abs(angrat2)) then
!-- method 1: calculate z from invariants
        z = acos(angrat1)
!-- check cos flip
        mux1 = muxf(z)
        mux2 = muxf(pc_pi2-z)
        if(abs(mux1-mux) > abs(mux2-mux)) z = pc_pi2-z
     elseif(abs(angrat2)<1d0) then
!-- method 2: calculate z from invariants
        z = asin(angrat2)
!-- check sin flip
        mux1 = muxf(z)
        mux2 = muxf(pc_pi-z)
        if(abs(mux1-mux) > abs(mux2-mux)) z = pc_pi-z
     else
!-- method 3: calculate om from xold and omold
        angrat3 = (sqrt(1d0-mu**2)*d + xold*cos(omold))/x
        om = acos(angrat3)
!-- check cos flip
        mux1 = muxf(zold+omold-om)
        mux2 = muxf(-muz-pc_pi+om)
        if(abs(mux1-mux) > abs(mux2-mux)) om = pc_pi2-om
        z = zold - (om - omold)
        if(z>pc_pi2) z = z - pc_pi2
     endif

     if(z<0d0) z = z + pc_pi2

!-- update om
     om = omold - (z - zold)
     if(om<0d0) om = om+pc_pi2
     if(om>pc_pi2) om = om-pc_pi2

     if(om/=om) then
!       write(0,*) 'omnan',d, xold, x, zold, z, omold, om, mu
!       stop 'transport2: om is nan'
        ierr = 6
        return
     endif!}}}
  endif !d>0d0

!-- tallying energy densities
  if(.not.trn_isimcanlog) then
!-- depositing nonanalog absorbed energy
     edep = e* &
          (1d0-exp(-grd_capgam(ic)* &
          elabfact*d*thelp))*elabfact
     if(edep/=edep) then
!       stop 'transport2_gamgrey: invalid energy deposition'
        ierr = 10
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgam(ic) * &
          elabfact*d*thelp)
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) then
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
  endif

!
!-- checking which event occurs from min distance

!-- common manipulations for collisions
  if(d==dcol) then
!-- resampling direction
     lredir = .true.
     call rnd_r(r1,rndstate)
     mu = 1d0 - 2d0*r1
     call rnd_r(r1,rndstate)
     om = pc_pi2*r1
!-- checking velocity dependence
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
!-- recalculating dirdotu
        dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     endif
  elseif(any([dbx,dby]==d)) then
!-- checking if escapted domain
     loutx = d==dbx.and.ixnext==grd_nx+1
     louty = d==dby.and.(iynext==grd_ny+1.or.iynext==0)
     if(loutx.or.louty) then
!-- ending particle
        ptcl2%done = .true.
        ptcl2%lflux = .true.
        return
     endif
  endif

!-- effective collision
  if(d==dcol) then
!-- checking if analog!{{{
     if(trn_isimcanlog) then
!-- effective absorption:
!-- ending particle
        ptcl2%done=.true.
!-- adding comoving energy to deposition energy
        edep = e*elabfact
     else
!-- effectively scattered:
!-- transforming to cmf, then to lab:
!-- energy weight
        e = e*elabfact/(1d0-dirdotu*cinv)
     endif!}}}

!
!-- x-bound
  elseif(d==dbx) then

     if(ixnext>ix) then
        x = grd_xarr(ix+1)
     else
        x = grd_xarr(ix)
     endif
!-- IMC in adjacent cell
     ix = ixnext
     ic = grd_icell(ix,iy,iz)
!
!-- y-bound
  elseif(d==dby) then

     if(iynext>iy) then
        y = grd_yarr(iy+1)
     else
        y = grd_yarr(iy)
     endif
!-- IMC in adjacent cell
     iy = iynext
     ic = grd_icell(ix,iy,iz)
  else
!    stop 'transport2_gamgrey: invalid distance'
     ierr = 12
     return
  endif


!-- update planar projections
  if(lredir) then
!-- planar projections (invariant until collision)
     mux = x*sin(om)/sin(z+om)  !-- intercept
     muy = x*sin(z)/sin(z+om)  !-- distance to intercept
     muz = pc_pi-(z+om)  !-- direction angle
     if(muz<0d0) muz = muz+pc_pi2
     if(muz>pc_pi2) muz = muz-pc_pi2
  endif

end subroutine transport2_gamgrey
