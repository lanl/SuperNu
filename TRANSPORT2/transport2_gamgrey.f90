subroutine transport2_gamgrey(ptcl)

  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time
  integer, external :: binsrch

  logical :: loutx,louty
  integer :: imu, ihelp
  real*8 :: elabfact, dirdotu
  real*8 :: dtinv, thelp, thelpinv 
  real*8 :: dcol,dbx,dby,d
  real*8 :: rold, zold, omold
  real*8 :: r1

  integer,pointer :: ix,iy
  real*8,pointer :: x,y,mu,om,e,e0
!-- statement functions
  integer :: l
  real*8 :: dx,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

  ix => ptcl%ix
  iy => ptcl%iy
  x => ptcl%x
  y => ptcl%y
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
!
!-- shortcut
  dtinv = 1d0/dt
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
!
!-- calculating distance to boundary:
!-- to x-bound
  if(abs(mu)==1d0) then
!-- making greater than dcen
     dbx = 2d0*pc_c*dt*thelpinv
  else
     if(abs(sin(om))<grd_xarr(ix)/x .and. cos(om)<0d0) then
!-- inner boundary
        dbx = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
        if(dbx/=dbx) then
           stop 'transport2_gamgrey: invalid inner dbx'
        endif
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbx = 0d0
     else
!-- outer boundary
        dbx = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
        if(dbx/=dbx) then
           stop 'transport2_gamgrey: invalid outer dbx'
        endif
     endif
  endif
  if(dbx/=dbx) stop 'transport2_gamgrey: dbx nan'

!-- to y-bound
  if(mu>0d0) then
     dby = (grd_yarr(iy+1)-y)/mu
     if(dby<0d0) then
        stop 'upward dby'
     if((grd_yarr(iy)-y)/mu>0d0) stop &
          'transport2_gamgrey: y below cell'
     endif
  elseif(mu<0d0) then
     dby = (grd_yarr(iy)-y)/mu
     if(dby<0d0) stop 'downward dby'
     if((grd_yarr(iy+1)-y)/mu>0d0) stop &
          'transport2_gamgrey: y above cell'
  else
!-- making greater than dcen
     dby = 2d0*pc_c*dt*thelpinv
  endif

!
!-- calculating distance to effective collision:
  if(grd_capgam(ix,iy,1)<=0d0) then
!-- making greater than dcen
     dcol = 2d0*pc_c*dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ix,iy,1))
  else
!-- making greater than dcen
     dcol = 2d0*pc_c*dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dbx,dby,dcol)
  if(d<0d0) stop 'transport2_gamgrey: negative distance'

!-- updating position
  rold = x
  zold = y
  x = sqrt(rold**2+(1d0-mu**2)*d**2+2d0*rold*sqrt(1d0-mu**2)*d*cos(om))
  y = zold+mu*d
!-- updating azimuthal direction
  omold = om
  if(om>pc_pi.and.om<pc_pi2) then
     om = pc_pi2 + &
          atan2(-sqrt(max(x**2-(rold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-mu**2))
  else
     om = atan2(sqrt(max(x**2-(rold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-mu**2))
  endif
  if(om/=om) then
     stop 'transport2_gamgrey: om is nan'
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) then
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
  endif

!
!-- tallying energy densities
  if(.not.prt_isimcanlog) then
!-- depositing nonanalog absorbed energy
     grd_edep(ix,iy,1) = grd_edep(ix,iy,1)+e* &
          (1d0-exp(-grd_capgam(ix,iy,1)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ix,iy,1)/=grd_edep(ix,iy,1)) then
        stop 'transport2_gamgrey: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgam(ix,iy,1) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!-- common manipulations for collisions
  if(d==dcol) then
!-- resampling direction
     r1 = rand()
     mu = 1d0 - 2d0*r1
     r1 = rand()
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- calculating transformation factors
        dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
        om = atan2(sqrt(1d0-mu**2)*sin(om), &
             sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming to lab:
!-- y-cosine
        mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
        if(mu>1d0) then
           mu = 1d0
        elseif(mu<-1d0) then
           mu = -1d0
        endif
!-- azimuthal direction angle
        if(om<0d0) om=om+pc_pi2
!-- recalculating dirdotu
        dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     endif
  elseif(any([dbx,dby]==d)) then
!-- checking if escapted domain
     loutx = d==dbx.and.(sqrt(1d0-mu**2)*cos(om)>=0d0.and.ix==grd_nx)
     louty = d==dby.and.((mu>=0d0.and.iy==grd_ny).or.(mu<0.and.iy==1))
     if(loutx.or.louty) then
!-- ending particle
        prt_done = .true.
!-- retrieving lab frame flux group, polar bin
        imu = binsrch(mu,flx_mu,flx_nmu+1,0)
!-- tallying outbound luminosity
        flx_gamluminos(imu,1) = flx_gamluminos(imu,1)+e*dtinv
        flx_gamlumdev(imu,1) = flx_gamlumdev(imu,1)+(e*dtinv)**2
        flx_gamlumnum(imu,1) = flx_gamlumnum(imu,1)+1
        return
     endif
  endif

!-- effective collision
  if(d==dcol) then
!-- checking if analog!{{{
     if(prt_isimcanlog) then
!-- effective absorption:
!-- ending particle
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,1) = grd_edep(ix,iy,1) + e*elabfact
     else
!-- effectively scattered:
!-- transforming to cmf, then to lab:
!-- energy weight
        e = e*elabfact/(1d0-dirdotu*cinv)
        e0 = e0*elabfact/(1d0-dirdotu*cinv)
     endif!}}}

!
!-- x-bound
  elseif(d==dbx) then

     if(cos(om)>=0d0) then
        ihelp = 1
        x = grd_xarr(ix+1)
     else
        if(ix==1) stop &
             'transport2_gamgrey: cos(om)<0 and ix=1'
        ihelp = -1
        x = grd_xarr(ix)
     endif
!-- IMC in adjacent cell
     ix = ix+ihelp
!
!-- y-bound
  elseif(d==dby) then

     if(mu>=0d0) then
        ihelp = 1
        y = grd_yarr(iy+1)
     else
        ihelp = -1
        y = grd_yarr(iy)
     endif
!-- IMC in adjacent cell
     iy = iy+ihelp
  else
     stop 'transport2_gamgrey: invalid distance'
  endif

end subroutine transport2_gamgrey
