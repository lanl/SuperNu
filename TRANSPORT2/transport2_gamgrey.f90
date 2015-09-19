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

  logical :: loutx,louty
  integer :: ihelp,iznext
  real*8 :: elabfact, dirdotu, gm, xi
  real*8 :: thelp, thelpinv, zhelp
  real*8 :: dcol,dbx,dby,dbz
  real*8 :: darr(4)
  real*8 :: xold, omold
  real*8 :: r1
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x,y,z,mu,om,e,d

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0

!-- azimuthal projection
  xi = sqrt(1d0-mu**2)*sin(om)
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
!-- calculating distance to boundary:
!-- to x-bound
  if(abs(mu)==1d0) then
!-- making greater than dcen
     dbx = far
  else
     if(abs(sin(om))<grd_xarr(ix)/x .and. cos(om)<0d0) then
!-- inner boundary
        dbx = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
        if(dbx/=dbx) then
!          stop 'transport2_gamgrey: invalid inner dbx'
           ierr = 1
           return
        endif
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbx = 0d0
     else
!-- outer boundary
        dbx = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
        if(dbx/=dbx) then
!          stop 'transport2_gamgrey: invalid outer dbx'
           ierr = 2
           return
        endif
     endif
  endif
  if(dbx/=dbx) then
!    stop 'transport2_gamgrey: dbx nan'
     ierr = 3
     return
  endif

!-- to y-bound
  if(mu>0d0) then
     dby = (grd_yarr(iy+1)-y)/mu
     if(dby<0d0) then
!       stop 'upward dby'
        ierr = 4
        return
     endif
     if((grd_yarr(iy)-y)/mu>0d0) then
!       stop 'transport2_gamgrey: y below cell'
        ierr = 5
        return
     endif
  elseif(mu<0d0) then
     dby = (grd_yarr(iy)-y)/mu
     if(dby<0d0) then
!       stop 'downward dby'
        ierr = 6
        return
     endif
     if((grd_yarr(iy+1)-y)/mu>0d0) then
!       stop 'transport2_gamgrey: y above cell'
        ierr = 7
        return
     endif
  else
!-- making greater than dcen
     dby = far
  endif

!-- azimuthal boundary distance
  if(xi==0d0 .or. grd_nz==1) then
     dbz = far
  elseif(xi>0d0 .and. z>grd_zarr(iz+1)-pc_pi) then
!-- counterclockwise
     iznext=iz+1
     zhelp = sqrt(1d0-mu**2)*sin(om+z-grd_zarr(iz+1))
     if(z==grd_zarr(iz+1)) then
        dbz = 0d0
     elseif(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sin(grd_zarr(iz+1)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  elseif(xi<0d0 .and. z<grd_zarr(iz)+pc_pi) then
!-- clockwise
     iznext=iz-1
     zhelp = sqrt(1d0-mu**2)*sin(om+z-grd_zarr(iz))
     if(z==grd_zarr(iz)) then
        dbz = 0d0
     elseif(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sin(grd_zarr(iz)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  else
     dbz = far
  endif

!
!-- calculating distance to effective collision:
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
  darr = [dbx,dby,dbz,dcol]
  if(any(darr/=darr) .or. any(darr<0d0)) then
!    write(0,*) darr
!    write(*,*) ix,iy,x,y,mu,om
!    stop 'transport2_gamgrey: invalid distance'
     ierr = 8
     return
  endif
  d = minval(darr)

!-- updating position
  y = y+mu*d
  xold = x
  x = sqrt(xold**2+(1d0-mu**2)*d**2+2d0*xold*sqrt(1d0-mu**2)*d*cos(om))
!-- updating azimuthal direction
  omold = om
  if(om>pc_pi.and.om<pc_pi2) then
     om = pc_pi2 + &
          atan2(-sqrt(max(x**2-(xold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          xold*cos(omold)+d*sqrt(1d0-mu**2))
  else
     om = atan2(sqrt(max(x**2-(xold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          xold*cos(omold)+d*sqrt(1d0-mu**2))
  endif
  if(om/=om) then
!    stop 'transport2_gamgrey: om is nan'
     ierr = 9
     return
  endif
!-- updating azimuthal position
  z = z+(om-omold)
  if(z<0d0) z=z+pc_pi2
  if(z>pc_pi2) z=z-pc_pi2

!
!-- updating transformation factors
  if(grd_isvelocity) then
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
  endif

!
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
!-- checking which event occurs from min distance

!-- common manipulations for collisions
  if(d==dcol) then
!-- resampling direction
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
     loutx = d==dbx.and.(sqrt(1d0-mu**2)*cos(om)>=0d0.and.ix==grd_nx)
     louty = d==dby.and.((mu>=0d0.and.iy==grd_ny).or.(mu<0.and.iy==1))
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

     if(cos(om)>=0d0) then
        ihelp = 1
        x = grd_xarr(ix+1)
     else
        if(ix==1) then
!          stop 'transport2_gamgrey: cos(om)<0 and ix=1'
           ierr = 11
           return
        endif
        ihelp = -1
        x = grd_xarr(ix)
     endif
!-- IMC in adjacent cell
     ix = ix+ihelp
     ic = grd_icell(ix,iy,iz)
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
     ic = grd_icell(ix,iy,iz)
!
!-- z-bound
  elseif(d==dbz) then
!-- sanity check
     if(grd_nz==1) then
!       stop 'transport2_gamgrey: invalid z crossing'
        ierr = 12
        return
     endif
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
!       stop 'transport1_gamgrey: invalid iznext'
        ierr = 13
        return
     endif
     iz = iznext
     ic = grd_icell(ix,iy,iz)
  else
!    stop 'transport2_gamgrey: invalid distance'
     ierr = 14
     return
  endif

end subroutine transport2_gamgrey
