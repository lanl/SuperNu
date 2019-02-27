!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
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
  integer :: ixnext,iynext,iznext
  real*8 :: elabfact, dirdotu, gm, xi
  real*8,pointer :: mux,muy,muz
  real*8 :: thelp, thelpinv, help
  real*8 :: dcol,dbx,dby,dbz
  real*8 :: darr(4)
  real*8 :: xold, omold, zold
  real*8 :: r1
!-- distance out of physical reach
  real*8 :: far

  real*8 :: mux1,mux2
  real*8 :: angrat1,angrat2,angrat3

  integer,pointer :: ix, iy, iz, ic
  real*8,pointer :: x,y,z,mu,om,e,d
!-- statement functions
  integer :: l,j
  real*8 :: muxf,zz,xm,ym,rsq
  muxf(zz) = ptcl%x*sin(ptcl2%muz+zz)/sin(ptcl2%muz)
  !muyf(zz) = ptcl%x*sin(zz)/sin(ptcl2%muz)
  xm(l) = .5*(grd_xarr(l+1) + grd_xarr(l))
  ym(l) = .5*(grd_yarr(l+1) + grd_yarr(l))
  rsq(l,j) = xm(l)**2 + ym(j)**2

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
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

!-- azimuthal projection
  xi = sqrt(1d0-mu**2)*sin(om)

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

!-- azimuthal boundary distance
  if(xi==0d0 .or. grd_nz==1) then
     dbz = far
  elseif(xi>0d0 .and. z>grd_zarr(iz+1)-pc_pi) then
!-- counterclockwise
     iznext=iz+1
     help = sqrt(1d0-mu**2)*sin(om+z-grd_zarr(iz+1))
     if(z==grd_zarr(iz+1)) then
        dbz = 0d0
     elseif(help==0d0) then
        dbz = far
     else
        dbz = x*sin(grd_zarr(iz+1)-z)/help
        if(dbz<=0d0) dbz = far
     endif
  elseif(xi<0d0 .and. z<grd_zarr(iz)+pc_pi) then
!-- clockwise
     iznext=iz-1
     help = sqrt(1d0-mu**2)*sin(om+z-grd_zarr(iz))
     if(z==grd_zarr(iz)) then
        dbz = 0d0
     elseif(help==0d0) then
        dbz = far
     else
        dbz = x*sin(grd_zarr(iz)-z)/help
        if(dbz<=0d0) dbz = far
     endif
  else
     dbz = far
  endif
  
!
!-- effective collision distance
  if(grd_capgam(ic)<=0d0 .or. .not.trn_isimcanlog) then
!-- making greater than dcen
     dcol = far
  else
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ic))
  endif
!
!-- finding minimum distance
  darr = [dcol,dbx,dby,dbz]
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
!-- update position
!
!-- store prev position
  xold = x
  zold = z
  omold = om
!
!-- update distance to intercept
  muy = muy + d*sqrt(1d0-mu**2)

!-- update y position
  if(d==dby) then
!-- on boundary
     if(iynext>iy) then
        y = grd_yarr(iy+1)
     else
        y = grd_yarr(iy)
     endif
  else
!-- in cell
     y = y + mu*d
  endif

!-- update x position
  if(d==dbx) then
!-- on boundary
     if(ixnext>ix) then
        x = grd_xarr(ix+1)
     else
        x = grd_xarr(ix)
     endif
  else
!-- in cell
     if(abs(abs(cos(muz))-1d0)<1d-2) then
!-- muz calculation unreliable
        x = sqrt(xold**2 + (1d0-mu**2)*d**2 + &
           2d0*xold*sqrt(1d0-mu**2)*d*cos(omold))
     else
        x = sqrt(mux**2 + muy**2 - 2d0*mux*muy*cos(muz))
     endif
!
!-- special case: particle close to corners
     if(abs(x-xold)<1d-15*(x+xold)) then
        if(xold==grd_xarr(ix)) x = xold
        if(xold==grd_xarr(ix+1)) x = xold
     endif

     if(d==0) x = xold
  endif

!-- update z position
  if(d==dbz) then
!-- on boundary
     if(iznext==iz+1) then
        z = grd_zarr(iznext)
        if(iznext==grd_nz+1) then
           iznext = 1
           z = 0d0
        endif
     elseif(iznext==iz-1) then
        z = grd_zarr(iz)
        if(iznext==0) then
           iznext = grd_nz
           z = pc_pi2
        endif
     else
        ierr = 15
        return
     endif
  else
!-- in cell
     if(x==0) then
!-- todo: implement exception for x==0
        ierr = 8
        return
     endif
!-- trigoniometric ratios
     angrat1 = (x**2 + mux**2 - muy**2)/(2*x*mux)
     angrat2 = sin(muz)*muy/x

     if(abs(angrat1)<1d0 .and. abs(angrat1)>1d-15 .and. &
           abs(angrat1)<abs(angrat2)) then
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
     if(d==0d0) z = zold
  endif


!
!-- update om
  om = omold - (z - zold)
  if(om<0d0) om = om+pc_pi2
  if(om>pc_pi2) om = om-pc_pi2

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

!-- checking if escaped domain
  loutx = .false.
  if(d==dbx) then
     if(ixnext==grd_nx+1) then
        loutx = .true. !domain edge is reached
     elseif(ixnext>ix .and. grd_icell(ixnext,iy,iz)==grd_ivoid) then
        loutx = rsq(ixnext,iy)>min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1))**2 !enter void corners
     endif
  endif
  louty = .false.
  if(d==dby) then
     if(iynext==grd_ny+1 .or. iynext==0) then
        louty = .true. !domain edge is reached
     elseif(grd_icell(ix,iynext,iz)==grd_ivoid) then
        if(iynext>iy.and.grd_yarr(iynext)>0d0 .or. iynext<iy.and.grd_yarr(iynext)<0d0) then
           louty = rsq(ix,iynext)>min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1))**2 !enter void corners
        endif
     endif
  endif
  if(loutx.or.louty) then
!-- ending particle
     ptcl2%stat = 'flux'
!-- redefine for flux tally
     om = muz
     return
  endif

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
  endif

!-- effective collision
  if(d==dcol) then
     ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
     edep = e*elabfact
!
!-- x-bound
  elseif(d==dbx) then
!-- IMC in adjacent cell
     ix = ixnext
     ic = grd_icell(ix,iy,iz)
!
!-- y-bound
  elseif(d==dby) then
!-- IMC in adjacent cell
     iy = iynext
     ic = grd_icell(ix,iy,iz)
!
!-- z-bound
  elseif(d==dbz) then
!-- IMC in adjacent cell
     iz = iznext
     ic = grd_icell(ix,iy,iz)
!
  else
!    stop 'transport2_gamgrey: invalid distance'
     ierr = 13
     return
  endif


  if(om/=om) then
!       write(0,*) 'omnan',d, xold, x, zold, z, omold, om, mu
!       stop 'transport2: om is nan'
     ierr = 6
     return
  endif


!-- update planar projections
  if(lredir) then
!-- planar projections (invariant until collision)
     mux = x*sin(om)/sin(z+om)  !-- intercept
     muy = x*sin(z)/sin(z+om)  !-- distance to intercept
     muz = pc_pi-(z+om)  !-- direction angle
     if(muz<0d0) muz = muz+pc_pi2
     if(muz<0d0) muz = muz+pc_pi2
  endif

end subroutine transport2_gamgrey
! vim: fdm=marker
