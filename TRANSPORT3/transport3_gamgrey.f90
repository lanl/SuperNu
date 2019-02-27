!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport3_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)

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

  logical :: loutx,louty,loutz
  integer :: iznext,iynext,ixnext
  real*8 :: elabfact, eta, xi
  real*8 :: thelp, thelpinv, help
  real*8 :: dcol,dbx,dby,dbz
  real*8 :: darr(4)
  real*8 :: r1
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix, iy, iz, ic
  real*8,pointer :: x,y,z,mu,om,e,d
!-- statement functions
  integer :: l,j,k
  real*8 :: xm,ym,zm,rsq
  xm(l) = .5*(grd_xarr(l+1) + grd_xarr(l))
  ym(l) = .5*(grd_yarr(l+1) + grd_yarr(l))
  zm(l) = .5*(grd_zarr(l+1) + grd_zarr(l))
  rsq(l,j,k) = xm(l)**2 + ym(j)**2 + zm(k)**2

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

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
!
!-- projections
  eta = sqrt(1d0-mu**2)*sin(om)
  xi = sqrt(1d0-mu**2)*cos(om)
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     elabfact=1d0-(mu*z+eta*y+xi*x)*cinv
     thelp = tsp_t
  else
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
  if(xi==0d0) then
     dbx = far
  else
     dbx = max((grd_xarr(ix)-x)/xi,(grd_xarr(ix+1)-x)/xi)
  endif
  if(eta==0d0) then
     dby = far
  else
     dby = max((grd_yarr(iy)-y)/eta,(grd_yarr(iy+1)-y)/eta)
  endif
  if(mu==0d0) then
     dbz = far
  else
     dbz = max((grd_zarr(iz)-z)/mu,(grd_zarr(iz+1)-z)/mu)
  endif
!
!-- effective collision distance
  if(grd_capgam(ic)<=0d0) then
     dcol = far
  elseif(trn_isimcanlog) then
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ic))
  else
     dcol = far
  endif
!
!-- finding minimum distance
  darr = [dbx,dby,dbz,dcol]
  ptcl2%idist = minloc(darr,dim=1)
  d = darr(ptcl2%idist)
  if(any(darr/=darr)) then
     ierr = 1
     return
  endif
  if(d<0d0) then
     ierr = 2
     return
  endif

!-- updating position
  if(d==dbx) then
     if(xi>=0d0) then
        ixnext = ix+1
        x = grd_xarr(ix+1)
     else
        ixnext = ix-1
        x = grd_xarr(ix)
     endif
  else
     ixnext = ix
     x = x + xi*d
  endif

  if(d==dby) then
     if(eta>=0d0) then
        iynext = iy+1
        y = grd_yarr(iy+1)
     else
        iynext = iy-1
        y = grd_yarr(iy)
     endif
  else
     iynext = iy
     y = y + eta*d
  endif

  if(d==dbz) then
     if(mu>=0d0) then
        iznext = iz + 1
        z = grd_zarr(iz+1)
     else
        iznext = iz - 1
        z = grd_zarr(iz)
     endif
  else
     iznext = iz
     z = z + mu*d
  endif

!-- tallying energy densities
  if(.not.trn_isimcanlog) then
!-- depositing nonanalog absorbed energy
     edep = e* &
          (1d0-exp(-grd_capgam(ic)* &
          elabfact*d*thelp))*elabfact
     if(edep/=edep) then
!       stop 'transport3_gamgrey: invalid energy deposition'
        ierr = 3
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgam(ic) * &
          elabfact*d*thelp)
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) then
     elabfact=1d0-(xi*x+eta*y+mu*z)*cinv
  endif

!
!-- checking which event occurs from min distance
!
!-- checking if escaped domain
  loutx = .false.
  if(d==dbx) then
     if(ixnext==grd_nx+1 .or. ixnext==0) then
        loutx = .true. !domain edge is reached
     elseif(grd_icell(ixnext,iy,iz)==grd_ivoid) then
        if(ixnext>ix.eqv.grd_xarr(ixnext)>0d0) then !away from center
           loutx = rsq(ixnext,iy,iz)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  louty = .false.
  if(d==dby) then
     if(iynext==grd_ny+1 .or. iynext==0) then
        louty = .true. !domain edge is reached
     elseif(grd_icell(ix,iynext,iz)==grd_ivoid) then
        if(iynext>iy.eqv.grd_yarr(iynext)>0d0) then !away from center
           louty = rsq(ix,iynext,iz)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  loutz = .false.
  if(d==dbz) then
     if(iznext==grd_nz+1 .or. iznext==0) then
        loutz = .true. !domain edge is reached
     elseif(grd_icell(ix,iy,iznext)==grd_ivoid) then
        if(iznext>iz.eqv.grd_zarr(iznext)>0d0) then !away from center
           loutz = rsq(ix,iy,iznext)>grd_rvoid**2 !enter void corners
        endif
     endif
  endif
  if(loutx.or.louty.or.loutz) then
!-- ending particle
     ptcl2%stat = 'flux'
     return
  endif

!-- common manipulations for collisions
  if(d==dcol) then
!-- resampling direction
     call rnd_r(r1,rndstate)
     mu = 1d0 - 2d0*r1
     call rnd_r(r1,rndstate)
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) then
        eta = sqrt(1d0-mu**2)*sin(om)
        xi = sqrt(1d0-mu**2)*cos(om)
!-- transforming mu
        mu = (mu+z*cinv)/(1d0+(mu*z+eta*y+xi*x)*cinv)
        if(mu>1d0) then
           mu = 1d0
        elseif(mu<-1d0) then
           mu = -1d0
        endif
!-- transforming om
        om = atan2(eta+y*cinv,xi+x*cinv)
        if(om<0d0) om=om+pc_pi2
!-- x,y lab direction cosines
        eta = sqrt(1d0-mu**2)*sin(om)
        xi = sqrt(1d0-mu**2)*cos(om)
     endif
  endif

!-- effective collision
  if(d==dcol) then
     call rnd_r(r1,rndstate)
!-- checking if analog
     if(trn_isimcanlog) then
!-- effective absorption
        ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
        edep = e*elabfact
        return
     else
!-- transforming to lab
        if(grd_isvelocity) then
           help = elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
!
!-- energy weight
           e = e*help
        endif
     endif
  endif

  ix = ixnext
  iy = iynext
  iz = iznext
  ic = grd_icell(ix,iy,iz)

end subroutine transport3_gamgrey
! vim: fdm=marker
