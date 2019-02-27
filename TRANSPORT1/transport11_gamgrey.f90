!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport11_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)

  use randommod
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
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time
!
  real*8 :: r1, thelp,thelpinv
  real*8 :: db, dcol
  real*8 :: darr(2)
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, muold
! real*8 :: x1, x2, xx0
  real*8 :: help
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix,ic
  integer,parameter :: iy=1,iz=1
  real*8,pointer :: x, mu, e, d

  ix => ptcl2%ix
  ic => ptcl2%ic
  d => ptcl2%dist
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0

  if(grd_isvelocity) then
     siglabfact = 1.0d0 - mu*x*cinv
     dcollabfact = tsp_t*(1d0-mu*x*cinv)
     thelp = tsp_t
  else
     siglabfact = 1d0
     dcollabfact = 1d0
     thelp = 1d0
  endif
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*dt*thelpinv) !> dcen

!
!== DISTANCE CALCULATIONS
!
!-- distance to boundary = db
  if (ix == 1) then
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  elseif (mu < -sqrt(1d0-(grd_xarr(ix)/x)**2)) then
     db = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
  else
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  endif
!-- sanity check
  if(db/=db) then
!   stop 'transport11_gamgrey: db/=db'
    ierr = 1
    return
  endif
!
!-- distance to fictitious collision = dcol
  if(grd_capgam(ic)<=0d0 .or. .not.trn_isimcanlog) then
     dcol = far
  else
     call rnd_r(r1,rndstate)
     dcol = abs(log(r1)/(grd_capgam(ic)*dcollabfact))
  endif
!
!-- minimum distance = d
!  if(tsp_it==29) write(*,*) dcol,dthm,db,dcen,ddop
  darr = [dcol,db]
  if(any(darr/=darr) .or. any(darr<0d0)) then
!    write(0,*) darr
!    write(*,*) ix,x,mu
!    stop 'transport11_gamgrey: invalid distance'
     ierr = 2
     return
  endif
  d = minval(darr)
!
!== END OF DISTANCE CALCULATIONS
!
!-- position, angle, time update  
  rold = x
  x = sqrt((1d0-mu**2)*x**2 + (d+x*mu)**2)
!  x = sqrt(x**2+d**2+2d0*d*x*mu)
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (rold*mu+d)/x
  endif

!-- transformation factor set
  if(grd_isvelocity) then
     elabfact = 1d0 - muold*rold*cinv
  else
     elabfact = 1d0
  endif
  !calculating energy deposition and density
  !
  if(.not.trn_isimcanlog) then
     edep = e*(1d0-exp( &
          -grd_capgam(ic)*siglabfact*d*thelp))*elabfact
     !--
     e = e*exp(-grd_capgam(ic)*siglabfact*d*thelp)

  endif

!-- transformation factor reset
  if(grd_isvelocity) then
     elabfact = 1d0 - mu*x*cinv
  else
     elabfact = 1d0
  endif

!
!-- fictitious scattering with implicit capture
  if (d == dcol) then
     !!{{{
     call rnd_r(r1,rndstate)
     if(r1<=1d0.and.trn_isimcanlog) then
        ptcl2%stat = 'dead'
        edep = e*elabfact
!-- velocity effects accounting
!
     else
        call rnd_r(r1,rndstate)
        mu = 1d0-2d0*r1
        if(abs(mu)<0.0000001d0) then
           mu = 0.0000001d0
        endif
        if(grd_isvelocity) then
           mu = (mu+x*cinv)/(1d0+x*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-mu*x*cinv)
!
           e = e*elabfact*help
           
        endif
!
        call rnd_r(r1,rndstate)
     endif
     !!}}}
!
!------boundary crossing ----
  elseif (d == db) then
     if (mu>=0d0) then!{{{
        if (ix == grd_nx) then
           ptcl2%stat = 'flux'
        else
           x = grd_xarr(ix+1)
           ix = ix+1
        endif
     else
        if (ix==1) then
           x = grd_xarr(ix+1)
           ix = ix+1
        else
           x = grd_xarr(ix)
           ix = ix-1
        endif
     endif
     ic = grd_icell(ix,iy,iz)!}}}
  endif

end subroutine transport11_gamgrey
! vim: fdm=marker
