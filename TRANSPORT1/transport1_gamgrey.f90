!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport1_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)

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
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: dt = pc_year !give grey transport infinite time
!
  logical :: lout
  real*8 :: elabfact, eta, xi
  real*8,pointer :: mux,muy,muz
  real*8 :: thelp, thelpinv, help
  real*8 :: dcol,dbx,dby,dbz
  real*8 :: darr(4)
  real*8 :: r1

  integer :: ixnext,iynext,iynext1,iynext2,iznext
  integer :: idby1,idby2
  real*8 :: yhelp1,yhelp2,yhelp3,yhelp4,dby1,dby2
  real*8 :: zhelp
  real*8 :: xold,yold,zold
  real*8 :: muold,omold
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix,iy,iz,ic
  real*8,pointer :: x,y,z,mu,om,e,d

  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e

  mux => ptcl2%mux
  muy => ptcl2%muy
  muz => ptcl2%muz
  d => ptcl2%dist
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0

!-- spherical projections
  xi = sqrt(1d0-mu**2)*sin(om)

  idby1=0
  idby2=0
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
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*dt*thelpinv) !> dcen

!-- radial boundary distance (x)
  if(ix==1 .or. mu>=-sqrt(1d0-(grd_xarr(ix)/x)**2)) then
     dbx = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
     ixnext = ix+1
  else
     dbx = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
     ixnext = ix-1
  endif
!-- sanity check
  if(dbx/=dbx) then
!    stop 'transport1_gamgrey: dbx/=dbx'
     ierr = 1
     return
  endif
!
!-- polar boundary distance (y)
!-- dby1: iy->iy-1
  if(iy/=1) then
     yhelp1=grd_yarr(iy)**2-muz**2
     yhelp2=mu*grd_yarr(iy)**2-muz*y
     yhelp3=grd_yarr(iy)**2-y**2
  endif
  if(y==grd_yarr(iy).or.(y==grd_yarr(iy+1).and. &
       abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9)) yhelp3=0d0
  if(iy==1) then
!-- don't stop at axis
     idby1 = 8
     dby1 = far
  elseif(yhelp1==0d0.and.yhelp3==0d0) then
!-- particle, direction on cone
     idby1=1
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
  if(iy/=grd_ny) then
     yhelp1=grd_yarr(iy+1)**2-muz**2
     yhelp2=mu*grd_yarr(iy+1)**2-muz*y
     yhelp3=grd_yarr(iy+1)**2-y**2
  endif
  if(y==grd_yarr(iy+1).or.(y==grd_yarr(iy).and. &
       abs(grd_yarr(iy)+grd_yarr(iy+1))<1d-9)) yhelp3=0d0
  if(iy==grd_ny) then
!-- don't stop at axis
     idby2 = 9
     dby2 = far
  elseif(yhelp1==0d0.and.yhelp3==0d0) then
!-- particle, direction on cone
     idby2=1
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
!    write(0,*) iy, y
!    write(0,*) idby1, dby1, idby2, dby2
!    stop 'transport1_gg: dby1<0 and dby2<0'
     ierr = 2
     return
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
  if(xi==0d0 .or. grd_nz==1) then
     dbz = far
     iznext=iz
  elseif(xi>0d0 .and. z>grd_zarr(iz+1)-pc_pi) then
!-- counterclockwise
     iznext=iz+1
     zhelp = muy*cos(grd_zarr(iz+1))-mux*sin(grd_zarr(iz+1))
!-- at boundary already
     if(z==grd_zarr(iz+1)) then
        dbz = 0d0
!-- parallel to plane
     elseif(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz+1)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  elseif(xi<0d0 .and. z<grd_zarr(iz)+pc_pi) then
!-- clockwise
     iznext=iz-1
     zhelp = muy*cos(grd_zarr(iz))-mux*sin(grd_zarr(iz))
!-- at boundary already
     if(z==grd_zarr(iz)) then
        dbz = 0d0
!-- parallel to plane
     elseif(zhelp==0d0) then
        dbz = far
     else
        dbz = x*sqrt(1d0-y**2)*sin(grd_zarr(iz)-z)/zhelp
        if(dbz<=0d0) dbz = far
     endif
  else
     dbz = far
     iznext=iz
  endif

!-- effective collision distance
  if(grd_capgam(ic)<=0d0 .or. .not.trn_isimcanlog) then
     dcol = far
  else
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ic))
  endif

!
!-- finding minimum distance
  darr = [dbx,dby,dbz,dcol]
  ptcl2%idist = minloc(darr,dim=1)
  d = darr(ptcl2%idist)
  if(any(darr/=darr)) then
     ierr = 3
     return
  endif
  if(d<0d0) then
     ierr = 4
     return
  endif


!
!-- update position
!-- storing old position
  xold = x
  yold = y
  zold = z
  muold = mu
  omold = om

!-- x position
  if(d==dbx) then
!-- on boundary
     if(ixnext>ix) then
        x = grd_xarr(ix+1)
     else
        x = grd_xarr(ix)
     endif
  else
!-- in cell
     ixnext = ix
     x = sqrt((1d0-mu**2)*x**2+(d+x*mu)**2)
     if(x/=x) then
   !    stop 'transport1: x/=x'
        ierr = 5
        return
     endif
     if(d==0d0) x = xold
  endif

!-- y position
  if(d==dby) then
!-- on boundary
     if(iynext==iy-1) then
        if(iynext<1) then
!          stop 'transport1: iynext<1'
           ierr = 11
           return
        endif
        y = grd_yarr(iy)
     elseif(iynext==iy+1) then
        if(iynext>grd_ny) then
!          stop 'transport1: iynext>ny'
           ierr = 12
           return
        endif
        y = grd_yarr(iy+1)
     else
!-- sanity check
!       write(0,*) dby
!       write(0,*) y,grd_yarr(iy),grd_yarr(iy+1),iy,iynext
!       stop 'transport1: invalid polar bound crossing'
        ierr = 13
        return
     endif
  else
!-- in cell
     iynext = iy
     y = (xold*yold + muz*d)/x
     y = max(y,-1d0)
     y = min(y,1d0)
     if(d==0d0) y = yold
  endif

!-- z position
  if(d==dbz) then
     if(iznext==iz-1) then
        z = grd_zarr(iz)
        if(iznext==0) then
           iznext = grd_nz
           z = pc_pi2
        endif
     elseif(iznext==iz+1) then
        z = grd_zarr(iz+1)
        if(iznext==grd_nz+1) then
           iznext = 1
           z = 0d0
        endif
     elseif(d/=dby) then
!       stop 'transport1: invalid iznext'
        ierr = 15
        return
     endif
  else
!-- in cell
     iznext = iz
     z = atan2(xold*sqrt(1d0-yold**2)*sin(z)+muy*d , &
          xold*sqrt(1d0-yold**2)*cos(z)+mux*d)
     if(abs(z)<1d-9.and.iz==1) then
        z = 0d0
     elseif(abs(z)<1d-9.and.iz==grd_nz) then
        z = pc_pi2
     elseif(z<0d0) then
        z = z+pc_pi2
     endif
     if(d==0d0) z = zold
  endif


!
!-- update direction
!-- radial projection of direction
  mu = (xold*mu+d)/x
  mu = max(mu,-1d0)
  mu = min(mu,1d0)
!-- azimuthal angle of direction (about radius)
  eta = y*(cos(z)*mux+sin(z)*muy)-sqrt(1d0-y**2)*muz
  xi = cos(z)*muy-sin(z)*mux
  om = atan2(xi,eta)
  if(om<0d0) om = om+pc_pi2
!-- warn about inaccurate result
  if(abs(mu**2+eta**2+xi**2-1d0)>1d-9) then
     ierr = -1 !warning
!       write(0,*) 'transport1: invalid mu,eta,xi',mu**2+eta**2+xi**2-1d0,mu,eta,xi
  endif
!
!-- normalize direction
  help = 1d0/sqrt(mu**2+eta**2+xi**2)
  mu = mu*help
  !eta = eta*help
  !xi = xi*help


!
!-- special case
  if(x<1d-15*grd_xarr(2).and.muold==-1d0) then
!-- sanity check
     if(d==dbx) then
!       stop 'transport1: x<1d-15*xarr(2),d==dbx,mu=-1'
        ierr = 6
        return
     endif
!-- excluding dbz, thomson and collision
     dbz = far
     dcol = far
!-- resetting direction
     mu = 1d0
     om = omold
     !eta = 0d0
     !xi = 0d0
!-- reflecting y
     y = -yold
     iynext = binsrch(y,grd_yarr,grd_ny+1,.false.)
!-- reflecting z
     z = zold+pc_pi !z is not updated with atan2 calculation
     if(z>pc_pi2) z = z-pc_pi2
     if(grd_nz>1) iznext=binsrch(z,grd_zarr,grd_nz+1,.false.)
  else
     if(x<1d-15*grd_xarr(2)) then
!       stop 'transport1: x=0 and muold/=-1'
        ierr = 7
        return
     endif
  endif

!-- depositing nonanalog absorbed energy
  if(.not.trn_isimcanlog) then
     edep = e* &
          (1d0-exp(-grd_capgam(ic)* &
          elabfact*d*thelp))*elabfact
     if(edep/=edep) then
!       stop 'transport1_gamgrey: invalid energy deposition'
        ierr = 8
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgam(ic)*elabfact*d*thelp)
  endif
!
!-- updating transformation factors
  if(grd_isvelocity) elabfact = 1d0-mu*x*cinv

!-- common manipulations for collisions
  if(d==dcol) then
!-- checking if escaped domain
  elseif(d==dbx) then
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- ending particle
        ptcl2%stat = 'flux'
!-- redefine for flux tally
        mu = muz
        om = atan2(muy,mux)
        if(om<0d0) om=om+pc_pi2
        return
     endif
  endif

  if(d==dcol) then
     ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
     edep = edep + e*elabfact
  endif

  ix = ixnext
  iy = iynext
  iz = iznext
  ic = grd_icell(ix,iy,iz)

end subroutine transport1_gamgrey
! vim: fdm=marker
