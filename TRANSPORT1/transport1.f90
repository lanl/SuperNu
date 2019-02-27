!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport1(ptcl,ptcl2,rndstate,edep,eraddens,eamp,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens, eamp
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event.  If
!the puretran boolean is set to false, this routine couples to the
!corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c

  logical :: lredir !direction resampled
  logical :: lout
  real*8 :: elabfact, eta, xi
  real*8,pointer :: mux,muy,muz
  real*8 :: thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop
  real*8 :: darr(7),darrold(7)
  real*8 :: r1,r2

  integer :: ixnext,iynext,iynext1,iynext2,iznext
! integer :: ixold,iyold,izold,idistold
  integer :: idby1,idby2
  real*8 :: yhelp1,yhelp2,yhelp3,yhelp4,dby1,dby2
  real*8 :: zhelp
  real*8 :: xold,yold,zold
  real*8 :: muold,omold
! real*8 :: etaold,xiold
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump

  integer,pointer :: ix,iy,iz,ic,ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl,d
!-- statement functions
  integer :: l
  real*8 :: dx,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
! dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5d0*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25d0*(grd_yarr(l+1)+grd_yarr(l))**2)

  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

  mux => ptcl2%mux
  muy => ptcl2%muy
  muz => ptcl2%muz
  d => ptcl2%dist
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
  eamp = 0d0

!-- direction resample flag
  lredir = .false.

!
!-- spherical projections
  eta = sqrt(1d0-mu**2)*cos(om)
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
  far = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen

!-- census distance
  dcen = abs(pc_c*(tsp_t1-ptcl%t)*thelpinv)
!
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
!    stop 'transport1: dbx/=dbx'
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
!  write(0,*) idby1,dby1,idby2,dby2
!  if(y<grd_yarr(iy+1).and.y>grd_yarr(iy)) write(0,*) idby1,idby2
!  if(dby1==0d0.and.idby1==4) write(0,*) '1: ',idby1, y, iy, dby1
!  if(dby2==0d0.and.idby2==4) write(0,*) '2: ',idby2, y, iy
  ! if(dby1==0d0.and.dby2==0d0) stop 'transport1: invalid dby[1,2]'
  if(dby1<0d0.and.dby2<0d0) then
!    write(0,*) iy, y
!    write(0,*) idby1, dby1, idby2, dby2
!    stop 'transport1: dby1<0 and dby2<0'
     ierr = 2
     return
  endif
  ! if(dby1==0d0) then
  !    write(0,*) '1: ', iy, y, dby2, eta
  ! endif
  ! if(dby2==0d0) then
  !    write(0,*) '2: ', iy, y, dby1, eta
  ! endif
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
     iznext = iz
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
     iznext = iz
  endif

!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
     dthm = far
  endif

!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = far
  elseif(trn_isimcanlog) then
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ic))
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
     dcol = far
  endif

!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = pc_c*(elabfact-wl*grp_wlinv(ig+1))
     if(ddop<0d0) then
        ddop = far
     endif
  else
     ddop = far
  endif
!
!-- finding minimum distance
  darr = [dcen,dcol,dthm,ddop,dbx,dby,dbz]
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
  if(d/=dthm .and. d/=dcol) then
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
  endif


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
     dthm = far
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

!
!-- updating time
  ptcl%t = ptcl%t + thelp*d*cinv

!-- tallying energy densities
  if(trn_isimcanlog) then
!-- analog energy density
     eraddens = e*elabfact**2 * &
          d*thelp*cinv*tsp_dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)* &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
          thelp>1d-6) then
        eraddens = e* &
             (1.0d0-exp(-grd_fcoef(ic)*elabfact* &
             grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*elabfact * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        eraddens = e*elabfact**2 * &
             d*thelp*cinv*tsp_dtinv
     endif
!-- depositing nonanalog absorbed energy
     edep = e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic)* &
          elabfact*d*thelp))*elabfact
     if(edep/=edep) then
!       stop 'transport1: invalid energy deposition'
        ierr = 8
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
          elabfact*d*thelp)
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) elabfact = 1d0-mu*x*cinv

!
!-- census
  if(d==dcen) then
     ptcl2%stat = 'cens'
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     lredir = .true.
     call rnd_r(r1,rndstate)
     mu = 1d0 - 2d0*r1
     call rnd_r(r1,rndstate)
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+mu*x*cinv)
!-- checking if escaped domain
  elseif(d==dbx) then
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- observer time correction
        ptcl%t=ptcl%t-mu*x*thelp*cinv
!-- ending particle
        ptcl2%stat = 'flux'
!-- redefine for flux tally
        mu = muz
        om = atan2(muy,mux)
        if(om<0d0) om = om+pc_pi2
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
        totevelo=totevelo+e*(1d0-help)
!-- energy weight
        e = e*help
        e0 = e0*help
     endif

!
!-- effective collision
  elseif(d==dcol) then
     call rnd_r(r1,rndstate)
!-- checking if analog
     if(trn_isimcanlog.and.r1<=grd_fcoef(ic)) then
!-- effective absorption
        ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
        edep = e*elabfact
!-- velocity effects accounting
        totevelo = totevelo+e*(1d0-elabfact)
        return
     else
!-- effective scattering
!-- transforming to lab
        if(grd_isvelocity) then
           help = elabfact/(1d0-mu*x*cinv)
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-help)
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!-- redistributing wavelength
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
        if(grp_ng==1) then
        elseif(emitlump<.99d0 .or. trn_nolumpshortcut .or. in_puretran) then
           call rnd_r(r1,rndstate)
           ig = emitgroup(r1,ic)
           if((grd_cap(ig,ic)+grd_sig(ic)) * &
              min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
              thelp>=trn_tauddmc .and. .not.in_puretran) ptcl2%itype = 2
!-- don't sample, it will end up in the lump anyway
        else
           ptcl2%itype = 2
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(9,ic))
        endif

!-- checking for DDMC in new group
        if(ptcl2%itype==2) then
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*mu*x*cinv
!-- energy weight
              e = e*(1d0-mu*x*cinv)
              e0 = e0*(1d0-mu*x*cinv)
           endif
           wl = 0d0 !workaround ifort 13.1.3 bug
        else
!-- uniformly in new group
           call rnd_r(r1,rndstate)
           wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
           if(grd_isvelocity) wl = wl*(1d0-mu*x*cinv)
        endif
     endif

!
!-- radial bound
  elseif(d==dbx .and. dbx<dby) then

     l = grd_icell(ixnext,iy,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ixnext),xm(ixnext)*dyac(iy),xm(ixnext) * &
          ym(iy)*dz(iz))*thelp<trn_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        ix = ixnext
        iy = iynext
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming x-cosine to cmf
           mu = (mu-x*cinv)/elabfact
!-- amplification factor
           if(.not.trn_noampfact .and. mu<0d0) then
              help = 1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              totevelo = totevelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*x*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
              e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ixnext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if(r1<help*(1d0+1.5d0*abs(mu))) then
           if(grd_isvelocity) then
              ptcl2%itype = 2
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ixnext
           iy = iynext
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
!-- resampling x-cosine
           lredir = .true.
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu=(ix-ixnext)*max(r1,r2)
!-- resampling azimuthal
           call rnd_r(r1,rndstate)
           om = pc_pi2*r1
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
        endif
     endif

!
!-- polar bound
  elseif(d==dby) then

     l = grd_icell(ix,iynext,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iynext),xm(ix)*ym(iynext) * &
          dz(iznext))*thelp<trn_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        ix = ixnext
        iy = iynext
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming y-cosine to cmf
           mu=(mu-x*cinv)/(1d0-x*mu*cinv)
        endif
        eta = sqrt(1d0-mu**2)*cos(om)
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             xm(ix)*dyac(iynext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if(r1 < help*(1d0+1.5d0*abs(eta))) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ixnext
           iy = iynext
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           if(iynext==iy) then
              if(y>grd_yarr(iy)) then
                 eta = -max(r1,r2)
              elseif(y<grd_yarr(iy+1)) then
                 eta = max(r1,r2)
              else
                 ierr = -4 !warning
!                write(0,*) 'transport1: no idea how to calculate new eta'
                 eta = 2d0*r1 - 1d0
              endif
           else
              eta = (iynext-iy)*max(r1,r2)
           endif

           lredir = .true.
           call rnd_r(r1,rndstate)
           xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
!-- resampling x-cosine
           mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(xi,eta)
           if(om<0d0) om = om+pc_pi2
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
        endif
     endif

!
!-- azimuthal bound
  elseif(d==dbz) then
!-- sanity check
     if(grd_nz==1) then
!       stop 'transport1: invalid z crossing'
        ierr = 14
        return
     endif

     l = grd_icell(ix,iy,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy) * &
          dz(iznext))*thelp<trn_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        ix = ixnext
        iy = iynext
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming z-cosine to cmf
           mu = (mu-x*cinv)/elabfact
           xi = sqrt(1d0-mu**2)*sin(om)
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             xm(ix)*ym(iy)*dz(iznext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5d0*abs(xi))) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ixnext
           iy = iynext
           iz = iznext
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
!-- resampling z-cosine
           if(grd_nz==2) then
!-- rejection for 2 cells can occur at zarr(iz) and zarr(iz+1)
              if(iznext==1) then
                 if(z==grd_zarr(2)) then
                    xi = max(r1,r2)
                 elseif(z==0d0) then
                    xi = -max(r1,r2)
                 endif
              elseif(iznext==2) then
                 if(z==grd_zarr(2)) then
                    xi = -max(r1,r2)
                 elseif(z==pc_pi2) then
                    xi = max(r1,r2)
                 endif
              endif
           elseif(iznext==iz+1.or.(iznext==1.and.iz==grd_nz)) then
              xi = -max(r1,r2)
           else
              xi = max(r1,r2)
           endif

           lredir = .true.
           call rnd_r(r1,rndstate)
           eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
!-- resampling x-cosine
           mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(xi,eta)
           if(om<0d0) om = om+pc_pi2
!-- transforming mu to lab
           if(grd_isvelocity) mu=(mu+x*cinv)/(1d0+x*mu*cinv)
!-- reverting z
           if(iznext==grd_nz.and.iz==1) then
              if(z==pc_pi2) z = 0d0
           elseif(iznext==1.and.iz==grd_nz) then
              if(z==0d0) z = pc_pi2
           endif
        endif
     endif

!
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) then
!       stop 'transport1: ddop and no velocity'
        ierr = 16
        return
     endif
     if(ig<grp_ng) then
!-- shifting group
        ig = ig+1
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if((grd_cap(ig,ic)+grd_sig(ic)) * &
          min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
          thelp>=trn_tauddmc .and. .not.in_puretran) then
        ptcl2%itype = 2
        if(grd_isvelocity) then
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  else
!    stop 'transport1: invalid distance'
     ierr = 17
     return
  endif


!-- update planar projections
  if(lredir) then
!-- spherical projections
     eta = sqrt(1d0-mu**2)*cos(om)
     xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
     mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
     muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
     muz = mu*y-eta*sqrt(1d0-y**2)
  endif

!-- sanity check
  if(abs(y)==1d0) then
     ierr = -7 !warning
  endif

!!-- verbose
!  if(ierr/=0) then
!     write(0,*) '(x): ',ix,xold,x,grd_xarr(ix),grd_xarr(ix+1)
!     write(0,*) '(y): ',iy,yold,y,grd_yarr(iy),grd_yarr(iy+1)
!     write(0,*) '(z): ',iz,zold,z,grd_zarr(iz),grd_zarr(iz+1)
!     write(0,*) 'y==1 ',abs(yold)-1d0, abs(y)-1d0
!     write(0,*) 'lredir',lredir
!     write(0,*) 'old ix,iy,iz:',ixold,iyold,izold
!     write(0,*) 'iynxt1|2:',iynext1,iynext2
!     write(0,*) 'dby1|2  :',dby1,dby2
!     write(0,*) 'idby1|2 : ',idby1,idby2
!     write(0,*) 'ptcl id:',ptcl2%ipart,ptcl2%istep
!     write(0,*) 'dirold: ', muold,etaold,xiold
!     write(0,*) 'dir: ', mux,muy,muz,mu,eta,xi
!     write(0,*) 'darrold:', idistold, darrold
!     write(0,*) 'darr   :', ptcl2%idist, darr
!     write(0,*) 'atan2:',sqrt(1d0-yold**2),sin(zold),cos(zold),xold*sqrt(1d0-yold**2),muy*d,mux*d
!!    z = atan2(xold*sqrt(1d0-yold**2)*sin(z)+muy*d , &
!!         xold*sqrt(1d0-yold**2)*cos(z)+mux*d)
!     write(0,*)
!  endif

!-- save
  darrold = darr

end subroutine transport1
! vim: fdm=marker
