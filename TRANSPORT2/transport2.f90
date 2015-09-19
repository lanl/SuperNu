pure subroutine transport2(ptcl,ptcl2,rndstate,edep,eraddens,eamp,totevelo,ierr)

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
  logical :: loutx,louty
  integer :: ixnext,iynext,iznext
  real*8 :: elabfact, dirdotu, mu0, gm, xi
  real*8,pointer :: mux,muy,muz
  real*8 :: thelp, thelpinv, help, zhelp
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop
  real*8 :: darr(7)
  real*8 :: xold, omold, zold
  real*8 :: r1, r2
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump

  real*8 :: mux1,mux2
  real*8 :: angrat1,angrat2,angrat3

  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x,y,z,mu,om,e,e0,wl,d
!-- statement functions
  integer :: l
  real*8 :: dx,dy,xm,dz,muxf,zz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  xm(l) = .5d0*(grd_xarr(l+1) + grd_xarr(l))
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  muxf(zz) = ptcl%x*sin(ptcl2%muz+zz)/sin(ptcl2%muz)
  !muyf(zz) = ptcl%x*sin(zz)/sin(ptcl2%muz)

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
  e0 => ptcl%e0
  wl => ptcl%wl

  mux => ptcl2%mux
  muy => ptcl2%muy
  muz => ptcl2%muz

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
  eamp = 0d0

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
!write(0,*) ptcl2%ipart,ptcl2%istep,ptcl2%idist,x,z,mu,om,mux,x*sin(om)/sin(z+om),muy,x*sin(z)/sin(z+om)
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen

!-- census distance
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
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
!    stop 'transport2: dbx nan'
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
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
!-- making greater than dcen
     dthm = far
  endif
  if(dthm/=dthm) then
!    stop 'transport2: dthm nan'
     ierr = 2
     return
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
!-- making greater than dcen
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
!-- making greater than dcen
     dcol = far
  endif
  if(dcol/=dcol) then
!    stop 'transport2: dthm nan'
     ierr = 3
     return
  endif
!
!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = pc_c*(elabfact-wl*grp_wlinv(ig+1))
     if(ddop<0d0) then
        ddop = far
     endif
  else
!-- making greater than dcen
     ddop = far
  endif
!
!-- finding minimum distance
  darr = [dcen,dcol,dthm,ddop,dbx,dby,dbz]
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
!
!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d

!-- tallying energy densities
  if(trn_isimcanlog) then
!-- analog energy density
     eraddens = e*elabfact**2 * &
          d*thelp*cinv*tsp_dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)* &
          min(dx(ix),dy(iy),xm(ix)*dz(iz))*thelp>1d-6) then
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
!       stop 'transport2: invalid energy deposition'
        ierr = 7
        return
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
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

!
!-- census
  if(d==dcen) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
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
!-- checking if escaped domain
     loutx = d==dbx.and.ixnext==grd_nx+1
     louty = d==dby.and.(iynext==grd_ny+1.or.iynext==0)
     if(loutx.or.louty) then
!-- ending particle
        ptcl2%isvacant = .true.
        ptcl2%done = .true.
        ptcl2%lflux = .true.
        return
     endif
  endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- wavelength
        wl = wl*(1d0-dirdotu*cinv)/elabfact
        help = elabfact/(1d0-dirdotu*cinv)
!-- velocity effects accounting
        totevelo=totevelo+e*(1d0-help)
!
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
        ptcl2%isvacant=.true.
        ptcl2%done=.true.
!-- adding comoving energy to deposition energy
        edep = e*elabfact
        return
     else
!-- effective scattering
!
!-- transforming to lab
        if(grd_isvelocity) then
           help = elabfact/(1d0-dirdotu*cinv)
!-- velocity effects accounting
           totevelo=totevelo+e*(1d0-help)
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!
!-- sample group
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
        if(grp_ng==1) then
        elseif(emitlump<.99d0 .or. trn_nolumpshortcut .or. in_puretran) then
           call rnd_r(r1,rndstate)
           ig = emitgroup(r1,ic)
!-- checking if DDMC in new group
           if((grd_cap(ig,ic)+grd_sig(ic)) * &
              min(dx(ix),dy(iy),xm(ix)*dz(iz))*thelp >= trn_tauddmc &
              .and. .not.in_puretran) ptcl2%itype = 2
!-- don't sample, it will end up in the lump anyway
        else
           ptcl2%itype = 2
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(9,ic))
        endif
!
!-- sample wavelength
        if(ptcl2%itype==2) then
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*dirdotu*cinv
!-- energy weight
              e = e*(1d0-dirdotu*cinv)
              e0 = e0*(1d0-dirdotu*cinv)
           endif
           wl = 0d0 !workaround ifort 13.1.3 bug
        else
!-- uniformly in new group
           call rnd_r(r1,rndstate)
           wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
           if(grd_isvelocity) wl = wl*(1d0-dirdotu*cinv)
        endif
     endif

!
!-- x-bound
  elseif(d==dbx) then

     if(ixnext>ix) then
        x = grd_xarr(ix+1)
     else
        x = grd_xarr(ix)
     endif

     l = grd_icell(ixnext,iy,iz)

     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ixnext),dy(iy),xm(ixnext)*dz(iz))*thelp < trn_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ixnext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)-gm*x*cinv * &
                (1d0-gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu-gm*y*cinv*(1d0-gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0-dirdotu*cinv))
!-- x-cosine
           mu0 = sqrt(1d0-mu**2)*cos(om)
!-- amplification factor
           if(.not.trn_noampfact .and. mu0<0d0) then
              help=1d0/abs(mu0)
              help = min(100d0, help) !-- truncate singularity
   !
   !-- velocity effects accounting
              totevelo=totevelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu0))*x*cinv
   !
   !-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*x*cinv
   !-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu0))*x*cinv)
              e=e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu0))*x*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ixnext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
        help = help*(1d0 + 1.5d0*abs(mu0))
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help) then
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
           ic = grd_icell(ix,iy,iz)
        else
           lredir = .true.
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu0 = (ix-ixnext)*max(r1,r2)
           call rnd_r(r1,rndstate)
!-- resampling y-cosine
           mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
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
           endif
        endif
     endif

!
!-- y-bound
  elseif(d==dby) then

     if(iynext>iy) then
        y = grd_yarr(iy+1)
     else
        y = grd_yarr(iy)
     endif

     l = grd_icell(ix,iynext,iz)

     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iynext),xm(ix)*dz(iz))*thelp < trn_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iynext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- y-projection
           mu = (mu-gm*y*cinv*(1d0-gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0-dirdotu*cinv))
!-- amplification factor
           if(.not.trn_noampfact .and. &
                 ((mu<0d0.and.y>0d0).or.(mu>0d0.and.y<0d0))) then
              help=1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              totevelo=totevelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              eamp = e*2d0*0.55d0*max(0d0,help-2d0)*abs(y)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dy(iynext)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
        help = help*(1d0 + 1.5d0*abs(mu))
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iynext
           ic = grd_icell(ix,iy,iz)
        else
           lredir = .true.
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = (iy-iynext)*max(r1,r2)
!-- resampling azimuthal
           call rnd_r(r1,rndstate)
           om = pc_pi2*r1
           if(grd_isvelocity) then
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
           endif
        endif
     endif

!
!-- z-bound
  elseif(d==dbz) then
!-- sanity check
     if(grd_nz==1) then
!       stop 'transport1: invalid z crossing'
        ierr = 14
        return
     endif
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
     else
!       stop 'transport1: invalid iznext'
        ierr = 15
        return
     endif

     l = grd_icell(ix,iy,iznext)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),xm(ix)*dz(iznext)) * &
          thelp<trn_tauddmc .or. in_puretran) then
!-- IMC in adjacent cell
        iz = iznext
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- azimuthal direction angle
           om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                sqrt(1d0-mu**2)*cos(om)-gm*x*cinv * &
                (1d0-gm*dirdotu*cinv/(gm+1d0)))
           if(om<0d0) om=om+pc_pi2
!-- y-projection
           mu = (mu-gm*y*cinv*(1d0-gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0-dirdotu*cinv))
!-- transforming z-cosine to cmf
           xi = sqrt(1d0-mu**2)*sin(om)
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             xm(ix)*dz(iznext)*thelp
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
!-- resampling y-cosing
           mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(xi,sqrt(1d0-xi**2)*cos(pc_pi2*r1))
           if(om<0d0) om=om+pc_pi2
!-- transforming om, mu to lab
           if(grd_isvelocity) then
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
           endif
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
!       stop 'transport2: ddop and no velocity'
        ierr = 11
        return
     endif
     if(ig<grp_ng) then
!-- shifting group
        ig = ig+1
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        call rnd_r(r1,rndstate)
        wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if ((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),dy(iy),xm(ix)*dz(iz))*thelp >= trn_tauddmc &
          .and..not.in_puretran) then
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
!    stop 'transport2: invalid distance'
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

end subroutine transport2
