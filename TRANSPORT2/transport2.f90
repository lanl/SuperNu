subroutine transport2(ptcl,ic,ig,isvacant)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  use totalsmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  integer,intent(inout) :: ic, ig
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer,external :: emitgroup

  logical :: loutx,louty
  integer :: imu, ihelp
  real*8 :: elabfact, dirdotu, mu0, gm
  real*8 :: dtinv, thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,ddop,d
  real*8 :: darr(6)
  real*8 :: rold, zold, omold
  real*8 :: r1, r2
!-- distance out of physical reach
  real*8 :: far

  integer,pointer :: ix, iy
  integer,parameter :: iz=1
  real*8,pointer :: x,y,mu,om,e,e0,wl
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
  wl => ptcl%wl
!
!-- shortcut
  dtinv = 1d0/tsp_dt
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
     if(abs(sin(om))<grd_xarr(ix)/x .and. cos(om)<0d0) then
!-- inner boundary
        dbx = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbx = 0d0
     else
!-- outer boundary
        dbx = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
     endif
  endif
  if(dbx/=dbx) stop 'transport2: dbx nan'

!-- to y-bound
  if(mu>0d0) then
     dby = (grd_yarr(iy+1)-y)/mu
  elseif(mu<0d0) then
     dby = (grd_yarr(iy)-y)/mu
  else
!-- making greater than dcen
     dby = far
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     r1 = rnd_r(rnd_state)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
!-- making greater than dcen
     dthm = far
  endif
  if(dthm/=dthm) stop 'transport2: dthm nan'
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
!-- making greater than dcen
     dcol = far
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ic))
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
!-- making greater than dcen
     dcol = far
  endif
  if(dcol/=dcol) stop 'transport2: dthm nan'
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
  darr = [dcen,dbx,dby,dthm,dcol,ddop]
  if(any(darr/=darr) .or. any(darr<0d0)) then
     write(0,*) darr
     write(*,*) ix,iy,x,y,mu,om
     stop 'transport2: invalid distance'
  endif
  d = minval(darr)
!
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
     write(*,*) d, x, rold, omold, om, mu
     stop 'transport2: om is nan'
  endif
!
!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     grd_eraddens(ic)=grd_eraddens(ic)+e*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)* &
          min(dx(ix),dy(iy))*thelp>1d-6) then
        grd_eraddens(ic) = grd_eraddens(ic)+e* &
             (1.0d0-exp(-grd_fcoef(ic)*elabfact* &
             grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*elabfact * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ic)=grd_eraddens(ic)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ic)=grd_edep(ic)+e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ic)/=grd_edep(ic)) then
        stop 'transport2: invalid energy deposition'
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
     prt_done = .true.
     grd_numcensus(ic) = grd_numcensus(ic)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rnd_r(rnd_state)
     mu = 1d0 - 2d0*r1
     r1 = rnd_r(rnd_state)
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
     loutx = d==dbx.and.(cos(om)>=0d0.and.ix==grd_nx)
     louty = d==dby.and.((mu>=0d0.and.iy==grd_ny).or.(mu<0.and.iy==1))
     if(loutx.or.louty) then
!-- ending particle
        isvacant = .true.
        prt_done = .true.
!-- retrieving lab frame flux group, polar bin
        imu = binsrch(mu,flx_mu,flx_nmu+1)
        ig = binsrch(wl,flx_wl,flx_ng+1)
!-- checking group bounds
        if(ig>flx_ng.or.ig<1) then
           if(ig>flx_ng) then
              ig=flx_ng
           else
              ig=1
           endif
        endif
!-- tallying outbound energy
        tot_eout = tot_eout+e
!-- tallying outbound luminosity
        flx_luminos(ig,imu,1) = flx_luminos(ig,imu,1)+e*dtinv
        flx_lumdev(ig,imu,1) = flx_lumdev(ig,imu,1)+(e*dtinv)**2
        flx_lumnum(ig,imu,1) = flx_lumnum(ig,imu,1)+1
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
        tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
        e = e*help
        e0 = e0*help
     endif

!
!-- effective collision
  elseif(d==dcol) then
     r1 = rnd_r(rnd_state)
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ic)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ic) = grd_edep(ic) + e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        r1 = rnd_r(rnd_state)
        if(grp_ng>1) ig = emitgroup(r1,ic)
!-- uniformly in new group
        r1 = rnd_r(rnd_state)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- transforming to lab
        if(grd_isvelocity) then
           wl = wl*(1d0-dirdotu*cinv)
           help = elabfact/(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!-- checking if DDMC in new group
        if((grd_cap(ig,ic)+grd_sig(ic)) * &
             min(dx(ix),dy(iy))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*dirdotu*cinv
!
              e = e*(1d0-dirdotu*cinv)
              e0 = e0*(1d0-dirdotu*cinv)
              wl = wl/(1d0-dirdotu*cinv)
           endif
        endif
     endif

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

     l = grd_icell(ix+ihelp,iy,iz)

     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix+ihelp),dy(iy))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ix+ihelp
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
        endif
!-- x-cosine
        mu0 = sqrt(1d0-mu**2)*cos(om)
!-- amplification factor
        if(grd_isvelocity.and.mu0<0d0) then
           help=1d0/abs(mu0)
           help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
           tot_evelo=tot_evelo-e*2d0 * &
                (0.55d0*help-1.25d0*abs(mu0))*x*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
           grd_eamp(ic) = grd_eamp(ic) + &
              e*2d0*0.55d0*max(0d0,help-2d0)*x*cinv
!-- apply limited correction to the particle
           help = min(2d0,help)
           e0=e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu0))*x*cinv)
           e=e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu0))*x*cinv)
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ix+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(mu0))) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ix + ihelp
           ic = grd_icell(ix,iy,iz)
        else
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
           mu0 = -ihelp*max(r1,r2)
           r1 = rnd_r(rnd_state)
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

     if(mu>=0d0) then
        ihelp = 1
        y = grd_yarr(iy+1)
     else
        ihelp = -1
        y = grd_yarr(iy)
     endif

     l = grd_icell(ix,iy+ihelp,iz)

     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy+ihelp))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iy+ihelp
        ic = grd_icell(ix,iy,iz)
     else
        if(grd_isvelocity) then
           gm = 1d0/sqrt(1d0-(x**2+y**2)*cinv**2)
!-- y-projection
           mu = (mu-gm*y*cinv*(1d0-gm*dirdotu*cinv/(1d0+gm))) / &
                (gm*(1d0-dirdotu*cinv))
!-- amplification factor
           if((mu<0d0.and.y>0d0).or.(mu>0d0.and.y<0d0)) then
              help=1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              grd_eamp(ic) = grd_eamp(ic) + &
                 e*2d0*0.55d0*max(0d0,help-2d0)*abs(y)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(y)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dy(iy+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl%itype = 2
           grd_methodswap(ic)=grd_methodswap(ic)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iy + ihelp
           ic = grd_icell(ix,iy,iz)
        else
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
           mu = -ihelp*max(r1,r2)
!-- resampling azimuthal
           r1 = rnd_r(rnd_state)
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
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) stop 'transport2: ddop and no velocity'
     if(ig<grp_ng) then
!-- shifting group
        ig = ig+1
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        r1 = rnd_r(rnd_state)
        wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if ((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),dy(iy))*thelp >= prt_tauddmc &
          .and..not.in_puretran) then
        ptcl%itype = 2
        grd_methodswap(ic)=grd_methodswap(ic)+1
        if(grd_isvelocity) then
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  else
     stop 'transport2: invalid distance'
  endif

end subroutine transport2
