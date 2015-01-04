subroutine transport3(ptcl,ic,ig,isvacant)

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

  logical :: loutx,louty,loutz
  integer :: imu, iom, ihelp
  real*8 :: elabfact, eta, xi
  real*8 :: dtinv, thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop,d
  real*8 :: r1, r2

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)

  ix => ptcl%ix
  iy => ptcl%iy
  iz => ptcl%iz
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl
!
!-- shortcut
  dtinv = 1d0/tsp_dt
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

!-- census distance
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- boundary distances
  if(xi==0d0) then
     dbx = 2d0*pc_c*tsp_dt*thelpinv
  else
     dbx = max((grd_xarr(ix)-x)/xi,(grd_xarr(ix+1)-x)/xi)
  endif
  if(eta==0d0) then
     dby = 2d0*pc_c*tsp_dt*thelpinv
  else
     dby = max((grd_yarr(iy)-y)/eta,(grd_yarr(iy+1)-y)/eta)
  endif
  if(mu==0d0) then
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  else
     dbz = max((grd_zarr(iz)-z)/mu,(grd_zarr(iz+1)-z)/mu)
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     r1 = rnd_r(rnd_state)
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ic))
  else
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ic))
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     r1 = rnd_r(rnd_state)
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = pc_c*(elabfact-wl*grp_wlinv(ig+1))
     if(ddop<0d0) then
        ddop = 2d0*pc_c*tsp_dt*thelpinv
     endif
  else
     ddop = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dcen,dbx,dby,dbz,dthm,dcol,ddop)
  if(d<0d0) then
     write(*,*) dcen,dbx,dby,dbz,dthm,dcol,ddop,ix,iy,iz,x,y,z,xi,eta,mu
     stop 'transport3: negative distance'
  endif

!-- updating position
  x = x + xi*d
  y = y + eta*d
  z = z + mu*d
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
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        grd_eraddens(ic) = grd_eraddens(ic)+e* &
             (1.0d0-exp(-grd_fcoef(ic)*elabfact* &
             grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*elabfact * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ic) = grd_eraddens(ic)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ic) = grd_edep(ic)+e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ic)/=grd_edep(ic)) then
!       write(0,*) e,grd_fcoef(ic),grd_cap(ig,ic),elabfact,d,thelp
        stop 'transport3: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
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
  elseif(any([dbx,dby,dbz]==d)) then
!-- checking if escaped domain
     loutx = d==dbx.and.((xi>=0d0.and.ix==grd_nx).or.(xi<0.and.ix==1))
     louty = d==dby.and.((eta>=0d0.and.iy==grd_ny).or.(eta<0.and.iy==1))
     loutz = d==dbz.and.((mu>=0d0.and.iz==grd_nz).or.(mu<0.and.iz==1))
     if(loutx.or.louty.or.loutz) then
!-- ending particle
        isvacant = .true.
        prt_done = .true.
!-- retrieving lab frame flux group, polar, azimuthal bin
        iom = binsrch(om,flx_om,flx_nom+1)
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
        flx_luminos(ig,imu,iom) = flx_luminos(ig,imu,iom)+e*dtinv
        flx_lumdev(ig,imu,iom) = flx_lumdev(ig,imu,iom)+(e*dtinv)**2
        flx_lumnum(ig,imu,iom) = flx_lumnum(ig,imu,iom)+1
        return
     endif
  endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- lab wavelength
        wl = wl*(1d0-(mu*z+eta*y+xi*x)*cinv)/elabfact
        help = elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
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
        grd_edep(ic) = grd_edep(ic)+e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        r1 = rnd_r(rnd_state)
        ig = emitgroup(r1,ic)
!-- uniformly in new group
        r1 = rnd_r(rnd_state)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- transforming to lab
        if(grd_isvelocity) then
           wl = wl*(1d0-(mu*z+eta*y+xi*x)*cinv)
           help = elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!-- checking for DDMC in new group
        if((grd_cap(ig,ic)+grd_sig(ic)) * &
             min(dx(ix),dy(iy),dz(iz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ic) = grd_methodswap(ic)+1
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(mu*z+eta*y+xi*x)*cinv
!
              e = e*(1d0-(mu*z+eta*y+xi*x)*cinv)
              e0 = e0*(1d0-(mu*z+eta*y+xi*x)*cinv)
              wl = wl/(1d0-(mu*z+eta*y+xi*x)*cinv)
           endif
        endif
     endif

!
!-- x-bound
  elseif(d==dbx) then

     if(xi>=0d0) then
        ihelp = 1
        x = grd_xarr(ix+1)
     else
        ihelp = -1
        x = grd_xarr(ix)
     endif

     l = grd_icell(ix+ihelp,iy,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix+ihelp),dy(iy),dz(iz))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ix+ihelp
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming x-cosine to cmf
           xi = (xi-x*cinv)/elabfact
           if(xi>1d0) then
              xi = 1
           elseif(xi<-1d0) then
              xi = -1
           endif
!-- amplification factor
           if((xi<0d0.and.x>0d0).or.(xi>0d0.and.x<0d0)) then
              help=1d0/abs(xi)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              grd_eamp(ic) = grd_eamp(ic) + &
                 e*2d0*0.55d0*max(0d0,help-2d0)*abs(x)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dx(ix+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(xi))) then
           ptcl%itype = 2
           grd_methodswap(ic) = grd_methodswap(ic)+1
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
           xi = -ihelp*max(r1,r2)
           r1 = rnd_r(rnd_state)
           eta = sqrt(1d0-xi**2)*cos(pc_pi2*r1)
!-- resampling z-cosine
           mu = sqrt(1d0-xi**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(eta,xi)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- y-bound
  elseif(d==dby) then

     if(eta>=0d0) then
        ihelp = 1
        y = grd_yarr(iy+1)
     else
        ihelp = -1
        y = grd_yarr(iy)
     endif

     l = grd_icell(ix,iy+ihelp,iz)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy+ihelp),dz(iz))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iy+ihelp
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming y-cosine to cmf
           eta = (eta-y*cinv)/elabfact
           if(eta>1d0) then
              eta = 1
           elseif(eta<-1d0) then
              eta = -1
           endif
!-- amplification factor
           if((eta<0d0.and.y>0d0).or.(eta>0d0.and.y<0d0)) then
              help=1d0/abs(eta)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              grd_eamp(ic) = grd_eamp(ic) + &
                 e*2d0*0.55d0*max(0d0,help-2d0)*abs(y)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dy(iy+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(eta))) then
           ptcl%itype = 2
           grd_methodswap(ic) = grd_methodswap(ic)+1
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
           eta = -ihelp*max(r1,r2)
           r1 = rnd_r(rnd_state)
           xi = sqrt(1d0-eta**2)*cos(pc_pi2*r1)
!-- resampling z-cosine
           mu = sqrt(1d0-eta**2)*sin(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(eta,xi)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- z-bound
  elseif(d==dbz) then

     if(mu>=0d0) then
        ihelp = 1
        z = grd_zarr(iz+1)
     else
        ihelp = -1
        z = grd_zarr(iz)
     endif

     l = grd_icell(ix,iy,iz+ihelp)
     if((grd_cap(ig,l)+grd_sig(l)) * &
          min(dx(ix),dy(iy),dz(iz+ihelp))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iz = iz+ihelp
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
!-- transforming z-cosine to cmf
           mu = (mu-z*cinv)/elabfact
           if(mu>1d0) then
              mu = 1
           elseif(mu<-1d0) then
              mu = -1
           endif
!-- amplification factor
           if((mu<0d0.and.z>0d0).or.(mu>0d0.and.z<0d0)) then
              help=1d0/abs(mu)
              help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
              grd_eamp(ic) = grd_eamp(ic) + &
                 e*2d0*0.55d0*max(0d0,help-2d0)*abs(x)*cinv
!-- apply limited correction to the particle
              help = min(2d0,help)
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
           endif
        endif
        help = (grd_cap(ig,l)+grd_sig(l)) * &
             dz(iz+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
        !-- sampling
        r1 = rnd_r(rnd_state)
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl%itype = 2
           grd_methodswap(ic) = grd_methodswap(ic)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iz = iz + ihelp
           ic = grd_icell(ix,iy,iz)
        else
           r1 = rnd_r(rnd_state)
           r2 = rnd_r(rnd_state)
!-- resampling z-cosine
           mu = -ihelp*max(r1,r2)
           r1 = rnd_r(rnd_state)
!-- resampling azimuthal
           om = pc_pi2*r1
           xi = sqrt(1d0-mu**2)*cos(om)
           eta = sqrt(1d0-mu**2)*sin(om)
           if(grd_isvelocity) then
!-- transforming mu to lab
              mu=(mu+z*cinv)/(1d0+(x*xi+y*eta+z*mu)*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              om = atan2(eta+y*cinv,xi+x*cinv)
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) stop 'transport3: ddop and no velocity'
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
     if((grd_sig(ic)+grd_cap(ig,ic)) * &
          min(dx(ix),dy(iy),dz(iz))*thelp >= prt_tauddmc &
          .and..not.in_puretran) then
        ptcl%itype = 2
        grd_methodswap(ic) = grd_methodswap(ic)+1
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
     stop 'transport3: invalid distance'
  endif

end subroutine transport3
