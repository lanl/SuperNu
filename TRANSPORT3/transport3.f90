subroutine transport3(ptcl,isvacant)

  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  use totalsmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch

  logical :: loutx,louty,loutz
  integer :: ig, imu, iom, ihelp
  real*8 :: elabfact, eta, xi
  real*8 :: dtinv, thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop,d
  real*8 :: r1, r2, denom2

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
!
!-- looking up initial group
  ig = binsrch(wl/elabfact,grd_wl,grd_ng+1,in_ng)
!-- checking group bounds
  if(ig>grd_ng.or.ig<1) then
     if(ig==grd_ng+1) then
        ig = grd_ng
     elseif(ig==0) then
        ig = 1
     else
        stop 'transport3 (1): particle group invalid'
     endif
  endif

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
  if(grd_sig(ix,iy,iz)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ix,iy,iz))
  else
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ix,iy,iz)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ix,iy,iz))
  elseif(grd_fcoef(ix,iy,iz)<1d0.and.grd_fcoef(ix,iy,iz)>=0d0) then
     r1 = rand()
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ix,iy,iz))*grd_cap(ig,ix,iy,iz))
  else
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grd_ng) then
     ddop = pc_c*(elabfact-wl/grd_wl(ig+1))
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
!
!-- updating transformation factors
  if(grd_isvelocity) then
     elabfact=1d0-(xi*x+eta*y+mu*z)*cinv
  endif

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     grd_eraddens(ix,iy,iz)=grd_eraddens(ix,iy,iz)+e*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz)* &
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        grd_eraddens(ix,iy,iz) = grd_eraddens(ix,iy,iz)+e* &
             (1.0d0-exp(-grd_fcoef(ix,iy,iz)*elabfact* &
             grd_cap(ig,ix,iy,iz)*d*thelp))* &
             elabfact/(grd_fcoef(ix,iy,iz)*elabfact * &
             grd_cap(ig,ix,iy,iz)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ix,iy,iz)=grd_eraddens(ix,iy,iz)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e* &
          (1d0-exp(-grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ix,iy,iz)/=grd_edep(ix,iy,iz)) then
        stop 'transport3: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ix,iy,iz)*grd_cap(ig,ix,iy,iz) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!
!-- census
  if(d==dcen) then
     prt_done = .true.
     grd_numcensus(ix,iy,iz)=grd_numcensus(ix,iy,iz)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rand()
     mu = 1d0 - 2d0*r1
     r1 = rand()
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
        iom = binsrch(om,flx_om,flx_nom+1,0)
        imu = binsrch(mu,flx_mu,flx_nmu+1,0)
        ig = binsrch(wl,flx_wl,flx_ng+1,0)
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
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ix,iy,iz)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        denom2 = 0d0
        r1 = rand()
        do ig=1,grd_ng-1
           if ((r1>=denom2).and.(r1<denom2+grd_emitprob(ig,ix,iy,iz))) exit
           denom2 = denom2+grd_emitprob(ig,ix,iy,iz)
        enddo
!-- uniformly in new group
        r1 = rand()
        wl = 1d0/((1d0-r1)/grd_wl(ig)+r1/grd_wl(ig+1))
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
        if((grd_cap(ig,ix,iy,iz)+grd_sig(ix,iy,iz)) * &
             min(dx(ix),dy(iy),dz(iz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
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

     if((grd_cap(ig,ix+ihelp,iy,iz)+grd_sig(ix+ihelp,iy,iz)) * &
          min(dx(ix+ihelp),dy(iy),dz(iz))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ix+ihelp
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
           if((xi<0d0.and.x>0d0).or.(xi>0d0.and.x<0d0)) then
              help=1d0/abs(xi)
!-- truncate singularity emperially determined
!-- value of 1000 starts to add noise to W7 spectra
              help = min(100d0, help)
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv
!
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(xi))*abs(x)*cinv)
           endif
        endif
        help = (grd_cap(ig,ix+ihelp,iy,iz)+grd_sig(ix+ihelp,iy,iz)) * &
             dx(ix+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rand()
        if (r1 < help*(1d0+1.5d0*abs(xi))) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ix + ihelp
        else
           r1 = rand()
           r2 = rand()
           xi = -ihelp*max(r1,r2)
           r1 = rand()
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

     if((grd_cap(ig,ix,iy+ihelp,iz)+grd_sig(ix,iy+ihelp,iz)) * &
          min(dx(ix),dy(iy+ihelp),dz(iz))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iy+ihelp
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
           if((eta<0d0.and.y>0d0).or.(eta>0d0.and.y<0d0)) then
              help=1d0/abs(eta)
!-- truncate singularity emperially determined
!-- value of 1000 starts to add noise to W7 spectra
              help = min(100d0, help)
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv
!
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(eta))*abs(y)*cinv)
           endif
        endif
        help = (grd_cap(ig,ix,iy+ihelp,iz)+grd_sig(ix,iy+ihelp,iz)) * &
             dy(iy+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rand()
        if (r1 < help*(1d0+1.5d0*abs(eta))) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iy + ihelp
        else
           r1 = rand()
           r2 = rand()
           eta = -ihelp*max(r1,r2)
           r1 = rand()
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

     if((grd_cap(ig,ix,iy,iz+ihelp)+grd_sig(ix,iy,iz+ihelp)) * &
          min(dx(ix),dy(iy),dz(iz+ihelp))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iz = iz+ihelp
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
           if((mu<0d0.and.z>0d0).or.(mu>0d0.and.z<0d0)) then
              help=1d0/abs(mu)
!-- truncate singularity emperially determined
!-- value of 1000 starts to add noise to W7 spectra
              help = min(100d0, help)
!-- velocity effects accounting
              tot_evelo=tot_evelo-e*2d0 * &
                   (0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv
!
              e0=e0*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
              e=e*(1d0 + 2d0*(0.55d0*help-1.25d0*abs(mu))*abs(z)*cinv)
           endif
        endif
        help = (grd_cap(ig,ix,iy,iz+ihelp)+grd_sig(ix,iy,iz+ihelp)) * &
             dz(iz+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
        !-- sampling
        r1 = rand()
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iz = iz + ihelp
        else
           r1 = rand()
           r2 = rand()
!-- resampling z-cosine
           mu = -ihelp*max(r1,r2)
           r1 = rand()
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
     if(ig<grd_ng) then
!-- shifting group
        ig = ig+1
        wl = (grd_wl(ig)+1d-6*(grd_wl(ig+1)-grd_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        r1 = rand()
        wl=1d0/(r1/grd_wl(grd_ng+1) + (1d0-r1)/grd_wl(grd_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz)) * &
          min(dx(ix),dy(iy),dz(iz))*thelp >= prt_tauddmc &
          .and..not.in_puretran) then
        ptcl%itype = 2
        grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
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
