subroutine transport3(ptcl,isvacant)

  use gasgridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
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

  integer :: ig, imu, iom
  real*8 :: elabfact
  real*8 :: dtinv, thelp, thelpinv
  real*8 :: dcen,dcol,dthm,dbx,dby,dbz,ddop,d
  real*8 :: r1, r2, denom2

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,xi,om,ep,ep0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)

  ix => ptcl%zsrc
  iy => ptcl%iy
  iz => ptcl%iz
  x => ptcl%rsrc
  y => ptcl%y
  z => ptcl%z
  xi => ptcl%musrc
  om => ptcl%om
  ep => ptcl%esrc
  ep0 => ptcl%ebirth
  wl => ptcl%wlsrc
!
!-- shortcut
  dtinv = 1d0/tsp_dt
!-- projections
  eta = sqrt(1d0-xi**2)*sin(om)
  mu = sqrt(1d0-xi**2)*cos(om)
!
!-- setting vel-grid helper variables
  if(gas_isvelocity) then
!-- calculating initial transformation factors
     elabfact=1d0-(xi*z+eta*y+mu*x)*cinv
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
  ig = binsrch(wl/elabfact,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
  if(ig>gas_ng.or.ig<1) then
     if(ig==gas_ng+1) then
        ig = gas_ng
     elseif(ig==0) then
        ig = 1
     else
        stop 'transport3 (1): particle group invalid'
     endif
  endif

!-- census distance
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%tsrc)*thelpinv)
!
!-- boundary distances
  if(mu==0d0) then
     dbx = 2d0*pc_c*tsp_dt*thelpinv
  else
     dbx = max((gas_xarr(ix)-x)/mu,(gas_xarr(ix+1)-x)/mu)
  endif
  if(eta==0d0) then
     dby = 2d0*pc_c*tsp_dt*thelpinv
  else
     dby = max((gas_yarr(iy)-y)/eta,(gas_yarr(iy+1)-y)/eta)
  endif
  if(xi==0d0) then
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  else
     dbz = max((gas_zarr(iz)-z)/xi,(gas_zarr(iz+1)-z)/xi)
  endif
!
!-- Thomson scattering distance
  if(gas_sig(ix,iy,iz)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*gas_sig(ix,iy,iz))
  else
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- effective collision distance
  if(gas_cap(ig,ix,iy,iz)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*gas_cap(ig,ix,iy,iz))
  elseif(gas_fcoef(ix,iy,iz)<1d0.and.gas_fcoef(ix,iy,iz)>=0d0) then
     r1 = rand()
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-gas_fcoef(ix,iy,iz))*gas_cap(ig,ix,iy,iz))
  else
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- Doppler shift distance
  if(gas_isvelocity.and.ig<gas_ng) then
     ddop = pc_c*(elabfact-wl/gas_wl(ig+1))
     if(ddop<0d0) then
        ddop = 2d0*pc_c*tsp_dt*thelpinv
     endif
  else
     ddop = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dcen,dbx,dby,dbz,dthm,dcol,ddop)
  if(any((/dcen,dbx,dby,dyz,dthm,dcol,ddop/)<0d0)) then
     write(*,*) dcen,dbz,dbr,dthm,dcol,ddop
     stop 'transport3: negative distance'
  endif

!-- updating position
  x = x + mu*d
  y = y + eta*d
  z = z + xi*d
!
!-- updating time
  ptcl%tsrc = ptcl%tsrc + thelp*cinv*d
!
!-- updating transformation factors
  if(gas_isvelocity) then
     elabfact=1d0-(mu*x+eta*y+xi*z)*cinv
  endif

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     gas_eraddens(ix,iy,iz)=gas_eraddens(ix,iy,iz)+ep*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(gas_fcoef(ix,iy,iz)*gas_cap(ig,ix,iy,iz)* &
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        gas_eraddens(ix,iy,iz) = gas_eraddens(ix,iy,iz)+ep* &
             (1.0d0-exp(-gas_fcoef(ix,iy,iz)*elabfact* &
             gas_cap(ig,ix,iy,iz)*d*thelp))* &
             elabfact/(gas_fcoef(ix,iy,iz)*elabfact * &
             gas_cap(ig,ix,iy,iz)*pc_c*tsp_dt)
     else
!-- analog energy density
        gas_eraddens(ix,iy,iz)=gas_eraddens(ix,iy,iz)+ep*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- reducing particle energy
     ep = ep*exp(-gas_fcoef(ix,iy,iz)*gas_cap(ig,ix,iy,iz) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!
!-- census
  if(d==dcen) then
     prt_done = .true.
     gas_numcensus(ix,iy,iz)=gas_numcensus(ix,iy,iz)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rand()
     xi = 1d0 - 2d0*r1
     r1 = rand()
     om = pc_pi2*r1
!-- checking velocity dependence
     if(gas_isvelocity) then
        eta = sqrt(1d0-xi**2)*sin(om)
        mu = sqrt(1d0-xi**2)*cos(om)
!-- transforming xi
        xi = (xi+z*cinv)/(1d0+(xi*z+eta*y+mu*x)*cinv)
        if(xi>1d0) then
           xi = 1d0
        elseif(xi<-1d0) then
           xi = -1d0
        endif
!-- transforming om
        om = atan2(eta+y*cinv,mu+x*cinv)
        if(om<0d0) om=om+pc_pi2
!-- x,y lab direction cosines
        eta = sqrt(1d0-xi**2)*sin(om)
        mu = sqrt(1d0-xi**2)*cos(om)
!-- energy weight
        ep = ep*elabfact/(1d0-(xi*z+eta*y+mu*x)*cinv)
        ep0 = ep0*elabfact/(1d0-(xi*z+eta*y+mu*x)*cinv)
     endif

!-- common manipulations for boundary crossings
     if(d==dbx.or.d==dby.or.d==dbz) then

     endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(gas_isvelocity) then
!-- lab wavelength
        wl = wl*(1d0-(xi*z+eta*y+mu*x)*cinv)/elabfact
     endif

!
!-- effective collision
  elseif(d==dcol) then
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=gas_fcoef(ix,iy,iz)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        gas_edep(ix,iy,iz)=gas_edep(ix,iy,iz)+ep*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        denom2 = 0d0
        r1 = rand()
        do ig = 1, gas_ng
           if ((r1>=denom2).and.(r1<denom2+gas_emitprob(ig,ix,iy,iz))) exit
           denom2 = denom2+gas_emitprob(ig,ix,iy,iz)
        enddo
!-- uniformly in new group
        r1 = rand()
        wl = 1d0/((1d0-r1)/gas_wl(ig)+r1/gas_wl(ig+1))
!-- transforming to lab
        if(gas_isvelocity) then
           wl = wl*(1d0-(xi*z+eta*y+mu*x)*cinv)
        endif
!-- checking for DDMC in new group
        if((gas_cap(ig,ix,iy,iz)+gas_sig(ix,iy,iz)) * &
             min(dx(ix),dy(iy),dz(iz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%rtsrc = 2
           gas_methodswap(ix,iy,iz)=gas_methodswap(ix,iy,iz)+1
!-- transforming to cmf
           if(gas_isvelocity) then
              ep = ep*(1d0-(xi*z+eta*y+mu*x)*cinv)
              ep0 = ep0*(1d0-(xi*z+eta*y+mu*x)*cinv)
              wl = wl/(1d0-(xi*z+eta*y+mu*x)*cinv)
           endif
        endif
     endif

!
!-- x-bound
  elseif(d==dbx) then

!
!-- y-bound
  elseif(d==dby) then

!
!-- z-bound
  elseif(d==dbz) then

!
!-- Doppler shift
  elseif(d==ddop) then

  else
     stop 'transport3: invalid distance'
  endif

end subroutine transport3
