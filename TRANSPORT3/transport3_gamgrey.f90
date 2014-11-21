subroutine transport3_gamgrey(ptcl)

  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch

  logical :: loutx,louty,loutz
  integer :: imu, iom, ihelp
  real*8 :: elabfact, eta, xi
  real*8 :: dtinv, thelp, thelpinv
  real*8 :: dcol,dbx,dby,dbz,d
  real*8 :: r1

  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,mu,om,e,e0
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
!-- boundary distances
  if(xi==0d0) then
     dbx = 2d0*pc_c*tsp_dt*thelpinv
  else
     if((grd_xarr(ix)-x)/xi>0d0.and.(grd_xarr(ix+1)-x)/xi>0d0) stop &
          'transport3_gamgrey: x val out of cell'
     dbx = max((grd_xarr(ix)-x)/xi,(grd_xarr(ix+1)-x)/xi)
  endif
  if(eta==0d0) then
     dby = 2d0*pc_c*tsp_dt*thelpinv
  else
     if((grd_yarr(iy)-y)/eta>0d0.and.(grd_yarr(iy+1)-y)/eta>0d0) stop &
          'transport3_gamgrey: y val out of cell'
     dby = max((grd_yarr(iy)-y)/eta,(grd_yarr(iy+1)-y)/eta)
  endif
  if(mu==0d0) then
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  else
     if((grd_zarr(iz)-z)/mu>0d0.and.(grd_zarr(iz+1)-z)/mu>0d0) stop &
          'transport3_gamgrey: z val out of cell'
     dbz = max((grd_zarr(iz)-z)/mu,(grd_zarr(iz+1)-z)/mu)
  endif
!
!-- effective collision distance
  if(grd_capgam(ix,iy,iz)<=0d0) then
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*grd_capgam(ix,iy,iz))
  else
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dbx,dby,dbz,dcol)
  if(d<0d0) stop 'transport3_gamgrey: negative distance'

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
     if(grd_capgam(ix,iy,iz)* &
          min(dx(ix),dy(iy),dz(iz))*thelp>1d-6) then
        grd_eraddens(ix,iy,iz) = grd_eraddens(ix,iy,iz)+e* &
             (1.0d0-exp(-elabfact* &
             grd_capgam(ix,iy,iz)*d*thelp))* &
             elabfact/(elabfact * &
             grd_capgam(ix,iy,iz)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ix,iy,iz)=grd_eraddens(ix,iy,iz)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e* &
          (1d0-exp(-grd_capgam(ix,iy,iz)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ix,iy,iz)/=grd_edep(ix,iy,iz)) then
        stop 'transport3_gamgrey: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_capgam(ix,iy,iz) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!-- common manipulations for collisions
  if(d==dcol) then
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
!-- energy weight
        e = e*elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
        e0 = e0*elabfact/(1d0-(mu*z+eta*y+xi*x)*cinv)
     endif
  elseif(any([dbx,dby,dbz]==d)) then
!-- checking if escaped domain
     loutx = d==dbx.and.((xi>=0d0.and.ix==grd_nx).or.(xi<0.and.ix==1))
     louty = d==dby.and.((eta>=0d0.and.iy==grd_ny).or.(eta<0.and.iy==1))
     loutz = d==dbz.and.((mu>=0d0.and.iz==grd_nz).or.(mu<0.and.iz==1))
     if(loutx.or.louty.or.loutz) then
!-- ending particle
        prt_done = .true.
!-- retrieving lab frame flux group, polar, azimuthal bin
        iom = binsrch(om,flx_om,flx_nom+1,0)
        imu = binsrch(mu,flx_mu,flx_nmu+1,0)
!-- tallying outbound luminosity
        flx_gamluminos(imu,iom) = flx_gamluminos(imu,iom)+e*dtinv
        flx_gamlumdev(imu,iom) = flx_gamlumdev(imu,iom)+(e0*dtinv)**2
        flx_gamlumnum(imu,iom) = flx_gamlumnum(imu,iom)+1
        return
     endif
  endif

!-- effective collision
  if(d==dcol) then
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog) then
!-- effective absorption
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,iz)=grd_edep(ix,iy,iz)+e*elabfact
        return
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
!-- IMC in adjacent cell
     ix = ix+ihelp
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
!-- IMC in adjacent cell
     iy = iy+ihelp
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
!-- IMC in adjacent cell
     iz = iz+ihelp
  else
     stop 'transport3_gamgrey: invalid distance'
  endif

end subroutine transport3_gamgrey
