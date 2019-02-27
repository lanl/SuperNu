!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine particle_advance_gamgrey(nmpi)

!$ use omp_lib
  use miscmod
  use randommod
  use sourcemod
  use particlemod
  use transportmod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use timestepmod
  use totalsmod
  use fluxmod
  implicit none
  integer,intent(in) :: nmpi
!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
!##################################################
  integer :: ierr
  integer :: nhere, nemit, ndmy
  real*8 :: r1, edep, help
  integer :: i,j,k,l, ii, iimpi
  integer :: imu, iom, icold
  integer,pointer :: ic
  integer,pointer :: ix, iy, iz
  real*8,pointer :: x,y,z,mu,om,e,e0
  real*8 :: eta, xi
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, mu1, mu2, gm
  real*8 :: etot,pwr
  real*8 :: om0, mu0, x0, y0, z0
!
  integer :: nvol(grd_ncell)
!
  type(rnd_t) :: rndstate
  integer,save :: iomp=0
!
  real*8,parameter :: basefrac=.1d0
  real*8 :: base,edone,einv,invn,en
  integer :: n, ndone, mpart, npart, ipart
  integer*8 :: nstot,nsavail,nsbase
!
  integer,allocatable :: ipospart(:,:) !(3,npart)
!
  type(packet),target :: ptcl
  type(packet2),target :: ptcl2
!
!-- start clock
  t0 = t_time()

  grd_tally = 0d0
  flx_gamlumtime = 0d0
  flx_gamluminos = 0d0
  flx_gamlumdev = 0d0
  flx_gamlumnum = 0

!-- initializing volume numbers
  nvol = 0

!-- shortcut
  pwr = in_srcepwr

!-- total particle number
  nstot = nmpi*int(src_ns,8)

!-- total energy (converted by pwr)
  etot = sum(grd_emitex**pwr)
  if(etot/=etot) stop 'particle_advance_gamgrey: etot nan'

!-- base (flat,constant) particle number per cell over ALL RANKS
  n = count(grd_emitex>0d0)  !number of cells that emit
  base = dble(nstot)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase

  mpart = nint(10 + 1.1d0*src_ns) !big enough
  allocate(ipospart(3,mpart))

  iimpi = 0
  npart = 0
!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/nstot
  einv = 1d0/etot
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     en = grd_emitex(l)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     nvol(l) = n
     edone = edone + en
     ndone = ndone + n
!-- particle count on this rank
     call sourcenumbers_roundrobin(iimpi,grd_emitex(l), &
        0d0,nvol(l),nemit,nhere,ndmy)
     do ii=1,nhere
        npart = npart + 1
        if(npart>mpart) stop 'particle_advance_gamgrey: npart>npartmax'
        ipospart(:,npart) = [i,j,k]
     enddo !ii
  enddo
  enddo
  enddo


!$omp parallel &
!$omp shared(nvol) &
!$omp private(ptcl,ptcl2,x0,y0,z0,mu0,om0,cmffact,gm,mu1,mu2,eta,xi,labfact,iom,imu, &
!$omp    rndstate,edep,ierr, iomp, &
!$omp    x,y,z,mu,om,e,e0,ix,iy,iz,ic,icold,r1, &
!$omp    i,j,k) &
!$omp reduction(+:grd_tally,flx_gamluminos,flx_gamlumnum, &
!$omp    flx_gamlumdev,flx_gamlumtime)

!-- thread id                                                               
!$ iomp = omp_get_thread_num()
  rndstate = rnd_states(iomp+1)

!
!-- primary particle properties
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic

!-- default values in bounds
  y = grd_yarr(1)
  z = grd_zarr(1)

!$omp do schedule(static,1) !round-robin
  do ipart=1,npart
     i = ipospart(1,ipart)
     j = ipospart(2,ipart)
     k = ipospart(3,ipart)
!-- adopt position (get rid of this copy after merged with wlT branch)
     ix = i
     iy = j
     iz = k
     ic = grd_icell(ix,iy,iz)

!-- x position
     select case(grd_igeom)
     case(1,11)
        call rnd_r(r1,rndstate)
        x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        x = min(x,grd_xarr(i+1))
        x = max(x,grd_xarr(i))
     case(2)
        call rnd_r(r1,rndstate)
        x = sqrt(r1*grd_xarr(i+1)**2 + (1d0-r1)*grd_xarr(i)**2)
!-- must be inside cell
        x = min(x,grd_xarr(i+1))
        x = max(x,grd_xarr(i))
     case(3)
        call rnd_r(r1,rndstate)
        x = r1*grd_xarr(i+1) + (1d0-r1) * grd_xarr(i)
     endselect

!-- y,z position
     call rnd_r(r1,rndstate)
     y = r1*grd_yarr(j+1) + (1d0-r1) * grd_yarr(j)
     call rnd_r(r1,rndstate)
     z = r1*grd_zarr(k+1) + (1d0-r1) * grd_zarr(k)

!-- direction cosine (comoving)
     call rnd_r(r1,rndstate)
     mu0 = 1d0-2d0*r1
     call rnd_r(r1,rndstate)
     om0 = pc_pi2*r1

!-- transform direction
     if(.not.grd_isvelocity) then
         mu = mu0
         om = om0
     else
        select case(grd_igeom)!{{{
        case(1,11)
           x0 = x
           cmffact = 1d0+mu0*x0/pc_c !-- 1+dir*v/c
           mu = (mu0+x0/pc_c)/cmffact
           om = om0
        case(2)
           x0 = x
           y0 = y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           gm = 1d0/sqrt(1d0-(x**2+y**2)/pc_c**2)
!-- om
           om = atan2(sqrt(1d0-mu0**2)*sin(om0) , &
                sqrt(1d0-mu0**2)*cos(om0)+(gm*x/pc_c) * &
                (1d0+gm*(cmffact-1d0)/(gm+1d0)))
           if(om<0d0) om = om+pc_pi2
!-- mu
           mu = (mu0+(gm*y/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                (gm*cmffact)
        case(3)
           x0 = x
           y0 = y
           z0 = z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           mu = (mu0+z0/pc_c)/cmffact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- om
           om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(om<0d0) om = om+pc_pi2
        endselect!}}}
     endif

!-- velocity components in cartesian basis
     if(grd_igeom==1) then
!-- spherical projections
        eta = sqrt(1d0-mu**2)*cos(om)
        xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
        ptcl2%mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
        ptcl2%muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
        ptcl2%muz = mu*y-eta*sqrt(1d0-y**2)
     endif
!-- update invariant direction quantities
     if(grd_igeom==2) then
        ptcl2%mux = x*sin(om)/sin(z+om)  !-- intercept
        ptcl2%muy = x*sin(z)/sin(z+om)  !-- distance to intercept
        ptcl2%muz = pc_pi-(z+om)  !-- direction angle
        if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
        if(ptcl2%muz<0d0) ptcl2%muz = ptcl2%muz+pc_pi2
     endif
!
!-- emission energy per particle
     e = grd_emitex(ic)/nvol(ic)
     if(grd_isvelocity) e = e*cmffact
     e0 = e

!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
     ptcl2%ipart = ipart
     ptcl2%istep = 0
     ptcl2%idist = 0

     ptcl2%stat = 'live'
!
     do while (ptcl2%stat=='live')
        ptcl2%istep = ptcl2%istep + 1
        icold = ic
        call transport_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
!-- tally
        grd_tally(1,icold) = grd_tally(1,icold) + edep

!-- Russian roulette for termination of exhausted particles
        if(e<1d-6*e0 .and. ptcl2%stat=='live' .and. grd_capgam(ic)>0d0) then
           call rnd_r(r1,rndstate)!{{{
           if(r1<0.5d0) then
!-- transformation factor
              if(grd_isvelocity) then
                 select case(grd_igeom)
                 case(1,11)
                    labfact = 1.0d0 - mu*x/pc_c
                 case(2)
                    labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                         cos(om)*x)/pc_c
                 case(3)
                    labfact = 1d0-(mu*z+sqrt(1d0-mu**2) * &
                         (cos(om)*x+sin(om)*y))/pc_c
                 endselect
              else
                 labfact = 1d0
              endif
!
              ptcl2%stat = 'dead'
              grd_tally(1,ic) = grd_tally(1,ic) + e*labfact
           else
              e = 2d0*e
              e0 = 2d0*e0
           endif!}}}
        endif

!-- verify position
        if(ptcl2%stat=='live') then
           if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix) .or. x/=x) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: x not in cell', &
                ix,x,grd_xarr(ix),grd_xarr(ix+1)
           endif
           if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy) .or. y/=y) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: y not in cell', &
                iy,y,grd_yarr(iy),grd_yarr(iy+1)
           endif
           if(z>grd_zarr(iz+1) .or. z<grd_zarr(iz) .or. z/=z) then
              if(ierr==0) ierr = -99
              write(0,*) 'prt_adv_ggrey: z not in cell', &
                iz,z,grd_zarr(iz),grd_zarr(iz+1)
           endif
        endif

!-- check for errors
        if(ierr/=0 .or. ptcl2%istep>1000) then
           write(0,*) 'pagg: ierr,ipart,istep,idist:',ierr,ptcl2%ipart,ptcl2%istep,ptcl2%idist
           write(0,*) 'dist:',ptcl2%dist
           write(0,*) 't:',ptcl%t
           write(0,*) 'ix,iy,iz,ic,ig:',ptcl2%ix,ptcl2%iy,ptcl2%iz,ptcl2%ic,ptcl2%ig
           write(0,*) 'x,y,z:',ptcl%x,ptcl%y,ptcl%z
           write(0,*) 'mu,om:',ptcl%mu,ptcl%om
           write(0,*) 'mux,muy,muz:',ptcl2%mux,ptcl2%muy,ptcl2%muz
           write(0,*)
           if(ierr>0) then
              if(trn_errorfatal) stop 'particle_advance_gg: fatal transport error'
              ptcl2%stat = 'dead'
              exit
           endif
        endif
     enddo
!
!-- outbound luminosity tally
     if(ptcl2%stat=='flux') then
!-- lab frame flux group, polar, azimuthal bin
        imu = binsrch(mu,flx_mu,flx_nmu+1,.false.)
        iom = binsrch(om,flx_om,flx_nom+1,.false.)
!-- observer corrected time
        help=tsp_t+.5d0*tsp_dt
        select case(grd_igeom)
        case(1,11)
           labfact = mu*x/pc_c
        case(2)
           labfact = (mu*y+sqrt(1d0-mu**2) * &
                cos(om)*x)/pc_c
        case(3)
           labfact = (mu*z+sqrt(1d0-mu**2) * &
                (cos(om)*x+sin(om)*y))/pc_c
        endselect
        if(grd_isvelocity) labfact=labfact*tsp_t
        help=help-labfact
!-- tally outbound luminosity        
        flx_gamlumtime(imu,iom) = flx_gamlumtime(imu,iom)+help
        flx_gamluminos(imu,iom) = flx_gamluminos(imu,iom)+e
        flx_gamlumdev(imu,iom) = flx_gamlumdev(imu,iom)+e**2
        flx_gamlumnum(imu,iom) = flx_gamlumnum(imu,iom)+1
     endif

  enddo !ipart
!$omp end do
!
!-- save state
  rnd_states(iomp+1) = rndstate
!$omp end parallel

  tot_sfluxgamma = -sum(flx_gamluminos)

!-- convert to flux per second
  help = 1d0/tsp_dt
  flx_gamluminos = flx_gamluminos*help
  flx_gamlumdev = flx_gamlumdev*help**2

  deallocate(ipospart)

  t1 = t_time()
  call timereg(t_pcktgam, t1-t0)

end subroutine particle_advance_gamgrey
! vim: fdm=marker
