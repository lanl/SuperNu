subroutine particle_advance(lmpi0)

!$ use omp_lib
  use randommod
  use transportmod
  use miscmod
  use particlemod
  use timestepmod
  use totalsmod
  use groupmod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use countersmod
  use fluxmod
  use sourcemod
  implicit none
  logical,intent(in) :: lmpi0
!
!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
!##################################################
  real*8,parameter :: cinv=1d0/pc_c
  integer :: i
  integer :: nstepddmc, nstepimc, nmethodswap, ncensimc, ncensddmc, npckt, nflux
  real*8 :: r1, x1, x2, thelp, help, tau
! integer :: irl,irr
! real*8 :: xx0, bmax
! real*8 :: uul, uur, uumax, r0,r2,r3
  integer :: ipart
  integer, pointer :: ig, ic, icorig
  integer, pointer :: ix, iy, iz
  real*8, pointer :: x,y,z, mu, e, e0, wl, om
  integer :: iom, imu
  real*8 :: eta, xi
  real*8 :: t0,t1  !timing
  real*8 :: labfact, mu1, mu2
!
  integer :: icold, ierr
  real*8 :: edep, eraddens, eamp
!-- specint cache
  integer :: nstepmax, ndist(-3:7)
!
  type(packet),target :: ptcl
  type(packet2),target :: ptcl2
  type(grp_t_cache),target :: cache
  real*8,target :: specarr(grp_ng)
  integer,target :: glumps(grp_ng)
  logical,target :: llumps(grp_ng)
!
  type(rnd_t) :: rndstate
  integer,save :: iomp=0
!
!-- statement function
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)

!
!-- init timers
  t_pckt_stat = (/1d30, 0d0, 0d0/) !min,mean,max

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- energy tallies
  tot_erad = 0d0

  grd_tally = 0d0
  if(.not.trn_noampfact) grd_eamp = 0d0

  flx_luminos = 0d0
  flx_lumdev = 0d0
  flx_lumnum = 0

!-- optional grid tally
  if(in_io_dogrdtally) then
     grd_methodswap = 0
     grd_numcensimc = 0
     grd_numcensddmc = 0
     grd_numfluxorig = 0
  endif

!-- Propagate all particles that are not considered vacant
  npckt = 0
  nflux = 0
  nstepddmc = 0
  nstepimc = 0
  nmethodswap = 0
  ncensimc = 0
  ncensddmc = 0

  nstepmax = 0
  ndist = 0

!$omp parallel &
!$omp shared(thelp) &
!$omp private(ptcl,ptcl2,cache, &
!$omp    specarr,glumps,llumps, &
!$omp    x,y,z,mu,om,wl,e,e0,icorig,ix,iy,iz,ic,ig,icold,r1, &
!$omp    i,t0,t1,x1,x2,help,tau, &
!$omp    mu1,mu2,eta,xi,labfact,iom,imu, &
!$omp    rndstate,edep,eraddens,eamp,ierr, iomp) &
!$omp reduction(+:grd_tally, &
!$omp    tot_evelo,tot_erad,tot_eout, &
!$omp    npckt,nflux,nstepddmc,nstepimc,nmethodswap,ncensimc,ncensddmc,ndist) &
!$omp reduction(max:nstepmax)

!-- thread id                                                               
!$ iomp = omp_get_thread_num()

!-- each thread uses its own rnd stream
  rndstate = rnd_states(iomp+1)

!$!-- thread timing
!$ if(.false.) then!{{{
      t0 = t_time() !serial version
!$ else
!$    t0 = omp_get_wtime()
!$ endif!}}}
!
!-- assigning pointers to corresponding particle properties
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  wl => ptcl%wl
  e => ptcl%e
  e0 => ptcl%e0
  icorig => ptcl%icorig
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
!-- alloc
  cache%specarr => specarr
  cache%glumps => glumps
  cache%llumps => llumps
  cache%ic = 0

!$omp do schedule(static,1) !round-robin
  do ipart=1,prt_npartmax
     ! Checking vacancy
     if(prt_isvacant(ipart)) cycle
     ptcl2%isvacant = .false.
!
!-- active particle
     ptcl = prt_particles(ipart) !copy properties out of array
     ptcl2%ipart = ipart
     npckt = npckt+1

!-- cell position
     ix = binsrch(x,grd_xarr,grd_nx+1,.false.)
     iy = binsrch(y,grd_yarr,grd_ny+1,.false.)
     iz = binsrch(z,grd_zarr,grd_nz+1,.false.)
!-- cell pointer
     ic = grd_icell(ix,iy,iz)

!-- group index
     ig = binsrch(wl,grp_wl,grp_ng+1,.false.) !co-moving frame

!
!-- determine particle type
     select case(grd_igeom)
     case(11)
        help = thelp*dx(ix)
     case(1)
        help = thelp*min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) 
     case(2)
        help = thelp*min(dx(ix),dy(iy))
     case(3)
        help = thelp*min(dx(ix),dy(iy),dz(iz))
     endselect
!-- IMC or DDMC
     tau = (grd_sig(ic)+grd_cap(ig,ic))*help
     if(in_puretran .or. tau<prt_tauddmc) then
        ptcl2%itype = 1 !IMC
     else
        ptcl2%itype = 2 !DDMC
        icorig = 0
     endif

!
!-- transform IMC particle into lab frame
     if(grd_isvelocity.and.ptcl2%itype==1) then
        select case(grd_igeom)
        case(1,11)
           labfact = 1d0-x*mu/pc_c
        case(2)
           labfact = 1d0-(mu*y + sqrt(1d0-mu**2) * cos(om)*x)/pc_c
        case(3)
           help = sqrt(1d0-mu**2)
           mu1 = help*cos(om)
           mu2 = help*sin(om)
           labfact = 1d0-(mu*z + mu1*x + mu2*y)/pc_c
        endselect
!-- transform into lab frame
        wl = wl*labfact
        e = e/labfact
        e0 = e0/labfact
     endif


!-- First portion of operator split particle velocity position adjustment
     if(grd_isvelocity.and.ptcl2%itype==1) then
        call advection(.true.,ptcl,ptcl2) !procedure pointer to advection[123]
     endif

!-- velocity components in cartesian basis
     if(grd_igeom==1 .and. ptcl2%itype==1) then
!-- spherical projections
        eta = sqrt(1d0-mu**2)*cos(om)
        xi = sqrt(1d0-mu**2)*sin(om)
!-- planar projections (invariant until collision)
        ptcl2%mux = mu*sqrt(1d0-y**2)*cos(z)+eta*y*cos(z)-xi*sin(z)
        ptcl2%muy = mu*sqrt(1d0-y**2)*sin(z)+eta*y*sin(z)+xi*cos(z)
        ptcl2%muz = mu*y-eta*sqrt(1d0-y**2)
     endif

!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
!Calling either diffusion or transport depending on particle type (ptcl2%itype)
     ptcl2%istep = 0
     ptcl2%idist = 0

     ptcl2%done = .false.
     ptcl2%lflux = .false.
     ptcl2%lcens = .false.

     do while (.not.ptcl2%done)
        ptcl2%istep = ptcl2%istep + 1
        icold = ic
        if(ptcl2%itype==1 .or. in_puretran) then
           nstepimc = nstepimc + 1
           call transport(ptcl,ptcl2,rndstate,edep,eraddens,eamp,tot_evelo,ierr)
           if(ptcl2%itype/=1) then
              nmethodswap = nmethodswap + 1
              if(in_io_dogrdtally) grd_methodswap(icold) = grd_methodswap(icold) + 1
              icorig = 0
           endif
!-- tally eamp
           if(.not.trn_noampfact) grd_eamp(icold) = grd_eamp(icold) + eamp
        else
           nstepddmc = nstepddmc + 1
           call diffusion(ptcl,ptcl2,cache,rndstate,edep,eraddens,tot_evelo,ierr)
           if(ptcl2%itype==1) then
              nmethodswap = nmethodswap + 1
              if(in_io_dogrdtally) grd_methodswap(icold) = grd_methodswap(icold) + 1
              icorig = ptcl2%ic
           endif
        endif
        i = max(-3,ptcl2%idist)
        ndist(i) = ndist(i) + 1
!-- tally rest
        grd_tally(:,icold) = grd_tally(:,icold) + [edep,eraddens]
!
!-- outbound luminosity tally
        if(ptcl2%lflux) then
           nflux = nflux + 1
           if(icorig>0 .and. in_io_dogrdtally) grd_numfluxorig(icorig) = grd_numfluxorig(icorig) + 1
           tot_eout = tot_eout+e
!-- retrieving lab frame flux group, polar, azimuthal bin
           ig = binsrch(wl,flx_wl,flx_ng+1,.false.)
           imu = binsrch(mu,flx_mu,flx_nmu+1,.false.)
           iom = binsrch(om,flx_om,flx_nom+1,.false.)
!-- tallying outbound luminosity
           flx_luminos(ig,imu,iom) = flx_luminos(ig,imu,iom)+e
           flx_lumdev(ig,imu,iom) = flx_lumdev(ig,imu,iom)+e**2
           flx_lumnum(ig,imu,iom) = flx_lumnum(ig,imu,iom)+1
        endif
!
!-- Russian roulette for termination of exhausted particles
        if(e<1d-6*e0 .and. .not.ptcl2%isvacant .and. &
              grd_capgrey(ic)+grd_sig(ic)>0d0) then
!-- transformation factor!{{{
           if(.not.grd_isvelocity .or. ptcl2%itype==2) then
              labfact = 1d0
           else
              select case(grd_igeom)
              case(1,11)
                 labfact = 1d0 - mu*x/pc_c
              case(2)
                 labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                    cos(om)*x)/pc_c
              case(3)
                 labfact = 1d0-(mu*z+sqrt(1d0-mu**2) * &
                    (cos(om)*x+sin(om)*y))/pc_c
              endselect
           endif
!
           call rnd_r(r1,rndstate)
           if(r1<0.5d0) then
              ptcl2%isvacant = .true.
              ptcl2%done = .true.
              grd_tally(1,ic) = grd_tally(1,ic) + e*labfact
!-- velocity effects accounting
              if(ptcl2%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
           else
!-- weight addition accounted for in external source
              tot_eext = tot_eext + e
!
              e = 2d0*e
              e0 = 2d0*e0
           endif!}}}
        endif
!
!-- tally census
        if(ptcl2%lcens) then
           if(ptcl2%itype==1) then
              ncensimc = ncensimc + 1
              if(in_io_dogrdtally) grd_numcensimc(ic) = grd_numcensimc(ic) + 1
           else
              ncensddmc = ncensddmc + 1
              if(in_io_dogrdtally) grd_numcensddmc(ic) = grd_numcensddmc(ic) + 1
           endif
        endif

!-- verify position
        if(ptcl2%itype==1 .and. .not.ptcl2%done) then
           if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then!{{{
              write(0,*) 'prt_adv: x not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu, &
                 ptcl2%ipart,ptcl2%istep,ptcl2%idist
           endif
           if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
              write(0,*) 'prt_adv: y not in cell',iy,y,grd_yarr(iy),grd_yarr(iy+1),mu, &
                 ptcl2%ipart,ptcl2%istep,ptcl2%idist
           endif
           if(z>grd_zarr(iz+1) .or. z<grd_zarr(iz)) then
              write(0,*) 'prt_adv: z not in cell',iz,z,grd_zarr(iz),grd_zarr(iz+1),mu, &
                 ptcl2%ipart,ptcl2%istep,ptcl2%idist
           endif!}}}
        endif

!-- check exit status
        if(ierr/=0) then
           write(0,*) 'pa: ierr,ipart,istep,idist:',ierr,ptcl2%ipart,ptcl2%istep,ptcl2%idist
           write(0,*) 'dist:',ptcl2%dist
           write(0,*) 'tddmc:',tau
           write(0,*) 'ix,iy,iz,ic,ig:',ptcl2%ix,ptcl2%iy,ptcl2%iz,ptcl2%ic,ptcl2%ig
           write(0,*) 'x,y,z:',ptcl%x,ptcl%y,ptcl%z
           write(0,*) 'mu,om:',ptcl%mu,ptcl%om
           if(ierr>0) stop 'particle_advance: fatal transport error'
        endif
     enddo

!-- max step counter
     nstepmax = max(nstepmax,ptcl2%istep)

!-- continue only if particle not vacant
     prt_isvacant(ipart) = ptcl2%isvacant
     if(ptcl2%isvacant) cycle

!
!-- Redshifting DDMC particle energy weights and wavelengths
     if(ptcl2%itype==2 .and. grd_isvelocity) then
!-- redshifting energy weight!{{{
        tot_evelo = tot_evelo + e*(1d0-exp(-tsp_dt/tsp_t))
        e = e*exp(-tsp_dt/tsp_t)
        e0 = e0*exp(-tsp_dt/tsp_t)
        !
!
!-- find group
        ig = binsrch(wl,grp_wl,grp_ng+1,.false.)
!
        call rnd_r(r1,rndstate)
        x1 = grd_cap(ig,ic)
        x2 = grp_wl(ig)/(pc_c*tsp_t*(grp_wl(ig+1)-grp_wl(ig)))
        if(r1<x2/(x1+x2)) then
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
           wl = wl*exp(tsp_dt/tsp_t)
        endif
        !!}}}
     endif

     if(grd_isvelocity.and.ptcl2%itype==1) then
        call advection(.false.,ptcl,ptcl2) !procedure pointer to advection[123]
     endif

!-- radiation energy at census
     tot_erad = tot_erad + e


!-- sample position, direction, and wavelength of a censused DDMC particle
!-- It can be an IMC particle in the next time step
     if(ptcl2%itype==2) then
!
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
!
!-- sample position and direction
        select case(grd_igeom)
!-- 1D spherical
        case(11)
!-- sampling position uniformly!{{{
           call rnd_r(r1,rndstate)
           x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
!-- sampling angle isotropically
           call rnd_r(r1,rndstate)
           mu = 1.0 - 2.0*r1 !}}}
!-- 3D spherical
        case(1)
!-- sampling position uniformly!{{{
           call rnd_r(r1,rndstate)
           x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
           call rnd_r(r1,rndstate)
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
           call rnd_r(r1,rndstate)
           z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           y = min(y,grd_yarr(iy+1))
           y = max(y,grd_yarr(iy))
           z = min(z,grd_zarr(iz+1))
           z = max(z,grd_zarr(iz))
!-- sampling angle isotropically
           call rnd_r(r1,rndstate)
           mu = 1.0 - 2.0*r1
           call rnd_r(r1,rndstate)
           om = pc_pi2*r1!}}}
!-- 2D
        case(2)
!-- sampling position uniformly!{{{
           call rnd_r(r1,rndstate)
           x = sqrt(r1*grd_xarr(ix+1)**2 + (1d0-r1)*grd_xarr(ix)**2)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           call rnd_r(r1,rndstate)
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
!-- sampling direction values
           call rnd_r(r1,rndstate)
           om = pc_pi2*r1
           call rnd_r(r1,rndstate)
           mu = 1d0 - 2d0*r1!}}}
!-- 3D
        case(3)
!-- sampling position uniformly !{{{
           call rnd_r(r1,rndstate)
           x = r1*grd_xarr(ix+1)+(1d0-r1)*grd_xarr(ix)
           call rnd_r(r1,rndstate)
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
           call rnd_r(r1,rndstate)
           z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           y = min(y,grd_yarr(iy+1))
           y = max(y,grd_yarr(iy))
           z = min(z,grd_zarr(iz+1))
           z = max(z,grd_zarr(iz))
!-- sampling direction values
           call rnd_r(r1,rndstate)
           om = pc_pi2*r1
           call rnd_r(r1,rndstate)
           mu = 1d0 - 2d0*r1 !}}}
        endselect
!
        if(grd_isvelocity) call direction2lab(x,y,z,mu,om)
     endif

!
!-- transform IMC particle energy to comoving frame for storage
     if(grd_isvelocity.and.ptcl2%itype==1) then
        select case(grd_igeom)
!-- [123]D spherical
        case(1,11)
           labfact = 1d0-x*mu/pc_c
!-- 2D
        case(2)
           labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * cos(om)*x)/pc_c
!-- 3D
        case(3)
           mu1 = sqrt(1d0-mu**2)*cos(om)
           mu2 = sqrt(1d0-mu**2)*sin(om)
           labfact = 1d0-(mu*z+mu1*x+mu2*y)/pc_c
        endselect

!-- apply inverse labfact for symmetry (since gamma factor is missing)
        wl = wl/labfact
        e = e*labfact
        e0 = e0*labfact
     endif

!
!-- save particle properties
     prt_particles(ipart) = ptcl

  enddo !ipart
!$omp end do nowait

! write(0,*) iomp,nstepmax, ndist

!-- save state
  rnd_states(iomp+1) = rndstate

!-- thread timing
!$ if(.false.) then
      t1 = t_time() !serial version
!$ else
!$    t1 = omp_get_wtime()
!$ endif

!$omp critical
  t1 = t1 - t0
  t_pckt_stat = (/min(t_pckt_stat(1),t1), t_pckt_stat(2)+t1/in_nomp, &
     max(t_pckt_stat(3),t1)/) !min,mean,max
!$omp end critical
!$omp end parallel

!-- print distance counters
! if(lmpi0) write(6,'(11(i6,"k"))',advance='no') ndist(-3:7)/1000

  src_nflux = nflux

  tot_sflux = -sum(flx_luminos)

!-- convert to flux per second
  help = 1d0/tsp_dt
  flx_luminos = flx_luminos*help
  flx_lumdev = flx_lumdev*help**2

  call counterreg(ct_nptransport, npckt)
  call counterreg(ct_npstepimc, nstepimc)
  call counterreg(ct_npstepddmc, nstepddmc)
  call counterreg(ct_npstepmax, nstepmax)
  call counterreg(ct_npcensimc, ncensimc)
  call counterreg(ct_npcensddmc, ncensddmc)
  call counterreg(ct_npmethswap, nmethodswap)
  call counterreg(ct_npflux, nflux)


end subroutine particle_advance
