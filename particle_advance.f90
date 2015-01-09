subroutine particle_advance

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
  use mpimod
  use fluxmod
  implicit none
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
  integer*8 :: nddmc, nimc, npckt
  real*8 :: r1, x1, x2, thelp, help
! integer :: irl,irr
! real*8 :: xx0, bmax
! real*8 :: uul, uur, uumax, r0,r2,r3
  integer :: ipart
  integer, pointer :: ig, ic
  integer, pointer :: ix, iy, iz
  real*8, pointer :: x,y,z, mu, e, e0, wl, om
  real*8 :: t0,t1  !timing
  real*8 :: labfact, mu1, mu2, om0, mu0
!-- specint cache
  real*8 :: specarr(grp_ng)
  integer :: icell(3)
!
  type(packet),target :: ptcl
  type(packet2),target :: ptcl2
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
!-- assigning pointers to corresponding particle properties
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  wl => ptcl%wl
  e => ptcl%e
  e0 => ptcl%e0
!-- secondary particle properties
  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- energy tallies
  grd_eamp = 0d0
  grd_edep = 0d0
  grd_eraddens = 0d0
  tot_erad = 0d0

  flx_luminos = 0d0
  flx_lumdev = 0d0
  flx_lumnum = 0
  grd_methodswap = 0
  grd_numcensus = 0

  t0 = t_time()
  ! Propagating all particles that are not considered vacant: loop
  npckt = 0
  nddmc = 0
  nimc = 0
  icell = 0
  do ipart=1,prt_npartmax
     ! Checking vacancy
     if(prt_isvacant(ipart)) cycle
     ptcl2%isvacant = .false.
!
!-- active particle
     ptcl = prt_particles(ipart) !copy properties out of array
     ptcl2%ipart = ipart
     npckt = npckt + 1

     ptcl2%done = .false.

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
     select case(in_igeom)
     case(1)
        help = min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) 
     case(2)
        help = min(dx(ix),dy(iy))
     case(3)
        help = min(dx(ix),dy(iy),dz(iz))
     case(11)
        help = dx(ix)
     endselect
!-- IMC or DDMC
     if(in_puretran .or. (grd_sig(ic)+grd_cap(ig,ic))*help*thelp<prt_tauddmc) then
        ptcl2%itype = 1 !IMC
!-- to be transformed
        mu0 = mu
        om0 = om
     else
        ptcl2%itype = 2 !DDMC
     endif

!
!-- transform IMC particle into lab frame
     if(grd_isvelocity.and.ptcl2%itype==1) then
        select case(in_igeom)
        case(1,11)
           labfact = 1d0-x*mu/pc_c
        case(2)
           labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * cos(om)*x)/pc_c
        case(3)
           help = sqrt(1d0-mu**2)
           mu1 = help*cos(om)
           mu2 = help*sin(om)
           labfact = 1d0-(mu*z+mu1*x+mu2*y)/pc_c
        endselect
!-- transform into lab frame
        wl = wl*labfact
        e = e/labfact
        e0 = e0/labfact
     endif


!-- First portion of operator split particle velocity position adjustment
     if((grd_isvelocity).and.(ptcl2%itype==1)) then
        call advection(.true.,ptcl,ptcl2) !procedure pointer to advection[123]
     endif

!     write(*,*) ipart
!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
!Calling either diffusion or transport depending on particle type (ptcl2%itype)
     select case(in_igeom)

!-- 3D spherical
     case(1)
        ptcl2%istep = 0!{{{
        do while ((.not.ptcl2%done).and.(.not.ptcl2%isvacant))
           ptcl2%istep = ptcl2%istep + 1
           if(ptcl2%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport1(ptcl,ptcl2)
           else
              nddmc = nddmc + 1
              call diffusion1(ptcl,ptcl2,icell,specarr)
           endif
!-- verify position
           if(ptcl2%itype==1 .and. .not.ptcl2%done) then
              if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then
                 write(0,*) 'prt_adv: r not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
              endif
              if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
                 write(0,*) 'prt_adv: theta not in cell',iy,y,grd_yarr(iy),grd_yarr(iy+1),mu
              endif
              if(z>grd_zarr(iz+1) .or. z<grd_zarr(iz)) then
                 write(0,*) 'prt_adv: phi not in cell',iz,z,grd_zarr(iz),grd_zarr(iz+1),mu,om
              endif
           endif
!-- Russian roulette for termination of exhausted particles
           if(e<1d-6*e0 .and. .not.ptcl2%isvacant .and. &
                 grd_capgrey(ic)+grd_sig(ic)>0d0) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl2%itype==1) then
                 labfact = 1d0 - mu*x/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 ptcl2%isvacant = .true.
                 ptcl2%done = .true.
                 grd_edep(ic) = grd_edep(ic) + e*labfact
!-- velocity effects accounting
                 if(ptcl2%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo !}}}

!-- 2D
     case(2)
        ptcl2%istep = 0!{{{
        do while ((.not.ptcl2%done).and.(.not.ptcl2%isvacant))
           ptcl2%istep = ptcl2%istep + 1
           if(ptcl2%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport2(ptcl,ptcl2)
           else
              nddmc = nddmc + 1
              call diffusion2(ptcl,ptcl2,icell,specarr)
           endif
!-- verify position
           if(ptcl2%itype==1 .and. .not.ptcl2%done) then
              if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then
                 write(0,*) 'prt_adv: r not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
              endif
              if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
                 write(0,*) 'prt_adv: z not in cell',iy,y,grd_yarr(iy),grd_yarr(iy+1),mu
              endif
           endif
!-- Russian roulette for termination of exhausted particles
           if(e<1d-6*e0 .and. .not.ptcl2%isvacant .and. &
                 grd_capgrey(ic)+grd_sig(ic)>0d0) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl2%itype==1) then
                 labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                      cos(om)*x)/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 ptcl2%isvacant = .true.
                 ptcl2%done = .true.
                 grd_edep(ic) = grd_edep(ic) + e*labfact
!-- velocity effects accounting
                 if(ptcl2%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo !}}}

!-- 3D
     case(3)
        ptcl2%istep = 0!{{{
        do while ((.not.ptcl2%done).and.(.not.ptcl2%isvacant))
           ptcl2%istep = ptcl2%istep + 1
           if(ptcl2%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport3(ptcl,ptcl2)
           else
              nddmc = nddmc + 1
              call diffusion3(ptcl,ptcl2,icell,specarr)
           endif
!-- verify position
           if(ptcl2%itype==1 .and. .not.ptcl2%done) then
              if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then
                 write(0,*) 'prt_adv: x not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu,om
              endif
              if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
                 write(0,*) 'prt_adv: y not in cell',iy,y,grd_yarr(iy),grd_yarr(iy+1),mu,om
              endif
              if(z>grd_zarr(iz+1) .or. z<grd_zarr(iz)) then
                 write(0,*) 'prt_adv: z not in cell',iz,z,grd_zarr(iz),grd_zarr(iz+1),mu,om
              endif
           endif
!-- Russian roulette for termination of exhausted particles
           if(e<1d-6*e0 .and. .not.ptcl2%isvacant .and. &
                 grd_capgrey(ic)+grd_sig(ic)>0d0) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl2%itype==1) then
                 labfact = 1d0-(mu*z+sqrt(1d0-mu**2) * &
                      (cos(om)*x+sin(om)*y))/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 ptcl2%isvacant = .true.
                 ptcl2%done = .true.
                 grd_edep(ic) = grd_edep(ic) + e*labfact
!-- velocity effects accounting
                 if(ptcl2%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo!}}}

!-- 1D
     case(11)
        ptcl2%istep = 0!{{{
        do while ((.not.ptcl2%done).and.(.not.ptcl2%isvacant))
           ptcl2%istep = ptcl2%istep + 1
           if(ptcl2%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport11(ptcl,ptcl2)
           else
              nddmc = nddmc + 1
              call diffusion11(ptcl,ptcl2,icell,specarr)
           endif
!-- verify position
           if(ptcl2%itype==1 .and. .not.ptcl2%done .and. &
                  (x>grd_xarr(ix+1) .or. x<grd_xarr(ix))) then
              write(0,*) 'prt_adv: not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
           endif
!-- Russian roulette for termination of exhausted particles
           if(e<1d-6*e0 .and. .not.ptcl2%isvacant .and. &
                 grd_capgrey(ic)+grd_sig(ic)>0d0) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl2%itype==1) then
                 labfact = 1d0 - mu*x/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 ptcl2%isvacant = .true.
                 ptcl2%done = .true.
                 grd_edep(ic) = grd_edep(ic) + e*labfact
!-- velocity effects accounting
                 if(ptcl2%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo!}}}
     endselect

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
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        x1 = grd_cap(ig,ic)
        x2 = grp_wl(ig)/(pc_c*tsp_t*(grp_wl(ig+1)-grp_wl(ig)))
        if(r1<x2/(x1+x2)) then
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
           wl = wl*exp(tsp_dt/tsp_t)
        endif
        !!}}}
     endif

     if((grd_isvelocity).and.(ptcl2%itype==1)) then
        call advection(.false.,ptcl,ptcl2) !procedure pointer to advection[123]
     endif

!-- radiation energy at census
     tot_erad = tot_erad + e


!-- sample position, direction, and wavelength of a censused DDMC particle
!-- It can be an IMC particle in the next time step
     if(ptcl2%itype==2) then
!
!-- sample wavelength
        r1 = rnd_r(rnd_state)
        wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
!
!-- sample position and direction
        select case(in_igeom)
!-- 3D spherical
        case(1)
!-- sampling position uniformly!{{{
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           y = min(y,grd_yarr(iy+1))
           y = max(y,grd_yarr(iy))
           z = min(z,grd_zarr(iz+1))
           z = max(z,grd_zarr(iz))
!-- sampling angle isotropically
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           mu = 1.0 - 2.0*r1
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           om = pc_pi2*r1!}}}
!-- 2D
        case(2)
!-- sampling position uniformly!{{{
           r1 = rnd_r(rnd_state)
           x = sqrt(r1*grd_xarr(ix+1)**2 + (1d0-r1)*grd_xarr(ix)**2)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           r1 = rnd_r(rnd_state)
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
!-- sampling direction values
           r1 = rnd_r(rnd_state)
           om = pc_pi2*r1
           r1 = rnd_r(rnd_state)
           mu = 1d0 - 2d0*r1!}}}
!-- 3D
        case(3)
!-- sampling position uniformly !{{{
           r1 = rnd_r(rnd_state)
           x = r1*grd_xarr(ix+1)+(1d0-r1)*grd_xarr(ix)
           r1 = rnd_r(rnd_state)
           y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
           r1 = rnd_r(rnd_state)
           z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
           y = min(y,grd_yarr(iy+1))
           y = max(y,grd_yarr(iy))
           z = min(z,grd_zarr(iz+1))
           z = max(z,grd_zarr(iz))
!-- sampling direction values
           r1 = rnd_r(rnd_state)
           om = pc_pi2*r1
           r1 = rnd_r(rnd_state)
           mu = 1d0 - 2d0*r1 !}}}
!-- 1D spherical
        case(11)
!-- sampling position uniformly!{{{
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           x = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
           x = min(x,grd_xarr(ix+1))
           x = max(x,grd_xarr(ix))
!-- sampling angle isotropically
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           mu = 1.0 - 2.0*r1 !}}}
        endselect
!
        if(grd_isvelocity) call direction2lab(x,y,z,mu,om)
     endif

!
!-- transform IMC particle energy to comoving frame for storage
     if(grd_isvelocity.and.ptcl2%itype==1) then
        select case(in_igeom)
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

  t1 = t_time()
  t_pckt_stat = t1-t0  !register timing
  call timereg(t_pcktnpckt, dble(npckt))
  call timereg(t_pcktnddmc, dble(nddmc))
  call timereg(t_pcktnimc, dble(nimc))


end subroutine particle_advance
