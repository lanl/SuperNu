subroutine particle_advance

  use randommod
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
  logical :: lhelp
  integer*8 :: nddmc, nimc, npckt
  integer :: ipart, ig
  real*8 :: r1, x1, x2, help
! integer :: irl,irr
! real*8 :: xx0, bmax
! real*8 :: uul, uur, uumax, r0,r2,r3
  logical :: isvacant
  integer, pointer :: ix, iy, iz
  real*8, pointer :: x,y,z, mu, e, wl, om
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, mu1, mu2, gm
!-- specint cache
  real*8 :: specarr(grp_ng)
  integer :: icell(3)
!
  type(packet),target :: ptcl
!
  logical,parameter :: isshift=.true.
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
  ix => ptcl%ix
  iy => ptcl%iy
  iz => ptcl%iz
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  wl => ptcl%wl
  e => ptcl%e

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

  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  npckt = 0
  nddmc = 0
  nimc = 0
  icell = 0
  do ipart = 1, prt_npartmax
     ! Checking vacancy
     if(prt_isvacant(ipart)) cycle
     isvacant = .false.
!
!-- active particle
     ptcl = prt_particles(ipart)
     prt_ipart = ipart
     npckt = npckt + 1

!-- default, recalculated for isvelocity and itype==1
     labfact = 1d0

     if(grd_isvelocity.and.ptcl%itype==1) then
        select case(in_igeom)
!-- [123]D spherical
        case(1,4)
           labfact = 1d0-x*mu/pc_c !-- 1-dir*v/c
!-- 2D
        case(2)
           labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * cos(om)*x)/pc_c !-- 1-dir*v/c
!-- 3D
        case(3)
           mu1 = sqrt(1d0-mu**2)*cos(om)
           mu2 = sqrt(1d0-mu**2)*sin(om)
           labfact = 1d0-(mu*z+mu1*x+mu2*y)/pc_c !-- 1-dir*v/c
        endselect
     endif

     prt_done=.false.

     ! Looking up group
     if(ptcl%itype==1) then
        if(grd_isvelocity) then!{{{
           ig = binsrch(wl/labfact,grp_wl,grp_ng+1)
        else
           ig = binsrch(wl,grp_wl,grp_ng+1)
        endif
        !
        if(ig>grp_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grp_ng) then
              ig=grp_ng
              if(grd_isvelocity) then
                 wl=grp_wl(grp_ng+1)*labfact
              else
                 wl=grp_wl(grp_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(grd_isvelocity) then
                 wl=grp_wl(1)*labfact
              else
                 wl=grp_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wl,grp_wl,grp_ng+1)!{{{
        !
        if(ig>grp_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>grp_ng) then
              ig=grp_ng
              wl=grp_wl(grp_ng+1)
           elseif(ig<1) then
              ig=1
              wl=grp_wl(1)
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     endif


     ! Checking if particle conversions are required since prior time step
     if(.not.in_puretran) then
        if(grd_isvelocity) then!{{{
           help = tsp_t
        else
           help = 1d0
        endif
!
!-- selecting geometry
        select case(in_igeom)

!-- 3D spherical
        case(1)
           lhelp = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz)) * &!{{{
                min(dx(ix),xm(ix)*dyac(iy),xm(ix)*ym(iy)*dz(iz)) * &
                help<prt_tauddmc).or.in_puretran
           if (lhelp) then
              if (ptcl%itype == 2) then
!-- DDMC -> IMC
                 grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
!-- sampling position uniformly
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 x = (r1*grd_xarr(ix+1)**3 + &
                      (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- must be inside cell
                 x = min(x,grd_xarr(ix+1))
                 x = max(x,grd_xarr(ix))
!-- sampling angle isotropically
                 r1 = rnd_r(rnd_state)
                 prt_tlyrand = prt_tlyrand+1
                 mu = 1.0 - 2.0*r1
                 r1 = rnd_r(rnd_state)
                 prt_tlyrand = prt_tlyrand+1
                 om = pc_pi2*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+x*mu/pc_c
!-- mu
                    mu = (mu+x/pc_c)/cmffact
!-- 1-dir*v/c
                    labfact = 1d0-x*mu/pc_c
                 endif
              endif
           else
              if(ptcl%itype==1) then
!-- IMC -> DDMC
                 grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
              endif
           endif!}}}

!-- 2D
        case(2)
           lhelp = ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1)) * & !{{{
                min(dx(ix),dy(iy))*help < prt_tauddmc) &
                .or.in_puretran
           if (lhelp) then
              if (ptcl%itype == 2) then
!-- DDMC -> IMC
                 grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- sampling position uniformly
                 r1 =  rand()
                 x = sqrt(r1*grd_xarr(ix+1)**2 + &
                      (1d0-r1)*grd_xarr(ix)**2)
!-- must be inside cell
                 x = min(x,grd_xarr(ix+1))
                 x = max(x,grd_xarr(ix))
                 r1 = rnd_r(rnd_state)
                 y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
!-- sampling direction values
                 r1 = rnd_r(rnd_state)
                 om = pc_pi2*r1
                 r1 = rnd_r(rnd_state)
                 mu = 1d0 - 2d0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+(mu*y+sqrt(1d0-mu**2) * &
                         cos(om)*x)/pc_c
                    gm = 1d0/sqrt(1d0-(x**2+y**2)/pc_c**2)
!-- om
                    om = atan2(sqrt(1d0-mu**2)*sin(om) , &
                         sqrt(1d0-mu**2)*cos(om)+(gm*x/pc_c) * &
                         (1d0+gm*(cmffact-1d0)/(gm+1d0)))
                    if(om<0d0) om=om+pc_pi2
!-- mu
                    mu = (mu+(gm*y/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                         (gm*cmffact)
!-- 1-dir*v/c
                    labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                         cos(om)*x)/pc_c
                 endif
              endif
           else
              if(ptcl%itype==1) then
!-- IMC -> DDMC
                 grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
              endif
           endif!}}}

!-- 3D
        case(3)
           lhelp = ((grd_sig(ix,iy,iz)+grd_cap(ig,ix,iy,iz)) * & !{{{
                min(dx(ix),dy(iy),dz(iz))*help < prt_tauddmc) &
                .or.in_puretran
           if (lhelp) then
              if (ptcl%itype == 2) then
!-- DDMC -> IMC
                 grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
!-- sampling position uniformly
                 r1 =  rand()
                 x = r1*grd_xarr(ix+1)+(1d0-r1)*grd_xarr(ix)
                 r1 = rnd_r(rnd_state)
                 y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
                 r1 = rnd_r(rnd_state)
                 z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- sampling direction values
                 r1 = rnd_r(rnd_state)
                 om = pc_pi2*r1
                 r1 = rnd_r(rnd_state)
                 mu = 1d0 - 2d0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    mu1 = sqrt(1d0-mu**2)*cos(om)
                    mu2 = sqrt(1d0-mu**2)*sin(om)
                    cmffact = 1d0+(mu*z+mu1*x+mu2*y)/pc_c
!-- mu
                    mu = (mu+z/pc_c)/cmffact
                    if(mu>1d0) then
                       mu = 1d0
                    elseif(mu<-1d0) then
                       mu = -1d0
                    endif
!-- om
                    om = atan2(mu2+y/pc_c,mu1+x/pc_c)
                    if(om<0d0) om = om+pc_pi2
!-- 1-dir*v/c
                    mu1 = sqrt(1d0-mu**2)*cos(om)
                    mu2 = sqrt(1d0-mu**2)*sin(om)
                    labfact = 1d0-(mu*z+mu1*x+mu2*y)/pc_c
                 endif
              endif
           else
              if(ptcl%itype==1) then
!-- IMC -> DDMC
                 grd_methodswap(ix,iy,iz)=grd_methodswap(ix,iy,iz)+1
              endif
           endif!}}}

!-- 1D spherical
        case(4)
           lhelp = ((grd_sig(ix,1,1)+grd_cap(ig,ix,1,1)) * &!{{{
                dx(ix)*help<prt_tauddmc) &
                .or.in_puretran
           if (lhelp) then
              if (ptcl%itype == 2) then
!-- DDMC -> IMC
                 grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
!-- sampling position uniformly
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 x = (r1*grd_xarr(ix+1)**3 + &
                      (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
                 x = min(x,grd_xarr(ix+1))
                 x = max(x,grd_xarr(ix))
!-- sampling angle isotropically
                 r1 = rnd_r(rnd_state)
                 prt_tlyrand = prt_tlyrand+1
                 mu = 1.0 - 2.0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+x*mu/pc_c
!-- mu
                    mu = (mu+x/pc_c)/cmffact
!-- 1-dir*v/c
                    labfact = 1d0-x*mu/pc_c
                 endif
              endif
           else
              if(ptcl%itype==1) then
!-- IMC -> DDMC
                 grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
              endif
           endif!}}}

        case default
!-- don't let the compiler believe lhelp may be used uninitialized
           lhelp = .true.
        endselect


        if (lhelp) then
           if (ptcl%itype == 2) then
!-- DDMC -> IMC
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
              if(grd_isvelocity) then
!-- velocity effects accounting
                 tot_evelo=tot_evelo+e*(1d0-1d0/labfact)
!
                 e = e/labfact
                 ptcl%e0 = ptcl%e0/labfact
                 wl = wl*labfact
              endif
              ptcl%itype = 1
           endif
        else
           if(ptcl%itype==1) then
!-- IMC -> DDMC
              if(grd_isvelocity) then
!-- velocity effects accounting
                 tot_evelo = tot_evelo+e*(1d0-labfact)
!
                 e = e*labfact
                 ptcl%e0 = ptcl%e0*labfact
                 wl = wl/labfact
              endif
              ptcl%itype = 2
           endif
        endif!}}}
     endif
!
!-- looking up group
     if(ptcl%itype==1) then
        if(grd_isvelocity) then!{{{
           ig = binsrch(wl/labfact,grp_wl,grp_ng+1)
        else
           ig = binsrch(wl,grp_wl,grp_ng+1)
        endif
        if(ig>grp_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grp_ng) then
              ig=grp_ng
              if(grd_isvelocity) then
                 wl=grp_wl(grp_ng+1)*labfact
              else
                 wl=grp_wl(grp_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(grd_isvelocity) then
                 wl=grp_wl(1)*labfact
              else
                 wl=grp_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wl,grp_wl,grp_ng+1)!{{{
        !
        if(ig>grp_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>grp_ng) then
              ig=grp_ng
              wl=grp_wl(grp_ng+1)
           elseif(ig<1) then
              ig=1
              wl=grp_wl(1)
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     endif

!-- First portion of operator split particle velocity position adjustment
     if(isshift) then
     if ((grd_isvelocity).and.(ptcl%itype==1)) then
        select case(in_igeom)
!-- [123]D spherical
        case(1,4)
           call advection1(.true.,ptcl,ig)
!-- 2D
        case(2)
           call advection2(.true.,ptcl,ig)
!-- 3D
        case(3)
           call advection3(.true.,ptcl,ig)
        endselect
     endif
     endif

!     write(*,*) ipart
!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
!Calling either diffusion or transport depending on particle type (ptcl%itype)
     select case(in_igeom)

!-- 1D
     case(1,4)
        prt_istep = 0
        do while ((.not.prt_done).and.(.not.isvacant))
           prt_istep = prt_istep + 1
           if (ptcl%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport11(ptcl,ig,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion11(ptcl,ig,isvacant,icell,specarr)
           endif
!-- verify position
           if(ptcl%itype==1 .and. .not.prt_done .and. &
                  (x>grd_xarr(ix+1) .or. x<grd_xarr(ix))) then
              write(0,*) 'prt_adv: not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
           endif
!-- Russian roulette for termination of exhausted particles
           if (e<1d-6*ptcl%e0 .and. .not.isvacant) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl%itype==1) then
                 labfact = 1d0 - mu*x/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
!-- velocity effects accounting
                 if(ptcl%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 ptcl%e0 = 2d0*ptcl%e0
              endif
           endif
        enddo

!-- 2D
     case(2)
        prt_istep = 0
        do while ((.not.prt_done).and.(.not.isvacant))
           prt_istep = prt_istep + 1
           if (ptcl%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport2(ptcl,ig,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion2(ptcl,ig,isvacant,icell,specarr)
           endif
!-- verify position
           if(ptcl%itype==1 .and. .not.prt_done) then
              if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then
                 write(0,*) 'prt_adv: r not in cell',ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
              endif
              if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
                 write(0,*) 'prt_adv: z not in cell',iy,y,grd_yarr(iy),grd_yarr(iy+1),mu
              endif
           endif
!-- Russian roulette for termination of exhausted particles
           if (e<1d-6*ptcl%e0 .and. .not.isvacant) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl%itype==1) then
                 labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                      cos(om)*x)/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
!-- velocity effects accounting
                 if(ptcl%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 ptcl%e0 = 2d0*ptcl%e0
              endif
           endif
        enddo

!-- 3D
     case(3)
        prt_istep = 0
        do while ((.not.prt_done).and.(.not.isvacant))
           prt_istep = prt_istep + 1
           if (ptcl%itype == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport3(ptcl,ig,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion3(ptcl,ig,isvacant,icell,specarr)
           endif
!-- verify position
           if(ptcl%itype==1 .and. .not.prt_done) then
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
           if (e<1d-6*ptcl%e0 .and. .not.isvacant) then
!-- transformation factor
              if(grd_isvelocity .and. ptcl%itype==1) then
                 labfact = 1d0-(mu*z+sqrt(1d0-mu**2) * &
                      (cos(om)*x+sin(om)*y))/pc_c
              else
                 labfact = 1d0
              endif
!
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
!-- velocity effects accounting
                 if(ptcl%itype==1) tot_evelo = tot_evelo + e*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + e
!
                 e = 2d0*e
                 ptcl%e0 = 2d0*ptcl%e0
              endif
           endif
        enddo

     endselect

!-----------------------------------------------------------------------


!-- Redshifting DDMC particle energy weights and wavelengths
     if(.not.isvacant) then
     if(ptcl%itype == 2.and.grd_isvelocity) then
!-- redshifting energy weight!{{{
        tot_evelo=tot_evelo+e*(1d0-exp(-tsp_dt/tsp_t))
        e = e*exp(-tsp_dt/tsp_t)
        ptcl%e0 = ptcl%e0*exp(-tsp_dt/tsp_t)
        !
!
!-- find group
        ig = binsrch(wl,grp_wl,grp_ng+1)
        !
        if(ig>grp_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grp_ng) then
              ig=grp_ng
           else
              ig=1
           endif
        endif
        !
        !
        if(ig<grp_ng) then
           r1 = rnd_r(rnd_state)
           prt_tlyrand = prt_tlyrand+1
           x1 = grd_cap(ig,ix,iy,iz)
           x2 = grp_wl(ig)/(pc_c*tsp_t*(grp_wl(ig+1)-grp_wl(ig)))
           if(r1<x2/(x1+x2)) then
              r1 = rnd_r(rnd_state)
              prt_tlyrand = prt_tlyrand+1
              wl = 1d0/(r1*grp_wlinv(ig+1)+(1d0-r1)*grp_wlinv(ig))
              wl = wl*exp(tsp_dt/tsp_t)
           endif
        endif
        !!}}}
     endif
     endif

     if(isshift) then
     if ((grd_isvelocity).and.(ptcl%itype==1)) then
        select case(in_igeom)
!-- [123]D
        case(1,4)
           call advection1(.false.,ptcl,ig)
!-- 2D
        case(2)
           call advection2(.false.,ptcl,ig)
!-- 3D
        case(3)
           call advection3(.false.,ptcl,ig)
        endselect
     endif
     endif

!-- radiation energy at census
     if(.not.isvacant) tot_erad = tot_erad + e

!
!-- save particle results
!------------------------
     prt_isvacant(ipart) = isvacant
     if(.not.isvacant) prt_particles(ipart) = ptcl

  enddo !ipart

  call time(t1)
  t_pckt_stat = t1-t0  !register timing
  call timereg(t_pcktnpckt, dble(npckt))
  call timereg(t_pcktnddmc, dble(nddmc))
  call timereg(t_pcktnimc, dble(nimc))


end subroutine particle_advance
