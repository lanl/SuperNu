subroutine particle_advance

  use particlemod
  use timestepmod
  use gasgridmod
  use physconstmod
  use inputparmod
  use timingmod
  use mpimod
  implicit none

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
  integer,external :: binsrch
  real*8 :: r1, x1, x2, help, elabfact
! integer :: irl,irr
! real*8 :: xx0, bmax
! real*8 :: uul, uur, uumax, r0,r2,r3
  logical,pointer :: isvacant
  integer, pointer :: zsrc, iy, iz
  real*8, pointer :: rsrc, musrc, esrc, wlsrc, y, z, om
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, azitrfm
!
  type(packet),pointer :: ptcl
!
  logical,parameter :: isshift=.true.
!-- statement function
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)

  gas_edep = 0.0
  gas_erad = 0.0
  gas_luminos = 0.0
  gas_lumdev = 0.0
  gas_lumnum = 0
  gas_methodswap = 0
!
!--(rev. 121)
  gas_eraddens =0d0
!--
  gas_numcensus = 0
  
  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  npckt = 0
  nddmc = 0
  nimc = 0
  do ipart = 1, prt_npartmax
     ! Checking vacancy
     if(prt_isvacant(ipart)) cycle
!
!-- active particle
     isvacant => prt_isvacant(ipart)
     ptcl => prt_particles(ipart)
     npckt = npckt + 1
!
!-- assigning pointers to corresponding particle properties
     select case(in_igeom)
!-- 1D
     case(1)
        zsrc => ptcl%zsrc
        wlsrc => ptcl%wlsrc
        rsrc => ptcl%rsrc
        musrc => ptcl%musrc
        esrc => ptcl%esrc
!-- 1-dir*v/c
        if(gas_isvelocity.and.ptcl%rtsrc==1) then
           labfact = 1d0-rsrc*musrc/pc_c
        endif

!-- 2D
     case(2)
        zsrc => ptcl%zsrc
        iy => ptcl%iy
        wlsrc => ptcl%wlsrc
        rsrc => ptcl%rsrc
        y => ptcl%y
        musrc => ptcl%musrc
        om => ptcl%om
        esrc => ptcl%esrc
!-- 1-dir*v/c
        if(gas_isvelocity.and.ptcl%rtsrc==1) then
           labfact = 1d0-(musrc*y+sqrt(1d0-musrc**2) * &
                cos(om)*rsrc)/pc_c
        endif

!-- 3D
     case(3)
        stop 'particle_advance: no 3D transport'
     endselect

     prt_done=.false.

     ! Looking up group
     if(ptcl%rtsrc==1) then
        if(gas_isvelocity) then!{{{
           ig = binsrch(wlsrc/labfact,gas_wl,gas_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)
        endif
        !
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>gas_ng) then
              ig=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*labfact
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(gas_isvelocity) then
                 wlsrc=gas_wl(1)*labfact
              else
                 wlsrc=gas_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)!{{{
        !
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>gas_ng) then
              ig=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=gas_wl(1)
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
        if(gas_isvelocity) then!{{{
           help = tsp_t
        else
           help = 1d0
        endif
!
!-- selecting geometry
        select case(in_igeom)

!-- 1D
        case(1)
           lhelp = (gas_sig(zsrc,1,1)+gas_cap(ig,zsrc,1,1)) * &
                dx(zsrc)*help<prt_tauddmc
           if (lhelp) then
              if (ptcl%rtsrc == 2) then
                 gas_methodswap(zsrc,1,1)=gas_methodswap(zsrc,1,1)+1
!-- sampling position uniformly
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 rsrc = (r1*gas_xarr(zsrc+1)**3 + &
                      (1.0-r1)*gas_xarr(zsrc)**3)**(1.0/3.0)
!-- sampling angle isotropically
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 musrc = 1.0 - 2.0*r1
                 if(gas_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+rsrc*musrc/pc_c
!-- mu
                    musrc = (musrc+rsrc/pc_c)/cmffact
                 endif
!-- 1-dir*v/c
                 labfact = 1d0-rsrc*musrc/pc_c
              else
                 if(ptcl%rtsrc==1) then
                    gas_methodswap(zsrc,1,1)=gas_methodswap(zsrc,1,1)+1
                 endif
              endif!}}}
           endif

!-- 2D
        case(2)

!-- 3D
        case(3)
           stop 'particle_advance: no 3D transport'
        endselect

        if (lhelp) then
           if (ptcl%rtsrc == 2) then
              if(gas_isvelocity) then
!-- velocity effects accounting
                 gas_evelo=gas_evelo+esrc*(1d0-1d0/labfact)
!
                 esrc = esrc/labfact
                 ptcl%ebirth = ptcl%ebirth/labfact
              endif
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wlsrc = 1d0/(r1/gas_wl(ig+1)+(1d0-r1)/gas_wl(ig))
              if(gas_isvelocity) then
                 wlsrc = wlsrc*labfact
              endif
           endif
           ptcl%rtsrc = 1
        else
           if(ptcl%rtsrc==1.and.gas_isvelocity) then
              gas_evelo = gas_evelo+esrc*(1d0-labfact)
              esrc = esrc*labfact
              ptcl%ebirth = ptcl%ebirth*labfact
              wlsrc = wlsrc/labfact
           endif
           ptcl%rtsrc = 2
        endif!}}}
     endif 
!
!-- looking up group
     if(ptcl%rtsrc==1) then
        if(gas_isvelocity) then!{{{
           ig = binsrch(wlsrc/(1.0d0-rsrc*musrc/pc_c),gas_wl,gas_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)
        endif
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>gas_ng) then
              ig=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(gas_isvelocity) then
                 wlsrc=gas_wl(1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)!{{{
        !
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>gas_ng) then
              ig=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=gas_wl(1)
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
     if ((gas_isvelocity).and.(ptcl%rtsrc==1)) then
        call advection1(.true.,ig,zsrc,rsrc)
     endif
        !
     endif

!     write(*,*) ipart
!-----------------------------------------------------------------------        
     ! Advancing particle until census, absorption, or escape from domain
     do while ((.not.prt_done).and.(.not.isvacant))
        !Calling either diffusion or transport depending on particle type (ptcl%rtsrc)
        if (ptcl%rtsrc == 1.or.in_puretran) then
           nimc = nimc + 1
           call transport1(ptcl,isvacant)
        else
           nddmc = nddmc + 1
           call diffusion1(ptcl,isvacant)
        endif

!-- transformation factor
        if(gas_isvelocity .and. ptcl%rtsrc==1) then
           elabfact = 1.0d0 - musrc*rsrc/pc_c
        else
           elabfact = 1d0
        endif

!-- Russian roulette for termination of exhausted particles
        if (esrc<1d-6*ptcl%ebirth .and. .not.isvacant) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(r1<0.5d0) then
              isvacant = .true.
              prt_done = .true.
              gas_edep(zsrc,1,1) = gas_edep(zsrc,1,1) + esrc*elabfact
!-- velocity effects accounting
              if(ptcl%rtsrc==1) gas_evelo = gas_evelo + esrc*(1d0-elabfact)
           else
!-- weight addition accounted for in external source
              gas_eext = gas_eext + esrc
!
              esrc = 2d0*esrc
              ptcl%ebirth = 2d0*ptcl%ebirth
           endif
        endif
     enddo
!-----------------------------------------------------------------------


     if(.not.isvacant) then

     ! Redshifting DDMC particle energy weights and wavelengths
     if(ptcl%rtsrc == 2.and.gas_isvelocity) then
!-- redshifting energy weight!{{{
        gas_evelo=gas_evelo+esrc*(1d0-exp(-tsp_dt/tsp_t))
        esrc = esrc*exp(-tsp_dt/tsp_t)
        ptcl%ebirth = ptcl%ebirth*exp(-tsp_dt/tsp_t)
        !
!
!-- find group
        ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)
        !
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>gas_ng) then
              ig=gas_ng
           else
              ig=1
           endif
        endif
        !
        !
        if(ig<gas_ng) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           x1 = gas_cap(ig,zsrc,1,1)
           x2 = gas_wl(ig)/(pc_c*tsp_t*(gas_wl(ig+1)-gas_wl(ig)))
           if(r1<x2/(x1+x2)) then
!            if((gas_sig(zsrc,1,1)+gas_cap(ig+1,zsrc,1,1))*dx(zsrc) &
!                 *tsp_t>=prt_tauddmc) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wlsrc = 1d0/(r1/gas_wl(ig+1)+(1d0-r1)/gas_wl(ig))
              wlsrc = wlsrc*exp(tsp_dt/tsp_t)
!            endif
           endif
        endif
        !!}}}
     endif

     endif

     ! Looking up group
     if(ptcl%rtsrc==1) then
        if(gas_isvelocity) then!{{{
           ig = binsrch(wlsrc/(1.0d0-rsrc*musrc/pc_c),gas_wl,gas_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)
        endif
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>gas_ng) then
              ig=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(gas_isvelocity) then
                 wlsrc=gas_wl(1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,gas_wl,gas_ng+1,in_ng)!{{{
        !
        if(ig>gas_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>gas_ng) then
              ig=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=gas_wl(1)
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     endif
     
     if(isshift) then
     if ((gas_isvelocity).and.(ptcl%rtsrc==1)) then
       call advection1(.false.,ig,zsrc,rsrc)
     endif
     endif

     if(.not.isvacant) then
!
!-- radiation energy at census
     if(gas_isvelocity) then
        if(ptcl%rtsrc==2) then
           gas_erad = gas_erad + esrc
        else
           gas_erad = gas_erad + esrc !*(1d0-musrc*rsrc/pc_c)
!-- velocity effects accounting
!           gas_evelo=gas_evelo+esrc*musrc*rsrc/pc_c
!
        endif
     else
        gas_erad = gas_erad + esrc
     endif

     endif

  enddo !ipart

  call time(t1)
  t_pckt_stat = t1-t0  !register timing
  call timereg(t_pckt, t1-t0)
  call timereg(t_pcktnpckt, dble(npckt))
  call timereg(t_pcktnddmc, dble(nddmc))
  call timereg(t_pcktnimc, dble(nimc))
  npckt = 0
  nddmc = 0
  nimc = 0
  !write(6,*) eleft, eright

  gas_eext = gas_eext-gas_eleft-gas_eright

end subroutine particle_advance


              !wlsrc = 0.5d0*(gas_wl(ig)+gas_wl(ig+1))
              !r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
              !wlsrc=gas_wl(ig)*(1d0-r1)+gas_wl(ig+1)*r1
              !
!               r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!               if(r1<gas_cap(ig,zsrc,1,1)/(gas_cap(ig,zsrc,1,1)+gas_sig(zsrc,1,1))) then
!                  x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_temp(zsrc,1,1))
!                  x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(zsrc,1,1))
!                  if (x2<pc_plkpk) then
!                     bmax = x2**3/(exp(x2)-1d0)
!                  elseif (x1>pc_plkpk) then
!                     bmax = x1**3/(exp(x1)-1d0)
!                  else
!                     bmax = pc_plkpk
!                  endif
!                  r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!                  r2 = rand()
!           prt_tlyrand = prt_tlyrand+1
!                  xx0 = (1d0-r1)*x1+r1*x2
!                  do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!                     r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!                     r2 = rand()
!           prt_tlyrand = prt_tlyrand+1
!                     xx0 = (1d0-r1)*x1+r1*x2
!                  enddo
!                  wlsrc = pc_h*pc_c/(xx0*pc_kb*gas_temp(zsrc,1,1))
!               else
