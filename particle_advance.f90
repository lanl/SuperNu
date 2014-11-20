subroutine particle_advance

  use particlemod
  use timestepmod
  use totalsmod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use mpimod
  use fluxmod
  implicit none
!
  integer,target :: one = 1
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
  real*8 :: r1, x1, x2, help
! integer :: irl,irr
! real*8 :: xx0, bmax
! real*8 :: uul, uur, uumax, r0,r2,r3
  logical,pointer :: isvacant
  integer, pointer :: zsrc, iy, iz
  real*8, pointer :: rsrc,y,z, musrc, esrc, wlsrc, om
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, azitrfm, mu1, mu2
!
  type(packet),pointer :: ptcl
!
  logical,parameter :: isshift=.true.
!-- statement function
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)

  grd_edep = 0.0
  tot_erad = 0.0
  flx_luminos = 0.0
  flx_lumdev = 0.0
  flx_lumnum = 0
  grd_methodswap = 0
!
!--(rev. 121)
  grd_eraddens =0d0
!--
  grd_numcensus = 0
  
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
     zsrc => ptcl%zsrc
     wlsrc => ptcl%wlsrc
     rsrc => ptcl%rsrc
     musrc => ptcl%musrc
     esrc => ptcl%esrc
     select case(in_igeom)

!-- 1D
     case(1)
        iy => one
        iz => one
!-- 1-dir*v/c
        if(grd_isvelocity.and.ptcl%rtsrc==1) then
           labfact = 1d0-rsrc*musrc/pc_c
        endif

!-- 2D
     case(2)
        iy => ptcl%iy
        iz => one
        y => ptcl%y
        om => ptcl%om
!-- 1-dir*v/c
        if(grd_isvelocity.and.ptcl%rtsrc==1) then
           labfact = 1d0-(musrc*y+sqrt(1d0-musrc**2) * &
                cos(om)*rsrc)/pc_c
        endif

!-- 3D
     case(3)
        iy => ptcl%iy
        iz => ptcl%iz
        y => ptcl%y
        z => ptcl%z
        om => ptcl%om
!-- 1-dir*v/c
        if(grd_isvelocity.and.ptcl%rtsrc==1) then
           mu1 = sqrt(1d0-musrc**2)*cos(om)
           mu2 = sqrt(1d0-musrc**2)*sin(om)
           labfact = 1d0-(musrc*z+mu1*rsrc+mu2*y)/pc_c
        endif
     endselect

     prt_done=.false.

     ! Looking up group
     if(ptcl%rtsrc==1) then
        if(grd_isvelocity) then!{{{
           ig = binsrch(wlsrc/labfact,grd_wl,grd_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)
        endif
        !
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grd_ng) then
              ig=grd_ng
              if(grd_isvelocity) then
                 wlsrc=grd_wl(grd_ng+1)*labfact
              else
                 wlsrc=grd_wl(grd_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(grd_isvelocity) then
                 wlsrc=grd_wl(1)*labfact
              else
                 wlsrc=grd_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)!{{{
        !
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>grd_ng) then
              ig=grd_ng
              wlsrc=grd_wl(grd_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=grd_wl(1)
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

!-- 1D
        case(1)
           lhelp = (grd_sig(zsrc,1,1)+grd_cap(ig,zsrc,1,1)) * &!{{{
                dx(zsrc)*help<prt_tauddmc
           if (lhelp) then
              if (ptcl%rtsrc == 2) then
!-- DDMC -> IMC
                 grd_methodswap(zsrc,1,1)=grd_methodswap(zsrc,1,1)+1
!-- sampling position uniformly
                 r1 =  rand()
                 prt_tlyrand = prt_tlyrand+1
                 rsrc = (r1*grd_xarr(zsrc+1)**3 + &
                      (1.0-r1)*grd_xarr(zsrc)**3)**(1.0/3.0)
!-- sampling angle isotropically
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 musrc = 1.0 - 2.0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+rsrc*musrc/pc_c
!-- mu
                    musrc = (musrc+rsrc/pc_c)/cmffact
!-- 1-dir*v/c
                    labfact = 1d0-rsrc*musrc/pc_c
                 endif
              endif
           else
              if(ptcl%rtsrc==1) then
!-- IMC -> DDMC
                 grd_methodswap(zsrc,1,1)=grd_methodswap(zsrc,1,1)+1
              endif
           endif!}}}

!-- 2D
        case(2)
           lhelp = ((grd_sig(zsrc,iy,1)+grd_cap(ig,zsrc,iy,1)) * & !{{{
                min(dx(zsrc),dy(iy))*help < prt_tauddmc) &
                .or.in_puretran
           if (lhelp) then
              if (ptcl%rtsrc == 2) then
!-- DDMC -> IMC
                 grd_methodswap(zsrc,iy,1)=grd_methodswap(zsrc,iy,1)+1
!-- sampling position uniformly
                 r1 =  rand()
                 rsrc = sqrt(r1*grd_xarr(zsrc+1)**2 + &
                      (1d0-r1)*grd_xarr(zsrc)**2)
                 r1 = rand()
                 y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
!-- sampling direction values
                 r1 = rand()
                 om = pc_pi2*r1
                 r1 = rand()
                 musrc = 1d0 - 2d0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    cmffact = 1d0+(musrc*y+sqrt(1d0-musrc**2) * &
                         cos(om)*rsrc)/pc_c
                    azitrfm = atan2(sqrt(1d0-musrc**2)*sin(om) , &
                         sqrt(1d0-musrc**2)*cos(om)+rsrc/pc_c)
!-- mu
                    musrc = (musrc+y/pc_c)/cmffact
                    if(musrc>1d0) then
                       musrc = 1d0
                    elseif(musrc<-1d0) then
                       musrc = -1d0
                    endif
!-- om
                    if(azitrfm >= 0d0) then
                       om = azitrfm
                    else
                       om = azitrfm+pc_pi2
                    endif
!-- 1-dir*v/c
                    labfact = 1d0-(musrc*y+sqrt(1d0-musrc**2) * &
                         cos(om)*rsrc)/pc_c
                 endif
              endif
           else
              if(ptcl%rtsrc==1) then
!-- IMC -> DDMC
                 grd_methodswap(zsrc,iy,1)=grd_methodswap(zsrc,iy,1)+1
              endif
           endif!}}}

!-- 3D
        case(3)
           lhelp = ((grd_sig(zsrc,iy,iz)+grd_cap(ig,zsrc,iy,iz)) * & !{{{
                min(dx(zsrc),dy(iy),dz(iz))*help < prt_tauddmc) &
                .or.in_puretran
           if (lhelp) then
              if (ptcl%rtsrc == 2) then
!-- DDMC -> IMC
                 grd_methodswap(zsrc,iy,iz)=grd_methodswap(zsrc,iy,iz)+1
!-- sampling position uniformly
                 r1 =  rand()
                 rsrc = r1*grd_xarr(zsrc+1)+(1d0-r1)*grd_xarr(zsrc)
                 r1 = rand()
                 y = r1*grd_yarr(iy+1)+(1d0-r1)*grd_yarr(iy)
                 r1 = rand()
                 z = r1*grd_zarr(iz+1)+(1d0-r1)*grd_zarr(iz)
!-- sampling direction values
                 r1 = rand()
                 om = pc_pi2*r1
                 r1 = rand()
                 musrc = 1d0 - 2d0*r1
                 if(grd_isvelocity) then
!-- 1+dir*v/c
                    mu1 = sqrt(1d0-musrc**2)*cos(om)
                    mu2 = sqrt(1d0-musrc**2)*sin(om)
                    cmffact = 1d0+(musrc*z+mu1*rsrc+mu2*y)/pc_c
!-- mu
                    musrc = (musrc+z/pc_c)/cmffact
                    if(musrc>1d0) then
                       musrc = 1d0
                    elseif(musrc<-1d0) then
                       musrc = -1d0
                    endif
!-- om
                    om = atan2(mu2+y/pc_c,mu1+rsrc/pc_c)
                    if(om<0d0) om = om+pc_pi2
!-- 1-dir*v/c
                    labfact = 1d0-(musrc*z+mu1*rsrc+mu2*y)/pc_c
                 endif
              endif
           else
              if(ptcl%rtsrc==1) then
!-- IMC -> DDMC
                 grd_methodswap(zsrc,iy,iz)=grd_methodswap(zsrc,iy,iz)+1
              endif
           endif!}}}

        endselect

        if (lhelp) then
           if (ptcl%rtsrc == 2) then
!-- DDMC -> IMC
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wlsrc = 1d0/(r1/grd_wl(ig+1)+(1d0-r1)/grd_wl(ig))
              if(grd_isvelocity) then
!-- velocity effects accounting
                 tot_evelo=tot_evelo+esrc*(1d0-1d0/labfact)
!
                 esrc = esrc/labfact
                 ptcl%ebirth = ptcl%ebirth/labfact
                 wlsrc = wlsrc*labfact
              endif
           endif
           ptcl%rtsrc = 1
        else
           if(ptcl%rtsrc==1.and.grd_isvelocity) then
!-- IMC -> DDMC
              tot_evelo = tot_evelo+esrc*(1d0-labfact)
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
        if(grd_isvelocity) then!{{{
           ig = binsrch(wlsrc/labfact,grd_wl,grd_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)
        endif
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grd_ng) then
              ig=grd_ng
              if(grd_isvelocity) then
                 wlsrc=grd_wl(grd_ng+1)*labfact
              else
                 wlsrc=grd_wl(grd_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(grd_isvelocity) then
                 wlsrc=grd_wl(1)*labfact
              else
                 wlsrc=grd_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)!{{{
        !
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>grd_ng) then
              ig=grd_ng
              wlsrc=grd_wl(grd_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=grd_wl(1)
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
     if ((grd_isvelocity).and.(ptcl%rtsrc==1)) then
        select case(in_igeom)
!-- 1D
        case(1)
           call advection1(.true.,ig,zsrc,rsrc)
!-- 2D
        case(2)
           call advection2(.true.,ig,zsrc,iy,rsrc,y)
!-- 3D
        case(3)
           call advection3(.true.,ig,zsrc,iy,iz,rsrc,y,z)
        endselect
     endif
     endif

!     write(*,*) ipart
!-----------------------------------------------------------------------        
!-- Advancing particle until census, absorption, or escape from domain
!Calling either diffusion or transport depending on particle type (ptcl%rtsrc)
     select case(in_igeom)

!-- 1D
     case(1)
        do while ((.not.prt_done).and.(.not.isvacant))
           if (ptcl%rtsrc == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport1(ptcl,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion1(ptcl,isvacant)
           endif
!-- transformation factor
           if(grd_isvelocity .and. ptcl%rtsrc==1) then
              labfact = 1.0d0 - musrc*rsrc/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (esrc<1d-6*ptcl%ebirth .and. .not.isvacant) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(zsrc,iy,iz) = grd_edep(zsrc,iy,iz) + esrc*labfact
!-- velocity effects accounting
                 if(ptcl%rtsrc==1) tot_evelo = tot_evelo + esrc*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + esrc
!
                 esrc = 2d0*esrc
                 ptcl%ebirth = 2d0*ptcl%ebirth
              endif
           endif
        enddo

!-- 2D
     case(2)
        do while ((.not.prt_done).and.(.not.isvacant))
           if (ptcl%rtsrc == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport2(ptcl,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion2(ptcl,isvacant)
           endif
!-- transformation factor
           if(grd_isvelocity .and. ptcl%rtsrc==1) then
              labfact = 1d0-(musrc*y+sqrt(1d0-musrc**2) * &
                   cos(om)*rsrc)/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (esrc<1d-6*ptcl%ebirth .and. .not.isvacant) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(zsrc,iy,iz) = grd_edep(zsrc,iy,iz) + esrc*labfact
!-- velocity effects accounting
                 if(ptcl%rtsrc==1) tot_evelo = tot_evelo + esrc*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + esrc
!
                 esrc = 2d0*esrc
                 ptcl%ebirth = 2d0*ptcl%ebirth
              endif
           endif
        enddo

!-- 3D
     case(3)
        do while ((.not.prt_done).and.(.not.isvacant))
           if (ptcl%rtsrc == 1.or.in_puretran) then
              nimc = nimc + 1
              call transport3(ptcl,isvacant)
           else
              nddmc = nddmc + 1
              call diffusion3(ptcl,isvacant)
           endif
!-- transformation factor
           if(grd_isvelocity .and. ptcl%rtsrc==1) then
              labfact = 1d0-(musrc*z+sqrt(1d0-musrc**2) * &
                   (cos(om)*rsrc+sin(om)*y))/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (esrc<1d-6*ptcl%ebirth .and. .not.isvacant) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 isvacant = .true.
                 prt_done = .true.
                 grd_edep(zsrc,iy,iz) = grd_edep(zsrc,iy,iz) + esrc*labfact
!-- velocity effects accounting
                 if(ptcl%rtsrc==1) tot_evelo = tot_evelo + esrc*(1d0-labfact)
              else
!-- weight addition accounted for in external source
                 tot_eext = tot_eext + esrc
!
                 esrc = 2d0*esrc
                 ptcl%ebirth = 2d0*ptcl%ebirth
              endif
           endif
        enddo

     endselect

!-----------------------------------------------------------------------


     if(.not.isvacant) then

     ! Redshifting DDMC particle energy weights and wavelengths
     if(ptcl%rtsrc == 2.and.grd_isvelocity) then
!-- redshifting energy weight!{{{
        tot_evelo=tot_evelo+esrc*(1d0-exp(-tsp_dt/tsp_t))
        esrc = esrc*exp(-tsp_dt/tsp_t)
        ptcl%ebirth = ptcl%ebirth*exp(-tsp_dt/tsp_t)
        !
!
!-- find group
        ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)
        !
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grd_ng) then
              ig=grd_ng
           else
              ig=1
           endif
        endif
        !
        !
        if(ig<grd_ng) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           x1 = grd_cap(ig,zsrc,iy,iz)
           x2 = grd_wl(ig)/(pc_c*tsp_t*(grd_wl(ig+1)-grd_wl(ig)))
           if(r1<x2/(x1+x2)) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              wlsrc = 1d0/(r1/grd_wl(ig+1)+(1d0-r1)/grd_wl(ig))
              wlsrc = wlsrc*exp(tsp_dt/tsp_t)
           endif
        endif
        !!}}}
     endif

     endif

     ! Looking up group
     if(ptcl%rtsrc==1) then
        if(grd_isvelocity) then!{{{
           ig = binsrch(wlsrc/labfact,grd_wl,grd_ng+1,in_ng)
        else
           ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)
        endif
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid energy bound
           if(ig>grd_ng) then
              ig=grd_ng
              if(grd_isvelocity) then
                 wlsrc=grd_wl(grd_ng+1)*labfact
              else
                 wlsrc=grd_wl(grd_ng+1)
              endif
           elseif(ig<1) then
              ig=1
              if(grd_isvelocity) then
                 wlsrc=grd_wl(1)*labfact
              else
                 wlsrc=grd_wl(1)
              endif
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     else
        ig = binsrch(wlsrc,grd_wl,grd_ng+1,in_ng)!{{{
        !
        if(ig>grd_ng.or.ig<1) then
           !particle out of wlgrid bound
           if(ig>grd_ng) then
              ig=grd_ng
              wlsrc=grd_wl(grd_ng+1)
           elseif(ig<1) then
              ig=1
              wlsrc=grd_wl(1)
           else
              write(*,*) 'domain leak!!'
              prt_done = .true.
              isvacant = .true.
           endif
        endif
        !!}}}
     endif
     
     if(isshift) then
     if ((grd_isvelocity).and.(ptcl%rtsrc==1)) then
        select case(in_igeom)
!-- 1D
        case(1)
           call advection1(.false.,ig,zsrc,rsrc)
!-- 2D
        case(2)
           call advection2(.false.,ig,zsrc,iy,rsrc,y)
!-- 3D
        case(3)
           call advection3(.false.,ig,zsrc,iy,iz,rsrc,y,z)
        endselect
     endif
     endif

     if(.not.isvacant) then
!
!-- radiation energy at census
     if(grd_isvelocity) then
        if(ptcl%rtsrc==2) then
           tot_erad = tot_erad + esrc
        else
           tot_erad = tot_erad + esrc !*(1d0-musrc*rsrc/pc_c)
!-- velocity effects accounting
!           tot_evelo=tot_evelo+esrc*musrc*rsrc/pc_c
!
        endif
     else
        tot_erad = tot_erad + esrc
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

  tot_eext = tot_eext-tot_eleft-tot_eright

end subroutine particle_advance


              !wlsrc = 0.5d0*(grd_wl(ig)+grd_wl(ig+1))
              !r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
              !wlsrc=grd_wl(ig)*(1d0-r1)+grd_wl(ig+1)*r1
              !
!               r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!               if(r1<grd_cap(ig,zsrc,1,1)/(grd_cap(ig,zsrc,1,1)+grd_sig(zsrc,1,1))) then
!                  x1 = pc_h*pc_c/(grd_wl(ig+1)*pc_kb*grd_temp(zsrc,1,1))
!                  x2 = pc_h*pc_c/(grd_wl(ig)*pc_kb*grd_temp(zsrc,1,1))
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
!                  wlsrc = pc_h*pc_c/(xx0*pc_kb*grd_temp(zsrc,1,1))
!               else
