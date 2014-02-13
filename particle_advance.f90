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

  integer*8 :: nddmc, nimc, npckt
  integer :: ipart, g, zholder, zfdiff, ir, binsrch
  real*8 :: r1, alph2, r2, x1, x2, xx0, bmax, help, P
  real*8 :: uul, uur, uumax, r0, r3, mu0, rold, dddd
  integer, pointer :: zsrc, rtsrc !, gsrc
  real*8, pointer :: rsrc, musrc, tsrc, esrc, ebirth, wlsrc
  logical, pointer :: isvacant
  real*8 :: t0,t1  !timing

  logical :: isshift=.true.
  logical :: partstopper=.true.
  logical :: showidfront=.false.

  gas_edep = 0.0
  gas_erad = 0.0
  gas_luminos = 0.0

!
!--(rev. 121)
  gas_eraddens =0d0
!--
  gas_numcensus(1:gas_nr) = 0
!-- advection split parameter
  alph2 = 0.5d0  !>=0,<=1
!--

  if(showidfront) then
     do ir = 1, gas_nr-1
        if(gas_isvelocity.and.(gas_sig(ir)+gas_cap(1,ir))*gas_drarr(ir) &
             *tsp_t>=prt_tauddmc &
             .and. &
             (gas_sig(ir+1)+gas_cap(1,ir+1))*gas_drarr(ir+1) &
             *tsp_t<prt_tauddmc) then
           write(*,*) ir, gas_cap(1,ir)*gas_drarr(ir)*tsp_t, &
                gas_cap(1,ir+1)*gas_drarr(ir+1)*tsp_t
        endif
     enddo
  endif
  
  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  npckt = 0
  nddmc = 0
  nimc = 0
  do ipart = 1, prt_npartmax
     ! Checking vacancy!{{{
     if (prt_particles(ipart)%isvacant) cycle
     npckt = npckt + 1

     ! Assigning pointers to corresponding particle properties
     zsrc => prt_particles(ipart)%zsrc
     !
     !Ryan W.: Replacing particle group with particle wavelength (rev. 120)
     wlsrc => prt_particles(ipart)%wlsrc
     !gsrc => prt_particles(ipart)%gsrc
     !
     rtsrc => prt_particles(ipart)%rtsrc
     rsrc => prt_particles(ipart)%rsrc
     musrc => prt_particles(ipart)%musrc
     tsrc => prt_particles(ipart)%tsrc
     esrc => prt_particles(ipart)%esrc
     ebirth => prt_particles(ipart)%ebirth
     isvacant => prt_particles(ipart)%isvacant

     prt_done=.false.

     ! Looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then!{{{
           g = binsrch(wlsrc/(1d0-rsrc*musrc/pc_c),gas_wl,gas_ng+1)
        else
           g = binsrch(wlsrc,gas_wl,gas_ng+1)
        endif
        !
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid energy bound
           if(g>gas_ng) then
              g=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(g<1) then
              g=1
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
        g = binsrch(wlsrc,gas_wl,gas_ng+1)!{{{
        !
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid bound
           if(g>gas_ng) then
              g=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(g<1) then
              g=1
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
        if ((gas_sig(zsrc)+gas_cap(g,zsrc))*gas_drarr(zsrc) &
             *help<prt_tauddmc*gas_curvcent(zsrc)) then
           !write(*,*) 'here', g, wlsrc, esrc
           if (rtsrc == 2) then
!-- sampling position uniformly
!              r1 =  rand()
!           prt_tlyrand = prt_tlyrand+1
!              rsrc = (r1*gas_rarr(zsrc+1)**3 + (1.0-r1)*gas_rarr(zsrc)**3)**(1.0/3.0)
!-- sampling position from source tilt
              r1 = 0d0
              r2 = 1d0
              uul = gas_tempb(zsrc)**4
              uur = gas_tempb(zsrc+1)**4
              uumax = max(uul,uur)
              do while (r2 > r1)
                 r3 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 r0 = (r3*gas_rarr(zsrc+1)**3+(1.0-r3)*gas_rarr(zsrc)**3)**(1.0/3.0)
                 r3 = (r0-gas_rarr(zsrc))/gas_drarr(zsrc)
                 r1 = (r3*uur+(1d0-r3)*uul)/uumax
                 r2 = rand()
                 prt_tlyrand = prt_tlyrand+1
              enddo
              rsrc = r0
!
!-- sampling angle isotropically
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              musrc = 1.0 - 2.0*r1
              if(gas_isvelocity) then
                 musrc = (musrc + rsrc/pc_c)/(1.0 + rsrc*musrc/pc_c)
!-- velocity effects accounting
                 gas_evelo=gas_evelo+esrc*(1d0-1d0/(1.0 - musrc*rsrc/pc_c))
!
                 esrc = esrc/(1.0 - musrc*rsrc/pc_c)
                 ebirth = ebirth/(1.0 - musrc*rsrc/pc_c)
              endif
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 wlsrc = 1d0/(r1/gas_wl(g+1)+(1d0-r1)/gas_wl(g))
!              endif
              !
              if(gas_isvelocity) then
                 wlsrc = wlsrc*(1.0-musrc*rsrc/pc_c)
              endif
           endif
           rtsrc = 1
        else
           rtsrc = 2
        endif!}}}
     endif 
!
!-- looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then!{{{
           g = binsrch(wlsrc/(1.0d0-rsrc*musrc/pc_c),gas_wl,gas_ng+1)
        else
           g = binsrch(wlsrc,gas_wl,gas_ng+1)
        endif
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid energy bound
           if(g>gas_ng) then
              g=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(g<1) then
              g=1
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
        g = binsrch(wlsrc,gas_wl,gas_ng+1)!{{{
        !
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid bound
           if(g>gas_ng) then
              g=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(g<1) then
              g=1
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
     if ((gas_isvelocity).and.(rtsrc==1)) then
        rold = rsrc!{{{
        rsrc = rsrc*tsp_t/(tsp_t+alph2*tsp_dt)
        !
        if (rsrc < gas_rarr(zsrc)) then
           !
           zholder = binsrch(rsrc,gas_rarr,gas_nr+1)
           !
           if(gas_isshell.and.zsrc==1) then
              prt_done = .true.
              isvacant = .true.
              gas_eleft = gas_eleft+esrc*(1d0-musrc*rsrc/pc_c)
!-- velocity effects accounting
              gas_evelo = gas_evelo+esrc*musrc*rsrc/pc_c
!
           elseif(.not.in_puretran.and.partstopper) then
              zfdiff = -1
              if(gas_isvelocity) then
                 help = tsp_t
              else
                 help = 1d0
              endif
              do ir = zsrc-1,zholder,-1
                 if((gas_sig(ir)+gas_cap(g,ir))*gas_drarr(ir) &
                      *help>=prt_tauddmc*gas_curvcent(ir)) then
                    zfdiff = ir
                    exit
                 endif
              enddo
              if(zfdiff.ne.-1) then
!--
                 zsrc = zfdiff+1
                 rsrc = gas_rarr(zsrc)
!--
              else
                 zsrc = zholder
              endif
           else
              zsrc = zholder
           endif
           !
        endif!}}}
     endif
        !
     endif

!     write(*,*) ipart
!-----------------------------------------------------------------------        
     ! Advancing particle until census, absorption, or escape from domain
     do while ((.not.prt_done).and.(.not.isvacant))
        !Calling either diffusion or transport depending on particle type (rtsrc)!{{{
        if (rtsrc == 1.or.in_puretran) then
           nimc = nimc + 1
           call transport1(zsrc,wlsrc,rsrc,musrc,tsrc, &
                esrc,ebirth,rtsrc,isvacant,ipart)
        else
           nddmc = nddmc + 1
           call diffusion1(zsrc,wlsrc,rsrc,musrc,tsrc, &
                esrc,ebirth,rtsrc,isvacant,ipart)
        endif!}}}
     enddo
!-----------------------------------------------------------------------
     !---------------
     !------------

     if(.not.isvacant) then

     ! Redshifting DDMC particle energy weights and wavelengths
     if(rtsrc == 2.and.gas_isvelocity) then
!-- redshifting energy weight
        gas_evelo=gas_evelo+esrc*(1d0-exp(-tsp_dt/tsp_t))
        esrc = esrc*exp(-tsp_dt/tsp_t)
        ebirth = ebirth*exp(-tsp_dt/tsp_t)
        !
!
!-- find group
        g = binsrch(wlsrc,gas_wl,gas_ng+1)
        !
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid energy bound
           if(g>gas_ng) then
              g=gas_ng
           else
              g=1
           endif
        endif
        !
        !
        if(g<gas_ng) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
!           if(r1<gas_sig(zsrc)/(gas_cap(g,zsrc)+gas_sig(zsrc))) then
              wlsrc = 1d0/(r1/gas_wl(g+1)+(1d0-r1)/gas_wl(g))
!           wlsrc = r1*gas_wl(g)+(1d0-r1)*gas_wl(g+1)
              wlsrc = wlsrc*exp(tsp_dt/tsp_t)
!           endif
        endif
        !
     endif

     endif

     ! Looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then!{{{
           g = binsrch(wlsrc/(1.0d0-rsrc*musrc/pc_c),gas_wl,gas_ng+1)
        else
           g = binsrch(wlsrc,gas_wl,gas_ng+1)
        endif
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid energy bound
           if(g>gas_ng) then
              g=gas_ng
              if(gas_isvelocity) then
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-rsrc*musrc/pc_c)
              else
                 wlsrc=gas_wl(gas_ng+1)
              endif
           elseif(g<1) then
              g=1
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
        g = binsrch(wlsrc,gas_wl,gas_ng+1)!{{{
        !
        if(g>gas_ng.or.g<1) then
           !particle out of wlgrid bound
           if(g>gas_ng) then
              g=gas_ng
              wlsrc=gas_wl(gas_ng+1)
           elseif(g<1) then
              g=1
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
     if ((gas_isvelocity).and.(rtsrc==1)) then
        !!{{{
        rold = rsrc
        rsrc = rsrc*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
        !
        if (rsrc < gas_rarr(zsrc)) then
           !
           zholder = binsrch(rsrc,gas_rarr,gas_nr+1)
           !
           if(gas_isshell.and.zsrc==1) then
              if(.not.isvacant) then
              prt_done = .true.
              isvacant = .true.
              gas_eleft=gas_eleft+esrc*(1d0-musrc*rsrc/pc_c)
!-- velocity effects accounting
              gas_evelo=gas_evelo+esrc*musrc*rsrc/pc_c
!
              endif
           elseif(.not.in_puretran.and.partstopper) then
              zfdiff = -1
              if(gas_isvelocity) then
                 help = tsp_t
              else
                 help = 1d0
              endif
              do ir = zsrc-1,zholder,-1                 
                 if((gas_sig(ir)+gas_cap(g,ir))*gas_drarr(ir) &
                      *help>=prt_tauddmc*gas_curvcent(ir)) then
                    zfdiff = ir
                    exit
                 endif
              enddo
              if(zfdiff.ne.-1) then
!--
                 zsrc = zfdiff+1
                 rsrc = gas_rarr(zsrc)
!--
              else
                 zsrc = zholder
              endif
           else
              zsrc = zholder
           endif
           !
        endif
        !!}}}
     endif
     endif

     if(.not.isvacant) then
!
!-- radiation energy at census
     if(gas_isvelocity) then
        if(rtsrc==2) then
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


              !wlsrc = 0.5d0*(gas_wl(g)+gas_wl(g+1))
              !r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
              !wlsrc=gas_wl(g)*(1d0-r1)+gas_wl(g+1)*r1
              !
!               r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!               if(r1<gas_cap(g,zsrc)/(gas_cap(g,zsrc)+gas_sig(zsrc))) then
!                  x1 = pc_h*pc_c/(gas_wl(g+1)*pc_kb*gas_temp(zsrc))
!                  x2 = pc_h*pc_c/(gas_wl(g)*pc_kb*gas_temp(zsrc))
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
!                  wlsrc = pc_h*pc_c/(xx0*pc_kb*gas_temp(zsrc))
!               else
