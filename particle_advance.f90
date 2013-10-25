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

  integer :: ipart, difs, transps, g, zholder, zfdiff, ir, binsrch
  real*8 :: r1, alph2, r2, x1, x2, xx0, bmax, help, P
  real*8 :: uul, uur, uumax, r0, r3, mu0, rold, dddd
  integer, pointer :: zsrc, rtsrc !, gsrc
  real*8, pointer :: rsrc, musrc, tsrc, esrc, ebirth, wlsrc
  logical, pointer :: isvacant
  real :: t0,t1  !timing

  logical :: isshift=.true.
  logical :: partstopper=.true.
  logical :: showidfront=.false.

  gas_edep = 0.0
  gas_erad = 0.0
  gas_eright = 0.0
  gas_luminos = 0.0
  gas_eleft = 0.0
!--(rev. 121)
  !gas_eraddens =0d0
!--
  difs = 0
  transps = 0
  gas_numcensus(1:gas_nr) = 0
!-- advection split parameter
  alph2 = 0.5d0  !>=0,<=1
!--

  if(showidfront) then
     do ir = 1, gas_nr-1
        if(gas_isvelocity.and.(gas_sig(ir)+gas_cap(1,ir))*gas_drarr(ir) &
             *tsp_texp>=prt_tauddmc &
             .and. &
             (gas_sig(ir+1)+gas_cap(1,ir+1))*gas_drarr(ir+1) &
             *tsp_texp<prt_tauddmc) then
           write(*,*) ir, gas_cap(1,ir)*gas_drarr(ir)*tsp_texp, &
                gas_cap(1,ir+1)*gas_drarr(ir+1)*tsp_texp
        endif
     enddo
  endif
  
  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  do ipart = 1, prt_npartmax
     ! Checking vacancy
     if (prt_particles(ipart)%isvacant) cycle

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


     ! Looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then
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
        !
     else
        g = binsrch(wlsrc,gas_wl,gas_ng+1)
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
        !
     endif

     
     ! Checking if particle conversions are required since prior time step
     if(.not.in_puretran) then
        if(gas_isvelocity) then
           help = tsp_texp
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
!-- converting to lab frame (doppler and aberration)
!                if(gas_isvelocity) then
!                   esrc = (1d0+0.7d0*rsrc/pc_c)*esrc
!                   ebirth = (1d0+0.7d0*rsrc/pc_c)*ebirth
!                endif
              if(gas_isvelocity) then
                 musrc = (musrc + rsrc/pc_c)/(1.0 + rsrc*musrc/pc_c)
                 esrc = esrc/(1.0 - musrc*rsrc/pc_c)
                 ebirth = ebirth/(1.0 - musrc*rsrc/pc_c)
              endif
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
        endif
     endif 
!
!-- looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then
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
        !
     else
        g = binsrch(wlsrc,gas_wl,gas_ng+1)
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
        !
     endif

!-- First portion of operator split particle velocity position adjustment
     if(isshift) then
     if ((gas_isvelocity).and.(rtsrc==1)) then
        rold = rsrc
        rsrc = rsrc*tsp_texp/(tsp_texp+alph2*tsp_dt)
        !
        if (rsrc < gas_rarr(zsrc)) then
           !
           zholder = binsrch(rsrc,gas_rarr,gas_nr+1)
           !
           if(gas_isshell.and.zsrc==1) then
              prt_done = .true.
              isvacant = .true.
           elseif(.not.in_puretran.and.partstopper) then
              zfdiff = -1
              if(gas_isvelocity) then
                 help = tsp_texp
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
!                  mu0 = (musrc-rsrc/pc_c)/(1d0-rsrc*musrc/pc_c)
!                  if(musrc>=0d0.and.mu0<0d0) then
! !--
!                     ebirth=ebirth*(1d0+2d0*(0.55d0/abs(mu0)-1.3d0*abs(mu0))*rsrc/pc_c)
!                     esrc=esrc*(1d0+2d0*(0.55d0/abs(mu0)-1.3d0*abs(mu0))*rsrc/pc_c)
! !--
!                     r1 = rand()
!                     prt_tlyrand = prt_tlyrand+1
!                     P=gas_ppr(g,zsrc-1)*(1d0+1.5d0*abs(mu0))
!                     if(r1<P) then
!                        zsrc = zfdiff
!                        rtsrc = 2
!                        esrc = esrc*(1d0-rsrc*musrc/pc_c)
!                        ebirth = ebirth*(1d0-rsrc*musrc/pc_c)
!                        wlsrc = wlsrc/(1d0-rsrc*musrc/pc_c)
!                     else
!                        r1 = rand()
!                        prt_tlyrand = prt_tlyrand+1
!                        r2 = rand()
!                        prt_tlyrand = prt_tlyrand+1
!                        mu0 = max(r1,r2)
!                        if(gas_isvelocity) then
!                           musrc = (mu0+rsrc/pc_c)/(1.0+rsrc*mu0/pc_c)
!                        endif
!                     endif
! !--
!                  endif
!--
              else
                 zsrc = zholder
              endif
           else
              zsrc = zholder
              if((gas_sig(zsrc)+gas_cap(g,zsrc))*gas_drarr(zsrc) &
                   *tsp_texp>=prt_tauddmc*gas_curvcent(zsrc)) then
                 rtsrc = 2
              endif
           endif
           !
        endif
     endif
        !
     endif

     !if(rtsrc==1) then
     !   write(*,*) g,zsrc,wlsrc,rsrc
     !endif
!-----------------------------------------------------------------------        
     ! Advancing particle until census, absorption, or escape from domain
     prt_done = .false.
     do while (.not.prt_done)
        !Calling either diffusion or transport depending on particle type (rtsrc)
        if (rtsrc == 1.or.in_puretran) then
           transps = transps + 1
           call transport1(zsrc,wlsrc,rsrc,musrc,tsrc, &
                esrc,ebirth,rtsrc,isvacant)
        else
           difs = difs + 1
           call diffusion1(zsrc,wlsrc,rsrc,musrc,tsrc, &
                esrc,ebirth,rtsrc,isvacant)
        endif
     enddo
!-----------------------------------------------------------------------
     !---------------
     !------------

     ! Redshifting DDMC particle energy weights and wavelengths
!      if(rtsrc == 2.and.gas_isvelocity) then
! !-- redshifting energy weight
!         esrc = esrc*exp(-tsp_dt/tsp_texp)
!         ebirth = ebirth*exp(-tsp_dt/tsp_texp)
!         !
! !
! !-- find group
! !         g = binsrch(wlsrc,gas_wl,gas_ng+1)
! !         !
! !         if(g>gas_ng.or.g<1) then
! !            !particle out of wlgrid energy bound
! !            if(g>gas_ng) then
! !               g=gas_ng
! !            else
! !               g=1
! !            endif
! !         endif
! !         !
! !         !
! !         r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
! !         wlsrc = 1d0/(r1/gas_wl(g+1)+(1d0-r1)/gas_wl(g))
! !         wlsrc = wlsrc*exp(tsp_dt/tsp_texp)
!         !
!      endif

     ! Looking up group
     if(rtsrc==1) then
        if(gas_isvelocity) then
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
        !
     else
        g = binsrch(wlsrc,gas_wl,gas_ng+1)
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
        !
     endif
     
     if(isshift) then
     if ((gas_isvelocity).and.(rtsrc==1)) then
        !
        rold = rsrc
        rsrc = rsrc*(tsp_texp+alph2*tsp_dt)/(tsp_texp+tsp_dt)
        !
        if (rsrc < gas_rarr(zsrc)) then
           !
           zholder = binsrch(rsrc,gas_rarr,gas_nr+1)
           !
           if(gas_isshell.and.zsrc==1) then
              prt_done = .true.
              isvacant = .true.
           elseif(.not.in_puretran.and.partstopper) then
              zfdiff = -1
              if(gas_isvelocity) then
                 help = tsp_texp
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
!                  mu0 = (musrc-rsrc/pc_c)/(1d0-rsrc*musrc/pc_c)
!                  if(musrc>=0d0.and.mu0<0d0) then
! !--
!                     ebirth=ebirth*(1d0+2d0*(0.55d0/abs(mu0)-1.3d0*abs(mu0))*rsrc/pc_c)
!                     esrc=esrc*(1d0+2d0*(0.55d0/abs(mu0)-1.3d0*abs(mu0))*rsrc/pc_c)
! !--
!                     r1 = rand()
!                     prt_tlyrand = prt_tlyrand+1
!                     P=gas_ppr(g,zsrc-1)*(1d0+1.5d0*abs(mu0))
!                     if(r1<P) then
!                        zsrc = zfdiff
!                        rtsrc = 2
!                        esrc = esrc*(1d0-rsrc*musrc/pc_c)
!                        ebirth = ebirth*(1d0-rsrc*musrc/pc_c)
!                        wlsrc = wlsrc/(1d0-rsrc*musrc/pc_c)
!                     else
!                        r1 = rand()
!                        prt_tlyrand = prt_tlyrand+1
!                        r2 = rand()
!                        prt_tlyrand = prt_tlyrand+1
!                        mu0 = max(r1,r2)
!                        if(gas_isvelocity) then
!                           musrc = (mu0+rsrc/pc_c)/(1.0+rsrc*mu0/pc_c)
!                        endif
!                     endif
! !--
!                 endif
!--
              else
                 zsrc = zholder
              endif
           else
              zsrc = zholder
           endif
           !
        endif
        !
     endif
     endif
  
  enddo

  call time(t1)
  call timereg(t_pckt, t1-t0)  !register timing
  write(6,*) transps, difs
  !write(6,*) eleft, eright

end subroutine particle_advance
