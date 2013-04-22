subroutine advance

  use particlemod
  use timestepmod
  use gasgridmod
  use physconstmod
  use inputparmod
  use timingmod
  implicit none

!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
!##################################################

  integer :: ipart, difs, transps, g, zholder, zfdiff, ir
  real*8 :: r1, alph2
  integer, pointer :: zsrc, rtsrc !, gsrc
  real*8, pointer :: rsrc, musrc, tsrc, esrc, ebirth, wlsrc
  logical, pointer :: isvacant
  real :: t0,t1  !timing

  logical :: isshift=.true.

  gas_edep = 0.0
  gas_erad = 0.0
  gas_eright = 0.0
  gas_eleft = 0.0
!--(rev. 121)
  !gas_eraddensg =0d0
!--
  difs = 0
  transps = 0
  gas_numcensus(1:gas_nr) = 0
  alph2 = 0.5d0  !>=0,<=1

  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  do ipart = 1, prt_npartmax
     ! Checking vacancy
     if (prt_particles(ipart)%isvacant.eqv..false.) then
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
           g = minloc(abs(gas_wl-wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)),1)
           if(wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)-gas_wl(g)<0d0) then
              g = g-1
           endif
        else
           g = minloc(abs(gas_wl-wlsrc),1)
           if(wlsrc-gas_wl(g)<0d0) then
              g = g-1
           endif
        endif
                        
        ! Checking if particle conversions are required since prior time step
        if (in_puretran.eqv..false.) then
           if ((gas_sig(zsrc)+gas_sigmapg(g,zsrc))*gas_drarr(zsrc) &
                *(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
              if (rtsrc == 2) then
                 r1 =  rand()
                 rsrc = (r1*gas_rarr(zsrc+1)**3 + (1.0-r1)*gas_rarr(zsrc)**3)**(1.0/3.0)
                 r1 = rand()
                 musrc = 1.0 - 2.0*r1
                 musrc = (musrc + gas_velyes*rsrc/pc_c)/(1.0 + gas_velyes*rsrc*musrc/pc_c)
                 esrc = esrc/(1.0 - gas_velyes*musrc*rsrc/pc_c)
                 ebirth = ebirth/(1.0 - gas_velyes*musrc*rsrc/pc_c)
                 wlsrc = wlsrc*(1.0-gas_velyes*musrc*rsrc/pc_c)
              endif
              rtsrc = 1
           else
              rtsrc = 2
           endif
        endif
        ! Looking up group
        if(rtsrc==1) then
           g = minloc(abs(gas_wl-wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)),1)
           if(wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)-gas_wl(g)<0d0) then
              g = g-1
           endif
           if(g>gas_ng.or.g<1) then
              !particle out of wlgrid energy bound
              if(g>gas_ng) then
                 g=gas_ng
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-gas_velyes*rsrc*musrc/pc_c)
              elseif(g<1) then
                 g=1
                 wlsrc=gas_wl(1)*(1.0d0-gas_velyes*rsrc*musrc/pc_c)
              else
                 write(*,*) 'domain leak!!'
                 prt_done = .true.
                 isvacant = .true.
              endif
           endif
           !
        else
           g = minloc(abs(gas_wl-wlsrc),1)
           if(wlsrc-gas_wl(g)<0d0) then
              g = g-1
           endif
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
        !if(tsp_tn==22) then
        !   write(*,*) rtsrc, wlsrc, g
        !   write(*,*) '--',gas_wl(g),gas_wl(g+1)
        !endif

        ! First portion of operator split particle velocity position adjustment
        if(isshift) then
        if ((in_isvelocity.eqv..true.).and.(rtsrc==1)) then
           rsrc = rsrc*tsp_texp/(tsp_texp+alph2*tsp_dt)
           !
           if (rsrc < gas_rarr(zsrc)) then
              !
              zholder = minloc(abs(gas_rarr-rsrc),1)
              if(rsrc<gas_rarr(zholder)) then
                 zholder = zholder-1
              endif
              if(gas_isshell.and.zsrc==1) then
                 prt_done = .true.
                 isvacant = .true.
              elseif(.not.in_puretran) then
                 zfdiff = -1
                 do ir = zsrc-1, zholder, -1
                    if((gas_sig(ir)+gas_sigmapg(g,ir))*gas_drarr(ir) &
                         *(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
                       zfdiff = ir
                       exit
                    endif
                 enddo
                 if(zfdiff==-1) then
                    zsrc = zholder
                 else
                    zsrc=zfdiff+1
                    rsrc=gas_rarr(zsrc)
                 endif
              else
                 zsrc = zholder
              endif
              !
           endif
           !
        endif
        endif
!-----------------------------------------------------------------------        
        ! Advancing particle until census, absorption, or escape from domain
        prt_done = .false.
        do while (prt_done .eqv. .false.)
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

        ! Looking up group
        if(rtsrc==1) then
           g = minloc(abs(gas_wl-wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)),1)
           if(wlsrc/(1.0d0-gas_velyes*rsrc*musrc/pc_c)-gas_wl(g)<0d0) then
              g = g-1
           endif
           if(g>gas_ng.or.g<1) then
              !particle out of wlgrid energy bound
              if(g>gas_ng) then
                 g=gas_ng
                 wlsrc=gas_wl(gas_ng+1)*(1.0d0-gas_velyes*rsrc*musrc/pc_c)
              elseif(g<1) then
                 g=1
                 wlsrc=gas_wl(1)*(1.0d0-gas_velyes*rsrc*musrc/pc_c)
              else
                 write(*,*) 'domain leak!!'
                 prt_done = .true.
                 isvacant = .true.
              endif
           endif
           !
        else
           g = minloc(abs(gas_wl-wlsrc),1)
           if(wlsrc-gas_wl(g)<0d0) then
              g = g-1
           endif
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
        if ((in_isvelocity.eqv..true.).and.(rtsrc==1)) then
           
           !
           rsrc = rsrc*tsp_texp/(tsp_texp+(1.0-alph2)*tsp_dt)
           !
           if (rsrc < gas_rarr(zsrc)) then
              !
              zholder = minloc(abs(gas_rarr-rsrc),1)
              if(rsrc<gas_rarr(zholder)) then
                 zholder = zholder-1
              endif
              if(gas_isshell.and.zsrc==1) then
                 prt_done = .true.
                 isvacant = .true.
              elseif(.not.in_puretran) then
                 zfdiff = -1
                 do ir = zsrc-1, zholder, -1
                    if((gas_sig(ir)+gas_sigmapg(g,ir))*gas_drarr(ir) &
                         *(gas_velno*1.0+gas_velyes*tsp_texp)>=5.0d0) then
                       zfdiff = ir
                       exit
                    endif
                 enddo
                 if(zfdiff==-1) then
                    zsrc = zholder
                 else
                    zsrc=zfdiff+1
                    rsrc=gas_rarr(zsrc)
                 endif
              else
                 zsrc = zholder
              endif
              !
           endif
           !
        endif
        endif
     endif
  
  enddo

  call time(t1)
  call timereg(t_pckt, t1-t0)  !register timing
  write(6,*) transps, difs
  !write(6,*) eleft, eright

end subroutine advance
