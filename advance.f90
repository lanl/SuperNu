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

  integer :: ipart, difs, transps
  real*8 :: r1, alph2
  integer, pointer :: zsrc, rtsrc, gsrc
  real*8, pointer :: rsrc, musrc, tsrc, esrc, ebirth
  logical, pointer :: isvacant
  real :: t0,t1  !timing

  gas_edep = 0.0
  gas_erad = 0.0
  gas_eright = 0.0
  gas_eleft = 0.0
  difs = 0
  transps = 0
  gas_numcensus(1:gas_nr) = 0

  call time(t0)
  ! Propagating all particles that are not considered vacant: loop
  do ipart = 1, prt_npartmax
     ! Checking vacancy
     if (prt_particles(ipart)%isvacant.eqv..false.) then
        ! Assigning pointers to corresponding particle properties
        zsrc => prt_particles(ipart)%zsrc
        gsrc => prt_particles(ipart)%gsrc
        rtsrc => prt_particles(ipart)%rtsrc
        rsrc => prt_particles(ipart)%rsrc
        musrc => prt_particles(ipart)%musrc
        tsrc => prt_particles(ipart)%tsrc
        esrc => prt_particles(ipart)%esrc
        ebirth => prt_particles(ipart)%ebirth
        isvacant => prt_particles(ipart)%isvacant

        ! Checking if particle conversions are required since prior time step
        if (in_puretran.eqv..false.) then
           if (gas_sigmapg(gsrc,zsrc)*gas_drarr(zsrc)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
              if (rtsrc == 2) then
                 r1 = rand()
                 rsrc = (r1*gas_rarr(zsrc + 1)**3 + (1.0 - r1)*gas_rarr(zsrc)**3)**(1.0/3.0)
                 r1 = rand()
                 musrc = 1.0 - 2.0*r1
                 musrc = (musrc + gas_velyes*rsrc/pc_c)/(1.0 + gas_velyes*rsrc*musrc/pc_c)
                 esrc = esrc/(1.0 - gas_velyes*musrc*rsrc/pc_c)
                 ebirth = ebirth/(1.0 - gas_velyes*musrc*rsrc/pc_c)
              endif
              rtsrc = 1
           else
              rtsrc = 2
           endif
        endif

        ! First portion of operator split particle velocity position adjustment
        alph2 = 0.75  !>=0,<=1
        if ((in_isvelocity.eqv..true.).and.(rtsrc==1)) then
           rsrc = rsrc*tsp_texp/(tsp_texp + alph2*tsp_dt)
           if (rsrc < gas_rarr(zsrc)) then
              zsrc = zsrc - 1
           endif
        endif

        ! Advancing particle until census, absorption, or escape from domain
        prt_done = .false.
        do while (prt_done .eqv. .false.)
           !Calling either diffusion or transport depending on particle type (rtsrc)
           if (rtsrc == 1) then
              transps = transps + 1
              call transport1(zsrc,gsrc,rsrc,musrc,tsrc,esrc,ebirth,rtsrc,isvacant)
           else
              difs = difs + 1
              call diffusion1(zsrc,gsrc,rsrc,musrc,tsrc,esrc,ebirth,rtsrc,isvacant)
           endif
        enddo

        if ((in_isvelocity.eqv..true.).and.(rtsrc==1)) then
           rsrc = rsrc*tsp_texp/(tsp_texp + (1.0 - alph2)*tsp_dt)
           if (rsrc < gas_rarr(zsrc)) then
              zsrc = zsrc - 1
           endif
        endif

     endif
  enddo
  call time(t1)
  call timereg(t_pckt, t1-t0)  !register timing
  !write(6,*) transps, difs
end subroutine advance
