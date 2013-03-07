SUBROUTINE advance

  USE particlemod
  USE timestepmod
  USE gasgridmod
  USE physconstmod
  USE inputparmod
  IMPLICIT NONE

!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
!##################################################
  
  INTEGER :: ipart, difs, transps
  REAL*8 :: r1, alph2
  INTEGER, POINTER :: zsrc, rtsrc, gsrc
  REAL*8, POINTER :: rsrc, musrc, tsrc, Esrc, Ebirth
  LOGICAL, POINTER :: isvacant 

  gas_edep = 0.0
  prt_erad = 0.0
  prt_eright = 0.0
  prt_eleft = 0.0
  difs = 0
  transps = 0
  gas_numcensus(1:gas_nr) = 0

  ! Propagating all particles that are not considered vacant: loop
  DO ipart = 1, prt_npartmax
     ! Checking vacancy
     IF (prt_particles(ipart)%isvacant.EQV..FALSE.) THEN
        ! Assigning pointers to corresponding particle properties
        zsrc => prt_particles(ipart)%zsrc
        gsrc => prt_particles(ipart)%gsrc
        rtsrc => prt_particles(ipart)%rtsrc
        rsrc => prt_particles(ipart)%rsrc
        musrc => prt_particles(ipart)%musrc
        tsrc => prt_particles(ipart)%tsrc
        Esrc => prt_particles(ipart)%Esrc
        Ebirth => prt_particles(ipart)%Ebirth
        isvacant => prt_particles(ipart)%isvacant

        ! Checking if particle conversions are required since prior time step
        IF (in_puretran.EQV..FALSE.) THEN
           IF (gas_sigmapg(gsrc,zsrc)*gas_drarr(zsrc)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) THEN
              IF (rtsrc == 2) THEN
                 r1 = RAND()
                 rsrc = (r1*gas_rarr(zsrc+1)**3+(1.0-r1)*gas_rarr(zsrc)**3)**(1.0/3.0)
                 r1 = RAND()
                 musrc = 1.0-2.0*r1
                 musrc = (musrc+gas_velyes*rsrc/pc_c)/(1.0+gas_velyes*rsrc*musrc/pc_c)
                 Esrc = Esrc/(1.0-gas_velyes*musrc*rsrc/pc_c)
                 Ebirth = Ebirth/(1.0-gas_velyes*musrc*rsrc/pc_c)
              ENDIF
              rtsrc = 1
           ELSE
              rtsrc = 2
           ENDIF
        ENDIF

        ! First portion of operator split particle velocity position adjustment
        alph2 = 0.75  !>=0,<=1
        IF ((in_isvelocity.EQV..TRUE.).AND.(rtsrc==1)) THEN
           rsrc = rsrc*tsp_texp/(tsp_texp+alph2*tsp_dt)
           IF (rsrc < gas_rarr(zsrc)) THEN
              zsrc = zsrc-1
           ENDIF
        ENDIF

        ! Advancing particle until census, absorption, or escape from domain
        prt_done = .FALSE.
        DO WHILE (prt_done .EQV. .FALSE.)
           !Calling either diffusion or transport depending on particle type (rtsrc)
           IF (rtsrc == 1) THEN
              transps=transps+1
              CALL transport1(zsrc,gsrc,rsrc,musrc,tsrc,Esrc,Ebirth,rtsrc,isvacant)
           ELSE
              difs=difs+1
              CALL diffusion1(zsrc,gsrc,rsrc,musrc,tsrc,Esrc,Ebirth,rtsrc,isvacant)
           ENDIF
        ENDDO

        IF ((in_isvelocity.EQV..TRUE.).AND.(rtsrc==1)) THEN
           rsrc = rsrc*tsp_texp/(tsp_texp+(1.0-alph2)*tsp_dt)
           IF (rsrc < gas_rarr(zsrc)) THEN
              zsrc = zsrc-1
           ENDIF
        ENDIF
   
     ENDIF
  ENDDO
  !WRITE(*,*) transps, difs
END SUBROUTINE advance
