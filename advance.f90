SUBROUTINE advance

  USE data_mod
  IMPLICIT NONE
  
  INTEGER(iknd) :: ipart, difs, transps
  REAL(rknd) :: r1, alph2
  INTEGER(iknd), POINTER :: zsrc, rtsrc, gsrc
  REAL(rknd), POINTER :: rsrc, musrc, tsrc, Esrc, Ebirth
  LOGICAL, POINTER :: isvacant 

  Edep = 0.0
  Erad = 0.0
  Eright = 0.0
  Eleft = 0.0
  difs = 0
  transps = 0
  numcensus(1:gas_nr) = 0
  
  DO ipart = 1, prt_Npartmax

     IF (particles(ipart)%isvacant.EQV..FALSE.) THEN

        zsrc => particles(ipart)%zsrc
        gsrc => particles(ipart)%gsrc
        rtsrc => particles(ipart)%rtsrc
        rsrc => particles(ipart)%rsrc
        musrc => particles(ipart)%musrc
        tsrc => particles(ipart)%tsrc
        Esrc => particles(ipart)%Esrc
        Ebirth => particles(ipart)%Ebirth
        isvacant => particles(ipart)%isvacant

        IF (puretran.EQV..FALSE.) THEN
           IF (sigmapg(gsrc,zsrc)*drarr(zsrc)*(velno*1.0+velyes*texp)<5.0_rknd) THEN
              IF (rtsrc == 2) THEN
                 r1 = RAND()
                 rsrc = (r1*rarr(zsrc+1)**3+(1.0-r1)*rarr(zsrc)**3)**(1.0/3.0)
                 r1 = RAND()
                 musrc = 1.0-2.0*r1
                 musrc = (musrc+velyes*rsrc/lspeed)/(1.0+velyes*rsrc*musrc/lspeed)
                 Esrc = Esrc/(1.0-velyes*musrc*rsrc/lspeed)
                 Ebirth = Ebirth/(1.0-velyes*musrc*rsrc/lspeed)
              ENDIF
              rtsrc = 1
           ELSE
              rtsrc = 2
           ENDIF
        ENDIF

        alph2 = 0.75  !>=0,<=1
        IF ((isvelocity.EQV..TRUE.).AND.(rtsrc==1)) THEN
           rsrc = rsrc*texp/(texp+alph2*dt)
           IF (rsrc < rarr(zsrc)) THEN
              zsrc = zsrc-1
           ENDIF
        ENDIF

        done = .FALSE.
        DO WHILE (done .EQV. .FALSE.)
           !Calling either diffusion or transport depending on particle type (rtsrc)
           IF (rtsrc == 1) THEN
              transps=transps+1
              CALL transport1(zsrc,gsrc,rsrc,musrc,tsrc,Esrc,Ebirth,rtsrc,isvacant)
           ELSE
              difs=difs+1
              CALL diffusion1(zsrc,gsrc,rsrc,musrc,tsrc,Esrc,Ebirth,rtsrc,isvacant)
           ENDIF
        ENDDO

        IF ((isvelocity.EQV..TRUE.).AND.(rtsrc==1)) THEN
           rsrc = rsrc*texp/(texp+(1.0-alph2)*dt)
           IF (rsrc < rarr(zsrc)) THEN
              zsrc = zsrc-1
           ENDIF
        ENDIF
   
     ENDIF
  ENDDO
  !WRITE(*,*) transps, difs
END SUBROUTINE advance
