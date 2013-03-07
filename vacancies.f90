SUBROUTINE vacancies

  USE particlemod
  IMPLICIT NONE

!##################################################
  !This subroutine creates an array of vacant index locations
  !in the particle array each time step.
!##################################################

  INTEGER :: ipart, ivac
  LOGICAL :: isfull

  !Initializing index counters and full array checking boolean
  isfull = .FALSE.
  ipart = 0
  ivac = 0

  !Filling prt_vacantarr with particle index of vacant particles: loop
  DO WHILE (isfull.EQV..FALSE.)
     ipart = ipart+1
     IF (prt_particles(ipart)%isvacant.EQV..TRUE.) THEN
        ivac = ivac+1
        prt_vacantarr(ivac) = ipart
     ENDIF
     IF (ivac == prt_nnew) THEN
        isfull = .TRUE.
     ELSEIF (ipart == prt_npartmax) THEN
        WRITE(*,*) 'Maximum number of prt_particles reached'
        isfull = .TRUE.
     ENDIF
  ENDDO
  
END SUBROUTINE vacancies
