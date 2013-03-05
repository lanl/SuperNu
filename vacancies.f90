SUBROUTINE vacancies

  USE data_mod
  IMPLICIT NONE

  INTEGER(iknd) :: ipart, ivac
  LOGICAL :: isfull

  isfull = .FALSE.
  ipart = 0
  ivac = 0

  !ALLOCATE(vacantarr(Nnew))
  DO WHILE (isfull.EQV..FALSE.)
     ipart = ipart+1
     IF (particles(ipart)%isvacant.EQV..TRUE.) THEN
        ivac = ivac+1
        vacantarr(ivac) = ipart
     ENDIF
     IF (ivac == Nnew) THEN
        isfull = .TRUE.
     ELSEIF (ipart == Npartmax) THEN
        WRITE(*,*) 'Maximum number of particles reached'
        isfull = .TRUE.
     ENDIF
  ENDDO
  
END SUBROUTINE vacancies
