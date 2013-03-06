SUBROUTINE vacancies

  USE particlemod
  IMPLICIT NONE

  INTEGER(iknd) :: ipart, ivac
  LOGICAL :: isfull

  isfull = .FALSE.
  ipart = 0
  ivac = 0

  !ALLOCATE(prt_vacantarr(prt_nnew))
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
