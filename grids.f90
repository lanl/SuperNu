SUBROUTINE grids
  USE data_mod
  IMPLICIT NONE

  ! This subroutine creates the spatial grid (or velocity grid if velyes=1)

  INTEGER(iknd) :: ir

  texp = 0.14
  rarr(1) = 0.0_rknd
  DO ir = 1, nr
     drarr(ir)=Lr/REAL(nr)
     rarr(ir+1)=rarr(ir)+drarr(ir)
     dr3arr(ir)=rarr(ir+1)**3-rarr(ir)**3
  ENDDO

  ! r/texp = velocity grid (calculated with initial spatial grid and 
  ! initial expansion time)

  IF (isvelocity.EQV..TRUE.) THEN
     rarr = rarr/texp
     drarr = drarr/texp
     dr3arr = dr3arr/texp**3
  ENDIF


END SUBROUTINE grids
