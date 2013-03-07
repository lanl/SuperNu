SUBROUTINE grids

  USE gasgridmod
  USE timestepmod
  IMPLICIT NONE

  ! This subroutine creates the spatial grid (or velocity grid if gas_velyes=1)

  INTEGER(iknd) :: ir

  tsp_texp = 0.14
  gas_rarr(1) = 0.0_rknd
  DO ir = 1, gas_nr
     gas_drarr(ir)=gas_lr/REAL(gas_nr)
     gas_rarr(ir+1)=gas_rarr(ir)+gas_drarr(ir)
     gas_dr3arr(ir)=gas_rarr(ir+1)**3-gas_rarr(ir)**3
  ENDDO

  ! r/tsp_texp = velocity grid (calculated with initial spatial grid and 
  ! initial expansion tsp_time)

  IF (gas_isvelocity.EQV..TRUE.) THEN
     gas_rarr = gas_rarr/tsp_texp
     gas_drarr = gas_drarr/tsp_texp
     gas_dr3arr = gas_dr3arr/tsp_texp**3
  ENDIF


END SUBROUTINE grids
