SUBROUTINE grids

  USE gasgridmod
  USE timestepmod
  IMPLICIT NONE

!##################################################
  ! This subroutine creates the spatial grid (or velocity grid if gas_velyes=1)
!##################################################

  INTEGER :: ir

  !Initial expansion time (to be set to start time)
  tsp_texp = 0.14
  !Initial inner most radius
  gas_rarr(1) = 0.0d0
  ! Initial grid, cell length, and cell volume generation loop
  DO ir = 1, gas_nr
     gas_drarr(ir)=gas_lr/REAL(gas_nr)
     gas_rarr(ir+1)=gas_rarr(ir)+gas_drarr(ir)
     gas_vals2(ir)%dr3_34pi=gas_rarr(ir+1)**3-gas_rarr(ir)**3
  ENDDO

  ! r/tsp_texp = velocity grid (calculated with initial spatial grid and 
  ! initial expansion tsp_time)

  IF (gas_isvelocity.EQV..TRUE.) THEN
     gas_rarr = gas_rarr/tsp_texp
     gas_drarr = gas_drarr/tsp_texp
     DO ir = 1, gas_nr
        gas_vals2(ir)%dr3_34pi=gas_vals2(ir)%dr3_34pi/tsp_texp**3
     ENDDO
  ENDIF


END SUBROUTINE grids
