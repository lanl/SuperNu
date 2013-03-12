SUBROUTINE write_output

  USE timestepmod
  USE gasgridmod
  USE particlemod
  IMPLICIT NONE

  INTEGER :: ir
  CHARACTER(16), SAVE :: pos='REWIND'

  OPEN(UNIT=4,FILE='output.tsp_time',STATUS='UNKNOWN',POSITION=pos)
  WRITE(4,*) tsp_time
  CLOSE(4)

  OPEN(UNIT=4,FILE='output.Lum',STATUS='UNKNOWN',POSITION=pos)
  WRITE(4,*) gas_eright/tsp_dt
  CLOSE(4)

  OPEN(UNIT=4,FILE='output.temp',STATUS='UNKNOWN',POSITION=pos)
  DO ir = 1, gas_nr
    WRITE(4,'(es16.8)',ADVANCE='NO') gas_vals2(ir)%tempkev
  ENDDO
  CLOSE(4)

  OPEN(UNIT=4,FILE='output.gas_fcoef',STATUS='UNKNOWN',POSITION=pos)
  DO ir = 1, gas_nr
    WRITE(4,'(es16.8)',ADVANCE='NO') gas_fcoef(ir)
  ENDDO
  CLOSE(4)

  pos='APPEND'

END SUBROUTINE write_output
