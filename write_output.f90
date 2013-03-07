SUBROUTINE write_output

  USE timestepmod
  USE gasgridmod
  USE particlemod
  IMPLICIT NONE

  INTEGER :: ir
  CHARACTER(16), SAVE :: pos='REWIND'

  OPEN(UNIT=2,FILE='tsp_time.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=3,FILE='temp.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=4,FILE='gas_fcoef.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=7,FILE='Lum.dat',STATUS='UNKNOWN',POSITION=pos)
  
  WRITE(2,*) tsp_time
  WRITE(7,*) gas_eright/tsp_dt

  DO ir = 1, gas_nr
     WRITE(3,'(es16.8)',ADVANCE='NO') gas_temp(ir)
     WRITE(4,'(es16.8)',ADVANCE='NO') gas_fcoef(ir)
  ENDDO

  CLOSE(2)
  CLOSE(3)
  CLOSE(4)
  CLOSE(7)

  pos='APPEND'

END SUBROUTINE write_output
