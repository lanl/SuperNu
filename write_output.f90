SUBROUTINE write_output

  USE data_mod
  IMPLICIT NONE

  INTEGER(iknd) :: ir
  CHARACTER(16), SAVE :: pos='REWIND'

  OPEN(UNIT=2,FILE='time.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=3,FILE='temp.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=4,FILE='fcoef.dat',STATUS='UNKNOWN',POSITION=pos)
  OPEN(UNIT=7,FILE='Lum.dat',STATUS='UNKNOWN',POSITION=pos)
  
  WRITE(2,*) time
  WRITE(7,*) Eright/dt

  DO ir = 1, nr
     WRITE(3,'(es16.8)',ADVANCE='NO') Temp(ir)
     WRITE(4,'(es16.8)',ADVANCE='NO') fcoef(ir)
  ENDDO

  CLOSE(2)
  CLOSE(3)
  CLOSE(4)
  CLOSE(7)

  pos='APPEND'

END SUBROUTINE write_output
