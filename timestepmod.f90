MODULE timestepmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: time_nt, tn
  REAL(rknd) :: texp, time, dt

  CONTAINS

    SUBROUTINE timestep_init(nt)
      time_nt = nt
    END SUBROUTINE timestep_init
    
END MODULE timestepmod
