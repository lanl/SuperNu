MODULE timestepmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: tsp_nt, tn
  REAL(rknd) :: texp, time, dt

  CONTAINS

    SUBROUTINE timestep_init(nt)
      INTEGER(iknd) nt
      tsp_nt = nt
    END SUBROUTINE timestep_init
    
END MODULE timestepmod
