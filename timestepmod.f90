MODULE timestepmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: tsp_nt, tsp_tn
  REAL(rknd) :: tsp_texp, tsp_time, tsp_dt

  CONTAINS

    SUBROUTINE timestep_init(nt)
      INTEGER(iknd) nt
      tsp_nt = nt
    END SUBROUTINE timestep_init
    
END MODULE timestepmod
