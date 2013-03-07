MODULE timestepmod

  IMPLICIT NONE

  INTEGER :: tsp_nt, tsp_tn
  REAL*8 :: tsp_texp, tsp_time, tsp_dt

  save

  CONTAINS

  SUBROUTINE timestep_init(nt)
!-----------------------------
    INTEGER,intent(in) :: nt
    tsp_nt = nt
  END SUBROUTINE timestep_init
    
END MODULE timestepmod
