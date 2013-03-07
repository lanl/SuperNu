MODULE timestepmod

  IMPLICIT NONE

  INTEGER :: tsp_nt
  INTEGER :: tsp_tn
  REAL*8 :: tsp_texp
  REAL*8 :: tsp_time 
  REAL*8 :: tsp_dt
  REAL*8 :: tsp_alpha

  save

  CONTAINS

  SUBROUTINE timestep_init(nt, alpha)
!-----------------------------
    INTEGER,intent(in) :: nt
    REAL*8,intent(in) :: alpha
    tsp_nt = nt
    tsp_alpha = alpha
  END SUBROUTINE timestep_init
    
END MODULE timestepmod
