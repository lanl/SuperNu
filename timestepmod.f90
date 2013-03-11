MODULE timestepmod

  IMPLICIT NONE

  INTEGER :: tsp_nt = 0  !total # of time steps
  INTEGER :: tsp_tn  !current time step number
  REAL*8 :: tsp_texp
  REAL*8 :: tsp_tcenter
  REAL*8 :: tsp_time
  REAL*8 :: tsp_dt
  REAL*8 :: tsp_alpha = 0d0

  save

  CONTAINS


  SUBROUTINE timestep_init(nt, alpha, tfirst, dt)
!------------------------------------------------
    use physconstmod
    INTEGER,intent(in) :: nt
    REAL*8,intent(in) :: alpha, tfirst, dt
!***********************************************************************
! set the timestep constants
!***********************************************************************
    tsp_nt = nt
    tsp_dt = dt
    tsp_alpha = alpha

!-- beginning of first time step
    tsp_time = 0d0
    !tsp_texp = tfirst*pc_day
    tsp_texp = 0.14
    tsp_tcenter = tsp_texp + .5d0*tsp_dt
  END SUBROUTINE timestep_init


  subroutine timestep_update(dt)
    implicit none
    REAL*8,intent(in) :: dt
!***********************************************************************
! update the timestep variables
!***********************************************************************
    tsp_dt = dt
    tsp_time = tsp_time+tsp_dt
    tsp_texp = tsp_texp+tsp_dt
    tsp_tcenter = tsp_texp + .5*tsp_dt
  end subroutine timestep_update
    
END MODULE timestepmod
