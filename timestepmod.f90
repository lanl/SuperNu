module timestepmod

  implicit none

  integer :: tsp_nt = 0  !total # of time steps
  integer :: tsp_tn  !current time step number
  real*8 :: tsp_texp
  real*8 :: tsp_tcenter
  real*8 :: tsp_time
  real*8 :: tsp_dt
  real*8 :: tsp_alpha = 0d0

  save

  contains


  subroutine timestep_init(nt, alpha, tfirst, dt)
!------------------------------------------------
    use physconstmod
    integer,intent(in) :: nt
    real*8,intent(in) :: alpha, tfirst, dt
!***********************************************************************
! set the timestep constants
!***********************************************************************
    tsp_nt = nt
    tsp_dt = dt
    tsp_alpha = alpha

!-- beginning of first time step
    tsp_time = 0d0
    tsp_texp = tfirst*pc_day
    tsp_tcenter = tsp_texp + .5d0*tsp_dt
  end subroutine timestep_init


  subroutine timestep_update(dt)
    implicit none
    real*8,intent(in) :: dt
!***********************************************************************
! update the timestep variables
!***********************************************************************
    tsp_dt = dt
    tsp_time = tsp_time+tsp_dt
    tsp_texp = tsp_texp+tsp_dt
    tsp_tcenter = tsp_texp + .5*tsp_dt
  end subroutine timestep_update
    
end module timestepmod
