module timestepmod

  implicit none

  integer :: tsp_nt = 0  !total # of time steps
  integer :: tsp_ntres = 0 !restart time step # (at beginning of time step)
  integer :: tsp_it  !current time step number
  real*8 :: tsp_t
  real*8 :: tsp_tcenter
  real*8 :: tsp_dt
  real*8 :: tsp_alpha = 0d0

  save

  contains


  subroutine timestep_init(nt, ntres, alpha, tfirst, dt, isbdf2)
!------------------------------------------------
    use physconstmod
    integer,intent(in) :: nt, ntres
    real*8,intent(in) :: alpha, tfirst, dt
    logical,intent(in) :: isbdf2
!***********************************************************************
! set the timestep constants
!***********************************************************************
    if(isbdf2) then
       tsp_nt = nt+1
       tsp_dt = dt/2d0
    else
       tsp_nt = nt
       tsp_dt = dt
    endif
    if(ntres<1) then
       tsp_ntres=1
    else
       tsp_ntres=ntres
    endif

    tsp_alpha = alpha

!-- beginning of first (restart) time step
    tsp_t = tfirst*pc_day+(tsp_ntres-1)*tsp_dt
    tsp_tcenter = tsp_t + .5d0*tsp_dt
!
  end subroutine timestep_init


  subroutine timestep_update(dt, isbdf2)
    implicit none
    real*8,intent(in) :: dt
    logical,intent(in) :: isbdf2
!***********************************************************************
! update the timestep variables
!***********************************************************************
    tsp_dt = dt
    if(.not.isbdf2.or.(isbdf2.and.tsp_it/=1)) then
       tsp_t = tsp_t+tsp_dt
       tsp_tcenter = tsp_t + .5*tsp_dt
    endif
  end subroutine timestep_update

end module timestepmod
