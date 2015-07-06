module timestepmod

  implicit none

  integer :: tsp_nt = 0  !total # of time steps
  integer :: tsp_ntres = 0 !restart time step # (at beginning of time step)
  integer :: tsp_it  !current time step number
  real*8 :: tsp_t
  real*8,allocatable :: tsp_tpreset(:)  !store preset time steps from input.tsp_time
  real*8 :: tsp_tcenter,tsp_tfirst,tsp_tlast
  real*8 :: tsp_dt,tsp_dtinv

  private read_timestep_preset

  save

  contains


  subroutine timestepmod_init
!----------------------------!{{{
    use physconstmod
!***********************************************************************
! set the timestep constants
!***********************************************************************
    tsp_ntres = max(tsp_ntres,1)
!
!-- read timestep configuration from file
    if(tsp_nt<0) then
       call read_timestep_preset(tsp_nt)
       tsp_t = tsp_tpreset(1)
    else
!-- configured by input parameters
       tsp_t = tsp_tfirst
    endif
!!}}}
  end subroutine timestepmod_init


  subroutine timestep_update(dt)
    implicit none!{{{
    real*8,intent(in) :: dt
!***********************************************************************
! update the timestep variables
!***********************************************************************
!-- preset time step sizes
    if(allocated(tsp_tpreset)) then
       tsp_dt = tsp_tpreset(tsp_it+1) - tsp_tpreset(tsp_it)
       tsp_dtinv = 1d0/tsp_dt
       tsp_t = tsp_tpreset(tsp_it)
       tsp_tcenter = .5*(tsp_tpreset(tsp_it+1) + tsp_tpreset(tsp_it))
       return
    endif
!
!-- update time step size
    tsp_dt = dt
    tsp_dtinv = 1d0/dt
!
!-- update time
    if(tsp_it==tsp_ntres) then
!-- first time step is exception
       tsp_t = tsp_t + max(tsp_ntres-1,0)*tsp_dt
    else
       tsp_t = tsp_t + tsp_dt
    endif
    tsp_tcenter = tsp_t + .5*tsp_dt
!}}}
  end subroutine timestep_update


  subroutine read_timestep_preset(nt)
    implicit none!{{{
    integer,intent(out) :: nt
!***********************************************************************
! read preset timestep values
!***********************************************************************
    integer :: istat
!
    open(4,file='input.tsp_time',status='old',iostat=istat)
    if(istat/=0) stop 'rd_tsp_preset: file missing: input.tsp_time'
!
!-- count lines
    nt = 0
    do
       read(4,*,end=9)
       nt = nt+1
    enddo
9   continue
    rewind(4)
!
!-- allocate
    if(nt<2) stop 'rd_tsp_preset: file too short: input.tsp_time'
    allocate(tsp_tpreset(nt))
!-- last time step is only for tsp_dt purpose
    nt = nt - 1
!
!-- read lines
    read(4,*,iostat=istat) tsp_tpreset
    if(istat/=0) stop 'rd_tsp_preset: file error: input.tsp_time'
!
!-- make sure no remaining data
    read(4,*,iostat=istat)
    if(istat==0) stop 'rd_tsp_preset: file too long: input.tsp_time'
!
    close(4)
!
    write(6,*) 'rd_tsp_preset: custom timesteps read successfully',nt
!!}}}
  end subroutine read_timestep_preset

end module timestepmod
