module timestepmod

  implicit none

  character(4) :: tsp_gridtype = 'none'  !total # of time steps
  integer :: tsp_nt = 0  !total # of time steps
  integer :: tsp_itrestart = 0 !restart time step # (at beginning of time step)
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
!-- read timestep configuration from file
    if(tsp_gridtype=='read') then
       call read_timestep_preset(tsp_nt)
       tsp_t = tsp_tpreset(1)
    else
!-- configured by input parameters
       tsp_t = tsp_tfirst
    endif
!!}}}
  end subroutine timestepmod_init


  subroutine timestep_update
    implicit none!{{{
!***********************************************************************
! update the timestep variables
!***********************************************************************
    real*8 :: help
!
!-- preset time step sizes
    select case(tsp_gridtype)
    case('read')
       if(.not.allocated(tsp_tpreset)) stop 'timestep_update: tpreset not allocated'
       tsp_dt = tsp_tpreset(tsp_it+1) - tsp_tpreset(tsp_it)
       tsp_t = tsp_tpreset(tsp_it)
    case('lin ')
!-- linear time grid
       tsp_dt = (tsp_tlast - tsp_tfirst)/tsp_nt
       tsp_t = tsp_tfirst + (tsp_it-1)*tsp_dt  !beginning of the time step
    case('expo')
!-- exponential time grid
       help = log(tsp_tlast/tsp_tfirst)/tsp_nt
       tsp_t = tsp_tfirst*exp((tsp_it-1)*help)  !beginning of the time step
       tsp_dt = tsp_tfirst*exp(tsp_it*help) - tsp_t
    end select

    tsp_tcenter = tsp_t + .5*tsp_dt
    tsp_dtinv = 1d0/tsp_dt
!}}}
  end subroutine timestep_update


  subroutine read_timestep_preset(nt)
    implicit none!{{{
    integer,intent(out) :: nt
!***********************************************************************
! read preset timestep values
!***********************************************************************
    integer :: istat
    character :: dmy
!
    open(4,file='input.tsp_time',status='old',iostat=istat)
    if(istat/=0) stop 'rd_tsp_preset: file missing: input.tsp_time'
!
!-- read header
    read(4,*,iostat=istat) dmy, nt
    if(istat/=0) stop 'rd_tsp_preset: file header error: input.tsp_time'
!
!-- allocate
    if(nt<1) stop 'rd_tsp_preset: header nt < 1'
    allocate(tsp_tpreset(nt))
!
!-- read lines
    read(4,*,iostat=istat) tsp_tpreset
    if(istat/=0) stop 'rd_tsp_preset: file error: input.tsp_time'
!
!-- make sure no remaining data
    read(4,*,iostat=istat) !-- last time step value is only for tsp_dt purpose
    if(istat/=0) stop 'rd_tsp_preset: file body too short'
    read(4,*,iostat=istat) !-- last time step value is only for tsp_dt purpose
    if(istat==0) stop 'rd_tsp_preset: file body too long'
!
    close(4)
!
    write(6,*) 'rd_tsp_preset: custom timesteps read successfully',nt
!!}}}
  end subroutine read_timestep_preset

end module timestepmod
