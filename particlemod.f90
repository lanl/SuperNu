module particlemod

  implicit none

  type packet
     real*8 :: x, y, z
     real*8 :: mu, om, t
     real*8 :: e, e0, wl
  end type packet
!
!-- secondary particle properties
  type packet2
     integer :: ix, iy, iz   !positional cell indices
     integer :: ic, ig       !index into compressed domain arrays, group index
     integer :: itype        !IMC or DDMC type
     integer :: ipart, istep !particle number and transport step number
     integer :: idist        !transport distance identifier
     logical :: done, isvacant !transport done and particle termination flags
  end type packet2
!
  type(packet),allocatable,target :: prt_particles(:)  !(prt_npartmax)
  logical,allocatable :: prt_isvacant(:)  !(prt_npartmax)
!
  integer :: prt_npartmax

  logical :: prt_isimcanlog !sets flux tally and energy deposition ...
  !to analog in IMC
  logical :: prt_isddmcanlog !sets flux tally and energy deposition ...
  !to analog in DDMC

  real*8 :: prt_tauddmc
  real*8 :: prt_taulump
  character(4) :: prt_tauvtime ! unif|incr
!
!
!-- rtw: random number counter added (rev. 262).
  integer :: prt_tlyrand
!-- rtw: array of rand counts from each rank
  integer, allocatable :: prt_tlyrandarr(:)
!-- particle property restart arrays:
   logical, allocatable :: prt_tlyvacant(:,:)
  integer, allocatable :: prt_tlyzsrc(:,:), prt_tlyrtsrc(:,:)
  real*8, allocatable :: prt_tlyrsrc(:,:), prt_tlymusrc(:,:), prt_tlytsrc(:,:)
  real*8, allocatable :: prt_tlyesrc(:,:), prt_tlyebirth(:,:), prt_tlywlsrc(:,:)

  save

  contains

  subroutine particlemod_init(npartmax,isimcanlog, &
       isddmcanlog,tauddmc,taulump,tauvtime)!{{{
!--------------------------------------
    implicit none
    integer,intent(in) :: npartmax
    logical,intent(in) :: isimcanlog, isddmcanlog
    real*8,intent(in) :: tauddmc, taulump
    character(4),intent(in) :: tauvtime
!***********************************************************************
! init particle module
!***********************************************************************
!-- adopt input values in module internal storage
    prt_npartmax = npartmax
    prt_isimcanlog = isimcanlog
    prt_isddmcanlog = isddmcanlog
    prt_tauddmc = tauddmc
    prt_taulump = taulump
    prt_tauvtime = tauvtime
!!}}}
  end subroutine particlemod_init


  subroutine particle_alloc(ltalk,norestart,nummespasint)
!--------------------------------------------------!{{{
    implicit none
    logical,intent(in) :: ltalk,norestart
    integer,intent(in) :: nummespasint

    integer :: n

!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax),prt_isvacant(prt_npartmax))
    prt_isvacant = .true.
!
!-- print size only on master rank
    if(ltalk) then
      n = int(sizeof(prt_particles)/1024) !kB
      write(6,*) 'ALLOC particles:',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif

!
!-- restart capabilities
    if(.not.norestart) then
!-- rand() count per rank allocation
       allocate(prt_tlyrandarr(nummespasint))
       prt_tlyrandarr = 0
!-- mpi gather arrays for particles
       allocate(prt_tlyvacant(prt_npartmax,nummespasint))
       allocate(prt_tlyzsrc(prt_npartmax,nummespasint))
       allocate(prt_tlyrtsrc(prt_npartmax,nummespasint))
       allocate(prt_tlyrsrc(prt_npartmax,nummespasint))
       allocate(prt_tlymusrc(prt_npartmax,nummespasint))
       allocate(prt_tlytsrc(prt_npartmax,nummespasint))
       allocate(prt_tlyesrc(prt_npartmax,nummespasint))
       allocate(prt_tlyebirth(prt_npartmax,nummespasint))
       allocate(prt_tlywlsrc(prt_npartmax,nummespasint))
    endif
!}}}
  end subroutine particle_alloc


  subroutine particle_dealloc
    deallocate(prt_particles,prt_isvacant)
!-- restart capabilities
    if(allocated(prt_tlyrandarr)) then
       deallocate(prt_tlyrandarr)
!-- mpi gather arrays for particles
       deallocate(prt_tlyvacant)
       deallocate(prt_tlyzsrc)
       deallocate(prt_tlyrtsrc)
       deallocate(prt_tlyrsrc)
       deallocate(prt_tlymusrc)
       deallocate(prt_tlytsrc)
       deallocate(prt_tlyesrc)
       deallocate(prt_tlyebirth)
       deallocate(prt_tlywlsrc)
    endif
  end subroutine particle_dealloc


  subroutine tau_update(t,tfirst,tlast)
    use physconstmod
    implicit none
    real*8,intent(in) :: t,tfirst,tlast
!-------------------------------------------------
! Update mfp thresholds for DDMC:
! - IMC-DDMC transition,
! - DDMC group lumping.
!-------------------------------------------------
    real*8 :: slp1,slp2

    if(prt_tauvtime=='unif') then
!-- constant thresholds
       return
    elseif(prt_tauvtime=='incr') then
!-- linear increase in thresholds (max mfp thresh = 5)
       slp1=0.667d0*prt_tauddmc/(pc_day*(tlast-tfirst))
       slp2=0.667d0*prt_taulump/(pc_day*(tlast-tfirst))

       prt_tauddmc = prt_tauddmc+(t-pc_day*tfirst)*slp1
       prt_taulump = prt_taulump+(t-pc_day*tfirst)*slp2
    else
       stop 'tau_update: prt_tauvtime invalid'
    endif


  end subroutine tau_update

end module particlemod
