module particlemod

  implicit none

  !Ryan W.: Changing group attribute to continuous wavelength (rev. 120)
  type packet
     real*8 :: rsrc, musrc, tsrc
     real*8 :: esrc, ebirth, wlsrc
     integer :: zsrc, rtsrc !,gsrc
     logical :: isvacant,ldummy
  end type packet
  type(packet),allocatable,target :: prt_particles(:)  !(prt_npartmax)
!
  integer :: prt_npartmax, prt_ns, prt_ninit
  integer :: prt_nsurf, prt_nexsrc, prt_nnew, prt_ninitnew
!-- rtw: random number counter added (rev. 262). associated with particle routines
  integer :: prt_tlyrand
!-- rtw: array of rand counts from each rank
  integer, allocatable :: prt_tlyrandarr(:)
!-- particle property restart arrays:
   logical, allocatable :: prt_tlyvacant(:,:)
  integer, allocatable :: prt_tlyzsrc(:,:), prt_tlyrtsrc(:,:)
  real*8, allocatable :: prt_tlyrsrc(:,:), prt_tlymusrc(:,:), prt_tlytsrc(:,:)
  real*8, allocatable :: prt_tlyesrc(:,:), prt_tlyebirth(:,:), prt_tlywlsrc(:,:)
!
  integer, allocatable :: prt_vacantarr(:) !array of vacant particle array locations

  logical :: prt_done
  logical :: prt_isimcanlog !sets flux tally and energy deposition ...
  !to analog in IMC
  logical :: prt_isddmcanlog !sets flux tally and energy deposition ...
  !to analog in DDMC

  real*8 :: prt_tauddmc
  real*8 :: prt_taulump
  character(4) :: prt_tauvtime ! unif|incr

  save

  contains

  subroutine particle_init(npartmax,ns,ninit,isimcanlog, &
       isddmcanlog,tauddmc,taulump,tauvtime,nummespasint,norestart)
!--------------------------------------
    integer,intent(in) :: npartmax, ns, ninit, nummespasint
    logical,intent(in) :: isimcanlog, isddmcanlog,norestart
    real*8,intent(in) :: tauddmc, taulump
    character(4),intent(in) :: tauvtime
!***********************************************************************
! init particle module
!***********************************************************************
!
!-- adopt input values in module internal storage
    prt_npartmax = npartmax
    prt_ns = ns
    prt_ninit = ninit
    prt_isimcanlog = isimcanlog
    prt_isddmcanlog = isddmcanlog
    prt_tauddmc = tauddmc
    prt_taulump = taulump
    prt_tauvtime = tauvtime
!
!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax))
    write(6,*) 'allocate prt_particles:',&
       sizeof(prt_particles)/1024**2,"MB"
    prt_particles%isvacant = .true.
    if(.not.norestart) then
!-- rand() count per rank allocation
       allocate(prt_tlyrandarr(nummespasint))
       prt_tlyrandarr = 0
!-- mpi gather arrays for particles
       allocate(prt_tlyvacant(npartmax,nummespasint))
       allocate(prt_tlyzsrc(npartmax,nummespasint))
       allocate(prt_tlyrtsrc(npartmax,nummespasint))
       allocate(prt_tlyrsrc(npartmax,nummespasint))
       allocate(prt_tlymusrc(npartmax,nummespasint))
       allocate(prt_tlytsrc(npartmax,nummespasint))
       allocate(prt_tlyesrc(npartmax,nummespasint))
       allocate(prt_tlyebirth(npartmax,nummespasint))
       allocate(prt_tlywlsrc(npartmax,nummespasint))
    endif
!
  end subroutine particle_init

end module particlemod
