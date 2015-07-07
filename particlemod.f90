module particlemod

  implicit none

  type packet
     real*8 :: x, y, z
     real*8 :: mu, om, t
     real*8 :: e, e0, wl
     integer :: icorig
  end type packet
!
!-- secondary particle properties
  type packet2
     real*8 :: mux,muy,muz   !particle direction in lab frame cartesian coordinates
     real*8 :: dist          !particle travel distance
     integer :: ix, iy, iz   !positional cell indices
     integer :: ic, ig       !index into compressed domain arrays, group index
     integer :: itype        !IMC or DDMC type
     integer :: ipart, istep !particle number and transport step number
     integer :: idist        !transport distance identifier
     logical :: done, lflux, lcens, isvacant !transport done and particle termination flags
  end type packet2
!
  type(packet),allocatable,target :: prt_particles(:)  !(prt_npartmax)
  logical,allocatable :: prt_isvacant(:)  !(prt_npartmax)
!
  integer :: prt_npartmax

  logical :: prt_isimcanlog  !sets flux tally and energy deposition ...
  !to analog in IMC
  logical :: prt_isddmcanlog !sets flux tally and energy deposition ...
  !to analog in DDMC

  real*8 :: prt_tauddmc
  real*8 :: prt_taulump
  character(4) :: prt_tauvtime ! unif|incr

  save

  contains

  subroutine particle_alloc(ltalk)
!--------------------------------------------------!{{{
    implicit none
    logical,intent(in) :: ltalk

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
!-- output
    if(ltalk) then
       write(6,*)
       write(6,*) 'particle array:'
       write(6,*) '===================='
       write(6,*) 'npart :',prt_npartmax
       write(6,*)
    endif
!}}}
  end subroutine particle_alloc


  subroutine particle_dealloc
    implicit none
    deallocate(prt_particles,prt_isvacant)
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
       slp1=0.667d0*prt_tauddmc/(tlast-tfirst)
       slp2=0.667d0*prt_taulump/(tlast-tfirst)

       prt_tauddmc = prt_tauddmc+(t-tfirst)*slp1
       prt_taulump = prt_taulump+(t-tfirst)*slp2
    else
       stop 'tau_update: prt_tauvtime invalid'
    endif
  end subroutine tau_update

end module particlemod
