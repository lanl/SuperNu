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


end module particlemod
