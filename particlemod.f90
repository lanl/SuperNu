module particlemod

  implicit none

  TYPE packet
     integer :: zsrc, gsrc, rtsrc
     real*8 :: rsrc, musrc, tsrc
     real*8 :: Esrc, Ebirth
     logical :: isvacant
  end TYPE packet
  TYPE(packet), dimension(:), pointer :: prt_particles  !(prt_npartmax)

  integer :: prt_npartmax, prt_ns
  integer :: prt_nsurf, prt_nnew
  integer, dimension(:), allocatable :: prt_vacantarr

  logical :: prt_done

  save

  contains

  subroutine particle_init(npartmax,ns)
!--------------------------------------
    integer,intent(in) :: npartmax, ns
!***********************************************************************
! init particle module
!***********************************************************************
    integer :: ipart
!
!-- adopt input values in module internal storage
    prt_npartmax = npartmax
    prt_ns = ns
!
!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax))

!-- Setting all entries of particle array to vacant: loop
    do ipart = 1, prt_npartmax
       prt_particles(ipart)%isvacant=.true.
    enddo
  end subroutine particle_init

end module particlemod
