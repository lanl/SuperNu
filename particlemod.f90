MODULE particlemod

  IMPLICIT NONE

  TYPE packet
     INTEGER :: zsrc, gsrc, rtsrc
     REAL*8 :: rsrc, musrc, tsrc
     REAL*8 :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet
  TYPE(packet), DIMENSION(:), POINTER :: prt_particles  !(prt_npartmax)

  INTEGER :: prt_npartmax, prt_ns
  INTEGER :: prt_nsurf, prt_nnew
  INTEGER, DIMENSION(:), ALLOCATABLE :: prt_vacantarr

  LOGICAL :: prt_done

  save

  CONTAINS

  SUBROUTINE particle_init(npartmax,ns)
!--------------------------------------
    INTEGER,intent(in) :: npartmax, ns
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
    ALLOCATE(prt_particles(prt_npartmax))

!-- Setting all entries of particle array to vacant: loop
    DO ipart = 1, prt_npartmax
       prt_particles(ipart)%isvacant=.TRUE.
    ENDDO
  END SUBROUTINE particle_init

END MODULE particlemod
