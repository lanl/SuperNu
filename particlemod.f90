MODULE particlemod

  IMPLICIT NONE

  TYPE packet
     INTEGER :: zsrc, gsrc, rtsrc
     REAL*8 :: rsrc, musrc, tsrc
     REAL*8 :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet
  TYPE(packet), DIMENSION(:), POINTER :: prt_particles

  INTEGER :: prt_npartmax, prt_ns
  INTEGER :: prt_nsurf, prt_nnew
  INTEGER, DIMENSION(:), ALLOCATABLE :: prt_vacantarr

  REAL*8 :: prt_eleft, prt_eright, prt_erad, prt_einit, prt_einp, prt_eint, prt_etot, prt_esurf

  LOGICAL :: prt_done

  save

  CONTAINS

  SUBROUTINE particle_init(npartmax,ns)
!--------------------------------------
    INTEGER,intent(in) :: npartmax, ns
    prt_npartmax = npartmax
    prt_ns = ns
!
!-- allocate permanent storage (dealloc in dealloc_all.f)
    ALLOCATE(prt_particles(prt_npartmax))
  END SUBROUTINE particle_init

END MODULE particlemod
