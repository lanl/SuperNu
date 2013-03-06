MODULE particlemod

  USE kindmod
  IMPLICIT NONE

  TYPE packet
     INTEGER(iknd) :: zsrc, gsrc, rtsrc
     REAL(rknd) :: rsrc, musrc, tsrc
     REAL(rknd) :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet
  TYPE(packet), DIMENSION(:), POINTER :: prt_particles

  INTEGER(iknd) :: prt_npartmax, prt_ns
  INTEGER(iknd) :: prt_nsurf, prt_nnew
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: prt_numcensus, prt_nvol, prt_vacantarr

  REAL(rknd) :: prt_eleft, prt_eright, prt_erad, prt_einit, prt_einp, prt_eint, prt_etot, prt_esurf

  LOGICAL :: prt_done

  CONTAINS

    SUBROUTINE particle_init(npartmax,ns)
      INTEGER(iknd) :: npartmax, ns
      prt_npartmax = Npartmax
      prt_ns = Ns
    END SUBROUTINE particle_init

END MODULE particlemod
