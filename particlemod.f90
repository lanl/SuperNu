MODULE particlemod

  USE kindmod
  IMPLICIT NONE

  TYPE packet
     INTEGER(iknd) :: zsrc, gsrc, rtsrc
     REAL(rknd) :: rsrc, musrc, tsrc
     REAL(rknd) :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet
  TYPE(packet), DIMENSION(:), POINTER :: particles

  INTEGER(iknd) :: prt_Npartmax, prt_Ns
  INTEGER(iknd) :: Nsurf, Nnew
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: numcensus, Nvol, vacantarr

  REAL(rknd) :: Eleft, Eright, Erad, Einit, Einp, Eint, Etot, Esurf

  LOGICAL :: done

  CONTAINS

    SUBROUTINE particle_init(Npartmax,Ns)
      prt_Npartmax = Npartmax
      prt_Ns = Ns
    END SUBROUTINE particle_init

END MODULE particlemod
