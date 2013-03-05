MODULE simparamod

  USE kindmod
  IMPLICIT NONE

  ! Parameters
  INTEGER(iknd) :: velno, velyes
  INTEGER(iknd) :: nt, nr, ng, advoption, tn
  INTEGER(iknd) :: Npartmax, Ns, seed
  INTEGER(iknd) :: Nsurf, Nnew
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: numcensus, Nvol, vacantarr

  REAL(rknd) :: texp, time, dt, Eleft, Eright, Erad, Einit, Einp, Eint, Etot, Esurf, Emat
  REAL(rknd) :: t_elapsed, Lr, alpha, nidecay

  LOGICAL :: done, isvelocity, puretran

  TYPE(packet), DIMENSION(:), POINTER :: particles

END MODULE simparamod
