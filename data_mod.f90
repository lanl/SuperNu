MODULE data_mod
  INTEGER, PARAMETER :: rknd = SELECTED_REAL_KIND(14,100)
  INTEGER, PARAMETER :: iknd = SELECTED_INT_KIND(8)

  ! Constants
  REAL(rknd) :: pi = 4.0*ATAN(1.0)
  REAL(rknd) :: lspeed = 2.998e10  !light speed (cm/s)
  REAL(rknd) :: a_coef = 1.371e14  !radiation constant (erg/Kev^4/cm^3)
  REAL(rknd) :: Nav = 6.022e23  ! Avogadro's number

  ! Parameters
  INTEGER(iknd) :: velno, velyes
  INTEGER(iknd) :: nt, nr, ng, advoption, tn
  INTEGER(iknd) :: Npartmax, Ns, seed
  INTEGER(iknd) :: Nsurf, Nnew
  INTEGER(iknd), DIMENSION(:), ALLOCATABLE :: numcensus, Nvol, vacantarr

  REAL(rknd) :: texp, time, dt, Eleft, Eright, Erad, Einit, Einp, Eint, Etot, Esurf, Emat
  REAL(rknd) :: t_elapsed, Lr, alpha, nidecay
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rarr, drarr, dr3arr
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Edep, Temp, sigmap, fcoef, Ur, Tempb
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rhoarr, bcoef, Emit, nisource
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmapg, sigmargleft, sigmargright, EmitProbg
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmaL, PPL, sigmaR, PPR

  LOGICAL :: done, isvelocity, puretran

  ! Derived data types
  TYPE packet
     INTEGER(iknd) :: zsrc, gsrc, rtsrc
     REAL(rknd) :: rsrc, musrc, tsrc
     REAL(rknd) :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet
  TYPE(packet), DIMENSION(:), POINTER :: particles

  ! Picket-fence probabilities
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Ppick

END MODULE data_mod
