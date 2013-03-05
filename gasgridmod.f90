MODULE gasgridmod

  USE kindmod

  INTEGER(iknd) :: nr, ng
  INTEGER(iknd) :: velno, velyes

  REAL(rknd) :: Emat

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rarr  !(nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: drarr, dr3arr  !(nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Edep, Temp, sigmap, fcoef, Ur !(nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Tempb  !(nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rhoarr, bcoef, Emit, nisource  !(nr)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmapg, sigmargleft, sigmargright, EmitProbg  !(nr,ng)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmaL, PPL, sigmaR, PPR  !(nr,ng)

END MODULE gasgridmod
