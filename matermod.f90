MODULE matermod

  USE kindmod

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rarr, drarr, dr3arr
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Edep, Temp, sigmap, fcoef, Ur, Tempb
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rhoarr, bcoef, Emit, nisource
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmapg, sigmargleft, sigmargright, EmitProbg
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmaL, PPL, sigmaR, PPR

END MODULE matermod
