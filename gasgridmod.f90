MODULE gasgridmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: gas_nr = 0
  INTEGER(iknd) :: gas_ng = 0
  INTEGER(iknd) :: velno, velyes

  REAL(rknd) :: Emat

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rarr  !(nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: drarr, dr3arr  !(nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Edep, Temp, sigmap, fcoef, Ur !(nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Tempb  !(nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rhoarr, bcoef, Emit, nisource  !(nr)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmapg, sigmargleft, sigmargright, EmitProbg  !(nr,ng)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmaL, PPL, sigmaR, PPR  !(nr,ng)

  CONTAINS

  SUBROUTINE gasgrid_init(nr,ng)
    gas_nr = nr
    gas_ng = ng
  END SUBROUTINE gasgrid_init

END MODULE gasgridmod
