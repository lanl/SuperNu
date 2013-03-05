MODULE gasgridmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: gas_nr = 0
  INTEGER(iknd) :: gas_ng = 0
  INTEGER(iknd) :: velno, velyes

  REAL(rknd) :: Emat

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rarr  !(gas_nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: drarr, dr3arr  !(gas_nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Edep, Temp, sigmap, fcoef, Ur !(gas_nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Tempb  !(gas_nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: rhoarr, bcoef, Emit, nisource  !(gas_nr)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmapg, sigmargleft, sigmargright, EmitProbg  !(gas_nr,gas_ng)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: sigmaL, PPL, sigmaR, PPR  !(gas_nr,gas_ng)

  CONTAINS

  SUBROUTINE gasgrid_init(nr,ng)
    gas_nr = nr
    gas_ng = ng
  END SUBROUTINE gasgrid_init

END MODULE gasgridmod
