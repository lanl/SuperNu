MODULE gasgridmod

  USE kindmod
  IMPLICIT NONE

  INTEGER(iknd) :: gas_nr = 0
  INTEGER(iknd) :: gas_ng = 0
  INTEGER(iknd) :: gas_velno, gas_velyes

  REAL(rknd) :: gas_emat

  REAL(rknd), DIMENSION(:), ALLOCATABLE :: gas_rarr  !(gas_nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: gas_drarr, gas_dr3arr  !(gas_nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: gas_edep, gas_temp, gas_sigmap, gas_fcoef, gas_ur !(gas_nr)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: gas_tempb  !(gas_nr+1)
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: gas_rhoarr, gas_bcoef, gas_emit, gas_nisource  !(gas_nr)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: gas_sigmapg, gas_sigmargleft, gas_sigmargright, gas_emitprobg  !(gas_nr,gas_ng)
  REAL(rknd), DIMENSION(:,:), ALLOCATABLE :: gas_sigmal, gas_ppl, gas_sigmar, gas_ppr  !(gas_nr,gas_ng)

  CONTAINS

  SUBROUTINE gasgrid_init(nr,ng)
    INTEGER(iknd) :: nr, ng
    gas_nr = nr
    gas_ng = ng
  END SUBROUTINE gasgrid_init

END MODULE gasgridmod
