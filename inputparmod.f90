MODULE inputparmod

  USE kindmod
  IMPLICIT NONE

  ! Parameters
  INTEGER(iknd) :: in_nt, in_nr, in_ng
  INTEGER(iknd) :: seed
  INTEGER(iknd) :: in_Ns, in_Npartmax

  REAL(rknd) :: t_elapsed, Lr, alpha
  REAL(rknd) :: nidecay

  LOGICAL :: isvelocity, puretran

END MODULE inputparmod
