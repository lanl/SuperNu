MODULE constmod

  USE kindmod
  IMPLICIT NONE

  REAL(rknd) :: pi = 4.0*ATAN(1.0)
  REAL(rknd) :: lspeed = 2.998e10  !light speed (cm/s)
  REAL(rknd) :: a_coef = 1.371e14  !radiation constant (erg/Kev^4/cm^3)
  REAL(rknd) :: Nav = 6.022e23  ! Avogadro's number

  ! Picket-fence probabilities
  REAL(rknd), DIMENSION(:), ALLOCATABLE :: Ppick


END MODULE constmod
