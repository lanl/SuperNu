!> @file rprocmod.f90
!!
!! module 'rprocmod' for calculating fitting formulae for the r-process
!!

!*****************************************************************************
!                                                                            *
!  Module 'hratelib' for computing heating rates in expanding ejecta         *
!  OK 29.08.2019                                                             *
!  OK 05.06.2022: updated with the new formula (as in the paper)             *
!  OK 25.04.2024: updated to use interpolation by epsdot                     *
!                                                                            *
!*****************************************************************************
module rprocmod
implicit none

! grid of velocity and Ye from which approximant is interpolated
DOUBLE PRECISION, PARAMETER :: &
   YE_GRID(10)= (/  .05d0,  .10d0,  .15d0,  .20d0,  .25d0,       &
                    .30d0,  .35d0,  .40d0,  .45d0,  .50d0 /),    &
     V_GRID(6)= (/  .05d0,  .1d0,  .2d0,  .3d0,  .4d0,  .5d0 /)

! approximant coefficients on the grid
!!   E0_GRID= RESHAPE( (/  10.0,   10.0,   10.0, 10.0, 10.0, 10.0, &
!!                         10.0,   10.0,   11.0, 11.0, 11.0, 11.0, &
!!                         14.0,   10.0,   11.0, 11.0, 11.0, 11.0, &
!!                         14.0,   10.0,   10.0, 10.0, 11.0, 11.0, &
!!                         20.0,   25.0,   40.0, 38.0, 58.0, 70.0, &
!!                         6.1,    18.0,   47.1, 47.1, 74.8, 74.8, &
!!                         7.3,    7.0,    16.3, 23.2, 43.2, 150., &
!!                         3.2e-3, 3.2e-3, 8e-3, 7e-3, 9e-3, 0.015,&
!!                         0.20,   0.20,   0.60, 1.50, 1.50, 1.50, &
!!                         0.40,   1.00,   2.00, 3.00, 3.00, 3.00 /),
DOUBLE PRECISION, PARAMETER, DIMENSION(SIZE(V_GRID), SIZE(YE_GRID)) :: &
   E0_GRID= RESHAPE( &
            (/  1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, 1.000d0, &
                1.000d0, 1.000d0, 1.041d0, 1.041d0, 1.041d0, 1.041d0, &
                1.146d0, 1.000d0, 1.041d0, 1.041d0, 1.041d0, 1.041d0, &
                1.146d0, 1.000d0, 1.000d0, 1.000d0, 1.041d0, 1.041d0, &
                1.301d0, 1.398d0, 1.602d0, 1.580d0, 1.763d0, 1.845d0, &
                0.785d0, 1.255d0, 1.673d0, 1.673d0, 1.874d0, 1.874d0, &
                0.863d0, 0.845d0, 1.212d0, 1.365d0, 1.635d0, 2.176d0, &
               -2.495d0,-2.495d0,-2.097d0,-2.155d0,-2.046d0,-1.824d0, &
               -0.699d0,-0.699d0,-0.222d0, 0.176d0, 0.176d0, 0.176d0, &
               -0.398d0, 0.000d0, 0.301d0, 0.477d0, 0.477d0, 0.477d0 /),&
               (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   ALP_GRID= RESHAPE( &
             (/ 1.37d0, 1.38d0, 1.41d0, 1.41d0, 1.41d0, 1.41d0,  &
                1.41d0, 1.38d0, 1.37d0, 1.37d0, 1.37d0, 1.37d0,  &
                1.41d0, 1.38d0, 1.37d0, 1.37d0, 1.37d0, 1.37d0,  &
                1.36d0, 1.25d0, 1.32d0, 1.32d0, 1.34d0, 1.34d0,  &
                1.44d0, 1.40d0, 1.46d0, 1.66d0, 1.60d0, 1.60d0,  &
                1.36d0, 1.33d0, 1.33d0, 1.33d0, 1.374d0,1.374d0, &
                1.40d0, 1.358d0,1.384d0,1.384d0,1.384d0,1.344d0, &
                1.80d0, 1.80d0, 2.10d0, 2.10d0, 1.90d0, 1.90d0, &
                8.00d0, 8.00d0, 7.00d0, 7.00d0, 7.00d0, 7.00d0, &
                1.40d0, 1.40d0, 1.40d0, 1.60d0, 1.60d0, 1.60d0 /), &
                (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   T0_GRID= RESHAPE( &
            (/  1.80d0, 1.40d0, 1.20d0, 1.20d0, 1.20d0,  1.20d0,  &
                1.40d0, 1.00d0, 0.85d0, 0.85d0, 0.85d0,  0.85d0,  &
                1.00d0, 0.80d0, 0.65d0, 0.65d0, 0.61d0,  0.61d0,  &
                0.85d0, 0.60d0, 0.45d0, 0.45d0, 0.45d0,  0.45d0,  &
                0.65d0, 0.38d0, 0.22d0, 0.18d0, 0.12d0,  0.095d0, &
                0.540d0,0.31d0, 0.18d0, 0.13d0, 0.095d0, 0.081d0, &
                0.385d0,0.235d0,0.1d0,  0.06d0, 0.035d0, 0.025d0, &
                26.0d0, 26.0d0, 0.4d0,  0.4d0,  0.12d0, -20.0d0,  &
                0.20d0, 0.12d0, 0.05d0, 0.03d0, 0.025d0, 0.021d0, &
                0.16d0, 0.08d0, 0.04d0, 0.02d0, 0.018d0, 0.016d0 /), &
             (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   SIG_GRID=RESHAPE( &
            (/  0.08d0, 0.08d0, 0.095d0,0.095d0,0.095d0,0.095d0, &
                0.10d0, 0.08d0, 0.070d0,0.070d0,0.070d0,0.070d0, &
                0.07d0, 0.08d0, 0.070d0,0.065d0,0.070d0,0.070d0, &
                0.040d0,0.030d0,0.05d0, 0.05d0, 0.05d0, 0.050d0, &
                0.05d0, 0.030d0,0.025d0,0.045d0,0.05d0, 0.05d0,  &
                0.11d0, 0.04d0, 0.021d0,0.021d0,0.017d0,0.017d0, &
                0.10d0, 0.094d0,0.068d0,0.05d0, 0.03d0, 0.01d0, &
                45.0d0, 45.0d0, 45.0d0, 45.0d0, 25.0d0, 40.0d0, &
                0.20d0, 0.12d0, 0.05d0, 0.03d0, 0.025d0,0.021d0, &
                0.03d0, 0.015d0,0.007d0,0.01d0, 0.009d0,0.007d0 /), &
             (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   ALP1_GRID=RESHAPE( &
             (/  7.50d0, 7.50d0, 7.50d0, 7.50d0, 7.50d0, 7.50d0, &
                 9.00d0, 9.00d0, 7.50d0, 7.50d0, 7.00d0, 7.00d0, &
                 8.00d0, 8.00d0, 7.50d0, 7.50d0, 7.00d0, 7.00d0, &
                 8.00d0, 8.00d0, 7.50d0, 7.50d0, 7.00d0, 7.00d0, &
                 8.00d0, 8.00d0, 5.00d0, 7.50d0, 7.00d0, 6.50d0, &
                 4.5d0,  3.8d0,  4.0d0,  4.0d0,  4.0d0,  4.0d0,  &
                 2.4d0,  3.8d0,  3.8d0,  3.21d0, 2.91d0, 3.61d0, &
                -1.55d0,-1.55d0,-0.75d0,-0.75d0,-2.50d0,-5.00d0, &
                -1.55d0,-1.55d0,-1.55d0,-1.55d0,-1.55d0,-1.55d0, &
                 3.00d0, 3.00d0, 3.00d0, 3.00d0, 3.00d0, 3.00d0 /),&
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   T1_GRID= RESHAPE( &
            (/ 0.040d0, 0.025d0, 0.014d0, 0.010d0, 0.008d0, 0.006d0, &
               0.040d0, 0.035d0, 0.020d0, 0.012d0, 0.010d0, 0.008d0, &
               0.080d0, 0.040d0, 0.020d0, 0.012d0, 0.012d0, 0.009d0, &
               0.080d0, 0.040d0, 0.030d0, 0.018d0, 0.012d0, 0.009d0, &
               0.080d0, 0.060d0, 0.065d0, 0.028d0, 0.020d0, 0.015d0, &
                0.14d0, 0.123d0, 0.089d0, 0.060d0, 0.045d0, 0.031d0, &
               0.264d0,   0.1d0,  0.07d0, 0.055d0, 0.042d0, 0.033d0, &
                 1.0d0,   1.0d0,   1.0d0,   1.0d0,  0.02d0,  0.01d0, &
                 1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0,   1.0d0, &
                0.04d0,  0.02d0,  0.01d0, 0.002d0, 0.002d0, 0.002d0 /), &
             (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   SIG1_GRID=RESHAPE( &
             (/ 0.250d0, 0.120d0, 0.045d0, 0.028d0, 0.020d0, 0.015d0, &
                0.250d0, 0.060d0, 0.035d0, 0.020d0, 0.016d0, 0.012d0, &
                0.170d0, 0.090d0, 0.035d0, 0.020d0, 0.012d0, 0.009d0, &
                0.170d0, 0.070d0, 0.035d0, 0.015d0, 0.012d0, 0.009d0, &
                0.170d0, 0.070d0, 0.050d0, 0.025d0, 0.020d0, 0.020d0, &
                0.065d0, 0.067d0, 0.053d0, 0.032d0, 0.032d0, 0.024d0, &
                0.075d0, 0.044d0,  0.03d0,  0.02d0,  0.02d0, 0.014d0, &
                 10.0d0,  10.0d0,  10.0d0,  10.0d0,  0.02d0,  0.01d0, &
                 10.0d0,  10.0d0,  10.0d0,  10.0d0,  10.0d0,  10.0d0, &
                 0.01d0, 0.005d0, 0.002d0,    1d-4,    1d-4,   1d-4/), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   C1_GRID=  RESHAPE( &
             (/  27.2d0, 27.8d0, 28.2d0, 28.2d0, 28.2d0, 28.2d0, &
                 28.0d0, 27.8d0, 27.8d0, 27.8d0, 27.8d0, 27.8d0, &
                 27.5d0, 27.0d0, 27.8d0, 27.8d0, 27.8d0, 27.8d0, &
                 28.8d0, 28.1d0, 27.8d0, 27.8d0, 27.5d0, 27.5d0, &
                 28.5d0, 28.0d0, 27.5d0, 28.5d0, 29.2d0, 29.0d0, &
                 25.0d0, 27.5d0, 25.8d0, 20.9d0, 29.3d0,  1.0d0, &
                 28.7d0, 27.0d0, 28.0d0, 28.0d0, 27.4d0, 25.3d0, &
                 28.5d0, 29.1d0, 29.5d0, 30.1d0, 30.4d0, 29.9d0, &
                 20.4d0, 20.6d0, 20.8d0, 20.9d0, 20.9d0, 21.0d0, &
                 29.9d0, 30.1d0, 30.1d0, 30.2d0, 30.3d0, 30.3d0 /), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   TAU1_GRID=RESHAPE( &
             (/  4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, &
                 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, &
                 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, &
                 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, 4.07d0, &
                 4.77d0, 4.77d0, 4.77d0, 4.77d0, 4.07d0, 4.07d0, &
                 4.77d0, 4.77d0, 28.2d0, 1.03d0, 0.613d0,1.0d0,  &
                 3.4d0,  14.5d0, 11.4d0, 14.3d0, 13.3d0, 13.3d0, &
                 2.52d0, 2.52d0, 2.52d0, 2.52d0, 2.52d0, 2.52d0, &
                 1.02d0, 1.02d0, 1.02d0, 1.02d0, 1.02d0, 1.02d0, &
                 0.22d0, 0.22d0, 0.22d0, 0.22d0, 0.22d0, 0.22d0 /), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   C2_GRID=  RESHAPE( &
             (/  21.5d0, 21.5d0, 22.1d0, 22.1d0, 22.1d0, 22.1d0, &
                 22.3d0, 21.5d0, 21.5d0, 21.8d0, 21.8d0, 21.8d0, &
                 22.0d0, 21.5d0, 21.5d0, 22.0d0, 21.8d0, 21.8d0, &
                 23.5d0, 22.5d0, 22.1d0, 22.0d0, 22.2d0, 22.2d0, &
                 22.0d0, 22.8d0, 23.0d0, 23.0d0, 23.5d0, 23.5d0, &
                 10.0d0, 0.0d0,  0.0d0,  19.8d0, 22.0d0, 21.0d0, &
                 26.2d0, 14.1d0, 18.8d0, 19.1d0, 23.8d0, 19.2d0, &
                 25.4d0, 25.4d0, 25.8d0, 26.0d0, 26.0d0, 25.8d0, &
                 18.4d0, 18.4d0, 18.6d0, 18.6d0, 18.6d0, 18.6d0, &
                 27.8d0, 28.0d0, 28.2d0, 28.2d0, 28.3d0, 28.3d0 /), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   TAU2_GRID=RESHAPE( &
             (/  4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, &
                 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, &
                 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, &
                 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, 4.62d0, &
                 5.62d0, 5.62d0, 5.62d0, 5.62d0, 4.62d0, 4.62d0, &
                 5.62d0, 5.18d0, 5.18d0, 34.7d0, 8.38d0, 22.6d0, &
                 0.15d0, 4.49d0, 95.0d0, 95.0d0, 0.95d0, 146.d0, &
                 0.12d0, 0.12d0, 0.12d0, 0.12d0, 0.12d0, 0.14d0, &
                 0.32d0, 0.32d0, 0.32d0, 0.32d0, 0.32d0, 0.32d0, &
                 0.02d0, 0.02d0, 0.02d0, 0.02d0, 0.02d0, 0.02d0 /), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   C3_GRID=  RESHAPE( &
             (/  19.4d0, 19.8d0, 20.1d0, 20.1d0, 20.1d0, 20.1d0, &
                 20.0d0, 19.8d0, 19.8d0, 19.8d0, 19.8d0, 19.8d0, &
                 19.9d0, 19.8d0, 19.8d0, 19.8d0, 19.8d0, 19.8d0, &
                  5.9d0,  9.8d0, 23.5d0, 23.5d0, 23.5d0, 23.5d0, &
                 27.3d0, 26.9d0, 26.6d0, 27.4d0, 25.8d0, 25.8d0, &
                 27.8d0, 26.9d0, 18.9d0, 25.4d0, 24.8d0, 25.8d0, &
                 22.8d0, 17.9d0, 18.9d0, 25.4d0, 24.8d0, 25.5d0, &
                 20.6d0, 20.2d0, 19.8d0, 19.2d0, 19.5d0, 18.4d0, &
                 12.6d0, 13.1d0, 14.1d0, 14.5d0, 14.5d0, 14.5d0, &
                 24.3d0, 24.2d0, 24.0d0, 24.0d0, 24.0d0, 23.9d0 /),&
              (/ SIZE(V_GRID), SIZE(YE_GRID) /)), &
   TAU3_GRID=RESHAPE( &
             (/  18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, &
                 18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, &
                 18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, 18.2d0, &
                 18.2d0, 18.2d0, 0.62d0, 0.62d0, 0.62d0, 0.62d0, &
                 0.18d0, 0.18d0, 0.18d0, 0.18d0, 0.32d0, 0.32d0, &
                 0.12d0, 0.18d0, 50.8d0, 0.18d0, 0.32d0, 0.32d0, &
                  2.4d0, 51.8d0, 50.8d0, 0.18d0, 0.32d0, 0.32d0, &
                  3.0d0,  2.5d0,  2.4d0,  2.4d0,  2.4d0, 60.4d0, &
                 200.d0, 200.d0, 200.d0, 200.d0, 200.d0, 200.d0, &
                 8.76d0, 8.76d0, 8.76d0, 8.76d0, 8.76d0, 8.76d0 /), &
              (/ SIZE(V_GRID), SIZE(YE_GRID) /))

CONTAINS

  !************************************************************************
  !                                                                       *
  !  Heating rates restricted to {v, ye}-grid.                            *
  !  OK 29.08.2019                                                        *
  !                                                                       *
  !************************************************************************
  ELEMENTAL DOUBLE PRECISION FUNCTION heating_rate_grid(iv, jye, t) RESULT(h)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: iv, jye
  DOUBLE PRECISION, INTENT(IN):: t
  !
  DOUBLE PRECISION:: e0,alp,t0,sig,alp1,t1,sig1,C1,C2,C3,tau1,tau2,tau3
  DOUBLE PRECISION:: a, b
  DOUBLE PRECISION, PARAMETER:: oneoverpi = .5d0/acos(0d0)

     e0=     E0_GRID(iv,jye)
     alp=   ALP_GRID(iv,jye)
     t0=     T0_GRID(iv,jye)
     sig=   SIG_GRID(iv,jye)
     alp1= ALP1_GRID(iv,jye)
     t1=     T1_GRID(iv,jye)
     sig1= SIG1_GRID(iv,jye)
     C1=     C1_GRID(iv,jye)
     tau1= TAU1_GRID(iv,jye)
     C2=     C2_GRID(iv,jye)
     tau2= TAU2_GRID(iv,jye)
     C3=     C3_GRID(iv,jye)
     tau3= TAU3_GRID(iv,jye)

     a= .5d0 - oneoverpi*atan((t - t0)/sig)
     b= .5d0 + oneoverpi*atan((t - t1)/sig1)
     h= 1d1**(e0+18d0)*(a**alp * b**alp1) &
      + exp(C1 - t/tau1*1d-3) &
      + exp(C2 - t/tau2*1d-5) &
      + exp(C3 - t/tau3*1d-5)

  END FUNCTION heating_rate_grid


  !************************************************************************
  !                                                                       *
  !  Heating rates for arbitrary {v, Ye} at arbitrary time t[s]           *
  !  OK  7.09.2019                                                        *
  !                                                                       *
  !************************************************************************
  DOUBLE PRECISION FUNCTION heating_rate_interpcoef(v, ye, t) RESULT(h)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN):: v  ! ejecta expansion velocity [c]
  DOUBLE PRECISION, INTENT(IN):: ye ! initial electron fraction
  DOUBLE PRECISION, INTENT(IN):: t  ! time [s]
  !
  INTEGER:: i1,i2,j1,j2
  DOUBLE PRECISION:: v1, v2, y1, y2, fv, fy, f11, f12, f21, f22
  DOUBLE PRECISION:: e0,alp,t0,sig,alp1,t1,sig1,C1,C2,C3,tau1,tau2,tau3
  DOUBLE PRECISION:: a, b
  DOUBLE PRECISION, PARAMETER:: oneoverpi = .5d0/acos(0d0)
  DOUBLE PRECISION, PARAMETER:: eps = 1d-13

     IF (v.LT.V_GRID(1) - eps .OR. v.GT.V_GRID(SIZE(V_GRID)) + eps) THEN
        STOP "ERROR: v outside the grid"
     ENDIF
     find_index_v: DO i1= 0, SIZE(V_GRID)-1
        IF (v.LE.V_GRID(i1+1)) EXIT find_index_v
     ENDDO find_index_v
     i2= i1 + 1

     IF (ye.LT.YE_GRID(1) - eps .OR. ye.GT.YE_GRID(SIZE(YE_GRID)) + eps) THEN
        STOP "ERROR: ye outside the grid"
     ENDIF
     find_index_ye: DO j1= 0, SIZE(YE_GRID)-1
        IF (ye.LE.YE_GRID(j1+1)) EXIT find_index_ye
     ENDDO find_index_ye
     j2= j1 + 1

     v1= V_GRID(i1)
     v2= V_GRID(i2)
     fv= (v - v1)/(v2 - v1)

     y1= YE_GRID(j1)
     y2= YE_GRID(j2)
     fy= (ye - y1)/(y2 - y1)

     f11= (1d0 - fv)*(1d0 - fy)
     f12= (1d0 - fv)*fy
     f21= fv*(1d0 - fy)
     f22= fv*fy

     e0=   f11*E0_GRID(i1,j1) + f12*E0_GRID(i1,j2) &
         + f21*E0_GRID(i2,j1) + f22*E0_GRID(i2,j2)

     alp=  f11*ALP_GRID(i1,j1) + f12*ALP_GRID(i1,j2) &
         + f21*ALP_GRID(i2,j1) + f22*ALP_GRID(i2,j2)

     t0=   f11*T0_GRID(i1,j1) + f12*T0_GRID(i1,j2) &
         + f21*T0_GRID(i2,j1) + f22*T0_GRID(i2,j2)

     sig=  f11*SIG_GRID(i1,j1) + f12*SIG_GRID(i1,j2) &
         + f21*SIG_GRID(i2,j1) + f22*SIG_GRID(i2,j2)

     alp1= f11*ALP1_GRID(i1,j1) + f12*ALP1_GRID(i1,j2) &
         + f21*ALP1_GRID(i2,j1) + f22*ALP1_GRID(i2,j2)

     t1=   f11*T1_GRID(i1,j1) + f12*T1_GRID(i1,j2) &
         + f21*T1_GRID(i2,j1) + f22*T1_GRID(i2,j2)

     sig1= f11*SIG1_GRID(i1,j1) + f12*SIG1_GRID(i1,j2) &
         + f21*SIG1_GRID(i2,j1) + f22*SIG1_GRID(i2,j2)

     C1=   f11*C1_GRID(i1,j1) + f12*C1_GRID(i1,j2) &
         + f21*C1_GRID(i2,j1) + f22*C1_GRID(i2,j2)

     tau1= f11*TAU1_GRID(i1,j1) + f12*TAU1_GRID(i1,j2) &
         + f21*TAU1_GRID(i2,j1) + f22*TAU1_GRID(i2,j2)

     C2=   f11*C2_GRID(i1,j1) + f12*C2_GRID(i1,j2) &
         + f21*C2_GRID(i2,j1) + f22*C2_GRID(i2,j2)

     tau2= f11*TAU2_GRID(i1,j1) + f12*TAU2_GRID(i1,j2) &
         + f21*TAU2_GRID(i2,j1) + f22*TAU2_GRID(i2,j2)

     C3=   f11*C3_GRID(i1,j1) + f12*C3_GRID(i1,j2) &
         + f21*C3_GRID(i2,j1) + f22*C3_GRID(i2,j2)

     tau3= f11*TAU3_GRID(i1,j1) + f12*TAU3_GRID(i1,j2) &
         + f21*TAU3_GRID(i2,j1) + f22*TAU3_GRID(i2,j2)

     a= .5d0 - oneoverpi*atan((t - t0)/sig)
     b= .5d0 + oneoverpi*atan((t - t1)/sig1)
     h= 1d1**(e0+18d0)*(a**alp * b**alp1) &
      + exp(C1 - t/tau1*1d-3) &
      + exp(C2 - t/tau2*1d-5) &
      + exp(C3 - t/tau3*1d-5)

  END FUNCTION heating_rate_interpcoef


  !************************************************************************
  !                                                                       *
  !  Heating rates for arbitrary {v, Ye} at arbitrary time t[s],          *
  !  by interpolating the function values instead of the variables        *
  !  AK  25.04.2024                                                       *
  !                                                                       *
  !************************************************************************
  ELEMENTAL DOUBLE PRECISION FUNCTION heating_rate(v, ye, t) RESULT(h)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN):: v  ! ejecta expansion velocity [c]
  DOUBLE PRECISION, INTENT(IN):: ye ! initial electron fraction
  DOUBLE PRECISION, INTENT(IN):: t  ! time [s]
  !
  INTEGER:: i1,i2,j1,j2
  DOUBLE PRECISION:: v1, v2, y1, y2, fv, fy, f11, f12, f21, f22
  DOUBLE PRECISION:: e0,alp,t0,sig,alp1,t1,sig1,C1,C2,C3,tau1,tau2,tau3
  DOUBLE PRECISION:: a, b
  DOUBLE PRECISION, PARAMETER:: oneoverpi = .5d0/acos(0d0)
  DOUBLE PRECISION, PARAMETER:: eps = 1d-13

     find_index_v: DO i1= 0, SIZE(V_GRID)-1
        IF (v.LE.V_GRID(i1+1)) EXIT find_index_v
     ENDDO find_index_v
     i2= i1 + 1

     find_index_ye: DO j1= 0, SIZE(YE_GRID)-1
        IF (ye.LE.YE_GRID(j1+1)) EXIT find_index_ye
     ENDDO find_index_ye
     j2= j1 + 1

     v1= V_GRID(i1)
     v2= V_GRID(i2)
     fv= (v - v1)/(v2 - v1)

     y1= YE_GRID(j1)
     y2= YE_GRID(j2)
     fy= (ye - y1)/(y2 - y1)

     f11= (1d0 - fv)*(1d0 - fy)
     f12= (1d0 - fv)*fy
     f21= fv*(1d0 - fy)
     f22= fv*fy

     h= dexp(f11*dlog(heating_rate_grid(i1,j1,t)) &
           + f12*dlog(heating_rate_grid(i1,j2,t)) &
           + f21*dlog(heating_rate_grid(i2,j1,t)) &
           + f22*dlog(heating_rate_grid(i2,j2,t)))

  END FUNCTION heating_rate


  !************************************************************************
  !                                                                       *
  !  Output heating rates for arbitrary {v, Ye} in SuperNu format         *
  !  OK  12.12.2019                                                       *
  !                                                                       *
  !************************************************************************
  SUBROUTINE output_supernu(v, ye)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN):: v  ! ejecta expansion velocity [c]
  DOUBLE PRECISION, INTENT(IN):: ye ! initial electron fraction
  !
  INTEGER, PARAMETER :: nmax = 1000
  DOUBLE PRECISION, PARAMETER:: &
    a=  0.05d0,      & !  5% alphas
    b=  0.20d0,      & ! 20% betas
    nu= 0.35d0,      & ! 35% neutrinos
    rad= a + b + nu, & ! gammas
    gam= 1d0 - rad,  & ! gammas
    tmin= 1d-3,      & ! [s] initial time
    tmax= 86400*100    ! [s] final time (100 days)
  INTEGER :: n
  DOUBLE PRECISION:: t, dtfac, hr

  PRINT ('("# Heating rates from Stephan Rosswog (WinNet),",I2,"% alpha'// &
         's, ",I2,"% betas, ",I2,"% nu")'),INT(a*100),INT(b*100),INT(nu*100)
  PRINT ('("# 1st block is thermalization parameters, 2nd block is heat'// &
         'ing rates")')
  PRINT *
  PRINT ('("# 1st line is thermalization option for columns 4-8 in the '// &
         'heating rates")')
  PRINT ('("# 2nd line is thermalization parameter (use depends on opti'// &
         'on in 1st line)")')
  PRINT ('("           1             1             0             1     '// &
         '        1")')
  PRINT ('("1.200000E-11  1.300000E-11  1.000000E+00  1.300000E-11  2.0'// &
         '00000E-12")')
  PRINT *
  PRINT ('("# ",I4)'), nmax
  PRINT ('("# 1:time[s], 2:total nuclear heating rate[erg/(g*s)] 3:tota'// &
         'l radiation[erg/g*s]")')
  PRINT ('("# 4:alpha 5:beta 6:gamma 7:electrons 8:fission products")')

  t= tmin
  dtfac= exp((log(tmax) - log(tmin))/(dble(nmax) - 1))
  DO n=1,1000
     hr= heating_rate(v, ye, t)
     PRINT '(12(ES12.5,1X))', t, hr, hr*rad, hr*a, hr*b, hr*gam, 0d0, 0d0
     t= t*dtfac
  ENDDO

  END SUBROUTINE output_supernu

endmodule rprocmod

