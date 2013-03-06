SUBROUTINE initials

  USE gasgridmod
  USE particlemod
  USE timestepmod
  USE physconstmod
  USE inputparmod

  IMPLICIT NONE

  INTEGER(iknd) :: ir, ipart, ig
  REAL(rknd) :: Um, t_elapsed

  Erad = 0.0   !Total radiation energy
  Einit = 0.0  !Total initial energy
  Einp = 0.0   !Input Energy
  Eint = 0.0   !Total internal energy
  !OPEN(UNIT=8,FILE='Tinit.dat',STATUS='UNKNOWN')
  !READ(8,*) Temp
  !DO ir = 1, gas_nr
  !   IF (Temp(ir)<1.e-6) THEN
  !      Temp(ir) = 1.e-6
  !   ENDIF
  !ENDDO

  DO ir = 1, gas_nr
     rhoarr(ir) = 2.4186e8 !g/cm^3
     Temp(ir) = 1.e3 !861.73
     !bcoef(ir) = 2.0*a_coef*Temp(ir)**3

     bcoef(ir) = 0.4*(1.e12*rhoarr(ir))*580.25_rknd

     Ur(ir) = a_coef*Temp(ir)**4
     Um = bcoef(ir)*Temp(ir)
     Einit = Einit + Um*4*pi*dr3arr(ir)*(velno*1.0+velyes*texp**3)/3.0
  ENDDO
  Einp = Einit
  DO ipart = 1, prt_Npartmax
     particles(ipart)%isvacant=.TRUE.
  ENDDO

  t_elapsed = in_tlast - in_tfirst
  dt = t_elapsed/tsp_nt
  time = RAND(in_seed)   !PRNG initial
  time = 0.0

END SUBROUTINE initials 
