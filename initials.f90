SUBROUTINE initials

  USE gasgridmod
  USE particlemod
  USE timestepmod
  USE physconstmod
  USE inputparmod

  IMPLICIT NONE

  INTEGER(iknd) :: ir, ipart, ig
  REAL(rknd) :: Um

  prt_erad = 0.0   !Total radiation energy
  prt_einit = 0.0  !Total initial energy
  prt_einp = 0.0   !Input Energy
  prt_eint = 0.0   !Total internal energy
  !OPEN(UNIT=8,FILE='Tinit.dat',STATUS='UNKNOWN')
  !READ(8,*) gas_temp
  !DO ir = 1, gas_nr
  !   IF (gas_temp(ir)<1.e-6) THEN
  !      gas_temp(ir) = 1.e-6
  !   ENDIF
  !ENDDO
  nidecay = 1.73*(1.6022e-6)  !erg/s/g

  DO ir = 1, gas_nr
     gas_rhoarr(ir) = 2.4186e8 !g/cm^3
     gas_temp(ir) = 1.e3 !861.73
     !gas_bcoef(ir) = 2.0*a_coef*gas_temp(ir)**3

     gas_bcoef(ir) = 0.4*(1.e12*gas_rhoarr(ir))*580.25_rknd

     gas_ur(ir) = a_coef*gas_temp(ir)**4
     Um = gas_bcoef(ir)*gas_temp(ir)
     prt_einit = prt_einit + Um*4*pi*gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**3)/3.0
  ENDDO
  prt_einp = prt_einit
  DO ipart = 1, prt_npartmax
     prt_particles(ipart)%isvacant=.TRUE.
  ENDDO

  tsp_dt = t_elapsed/tsp_nt
  tsp_time = RAND(seed)   !PRNG initial
  tsp_time = 0.0

END SUBROUTINE initials 
