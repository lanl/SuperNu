SUBROUTINE initials

  USE gasgridmod
  USE particlemod
  USE timestepmod
  USE physconstmod
  USE inputparmod

  IMPLICIT NONE

  INTEGER(iknd) :: ir, ipart, ig
  REAL(rknd) :: Um, t_elapsed

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

  DO ir = 1, gas_nr
     gas_rhoarr(ir) = 2.4186e8 !g/cm^3
     gas_temp(ir) = 1.e3 !861.73
     !gas_bcoef(ir) = 2.0*pc_acoef*gas_temp(ir)**3

     gas_bcoef(ir) = 0.4*(1.e12*gas_rhoarr(ir))*580.25_rknd

     gas_ur(ir) = pc_acoef*gas_temp(ir)**4
     Um = gas_bcoef(ir)*gas_temp(ir)
     prt_einit = prt_einit + Um*4*pc_pi*gas_dr3arr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**3)/3.0
  ENDDO
  prt_einp = prt_einit
  DO ipart = 1, prt_npartmax
     prt_particles(ipart)%isvacant=.TRUE.
  ENDDO

  t_elapsed = (in_tlast - in_tfirst) * pc_day  !convert input from days to seconds
  tsp_dt = t_elapsed/tsp_nt
  tsp_time = RAND(in_seed)   !PRNG initial
  tsp_time = 0.0

END SUBROUTINE initials 
