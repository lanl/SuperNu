SUBROUTINE initials

  USE gasgridmod
  USE particlemod
  USE timestepmod
  USE physconstmod
  USE inputparmod

  IMPLICIT NONE
!##################################################
  !This routine initializes material (gas) properties and time step size.
!##################################################

  INTEGER :: ir, ipart
  REAL*8 :: Um, t_elapsed

  gas_erad = 0.0   !Total radiation energy
  gas_einit = 0.0  !Total initial energy
  gas_einp = 0.0   !Input Energy
  gas_eint = 0.0   !Total internal energy
  !OPEN(UNIT=8,FILE='Tinit.dat',STATUS='UNKNOWN')
  !READ(8,*) gas_temp
  !DO ir = 1, gas_nr
  !   IF (gas_vals2(ir)%tempkev<1.e-6) THEN
  !      gas_vals2(ir)%tempkev = 1.e-6
  !   ENDIF
  !ENDDO

  ! Initial gas temperature, density, and heat capacity generation loop
  DO ir = 1, gas_nr
     gas_vals2(ir)%rho = 2.4186e8 !g/cm^3
     gas_vals2(ir)%tempkev = 1.e3 !861.73
     !gas_vals2(ir)%bcoef = 2.0*pc_acoef*gas_vals2(ir)%tempkev**3

     gas_vals2(ir)%bcoef = 0.4*(1.e12*gas_vals2(ir)%rho)*580.25d0

     gas_vals2(ir)%ur = pc_acoef*gas_vals2(ir)%tempkev**4
     Um = gas_vals2(ir)%bcoef*gas_vals2(ir)%tempkev
     gas_einit = gas_einit + Um*4*pc_pi*gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**3)/3.0
  ENDDO
  gas_einp = gas_einit
  ! Setting all entries of particle array to vacant: loop
  DO ipart = 1, prt_npartmax
     prt_particles(ipart)%isvacant=.TRUE.
  ENDDO

  ! Calculating elapsed physical time and time step size of problem (may be coded to loop if
  ! time step is not uniform)
  t_elapsed = (in_tlast - in_tfirst) * pc_day  !convert input from days to seconds
  tsp_dt = t_elapsed/tsp_nt
  tsp_time = RAND(in_seed)   !PRNG initial
  tsp_time = 0.0

END SUBROUTINE initials 
