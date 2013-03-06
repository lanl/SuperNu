PROGRAM driver

  USE inputparmod
  USE timestepmod
  USE gasgridmod
  USE particlemod
  IMPLICIT NONE

  INTEGER(iknd) :: it, ipart
  REAL(rknd) :: time_begin, time_end
  
  ! Namelist groups
  NAMELIST / scalarins / t_elapsed, Lr
  NAMELIST / arrayins / in_nt, in_nr, in_ng, isvelocity
  NAMELIST / mcins / in_Ns, in_Npartmax, alpha, seed
  ! Reading namelist (from file "input")
  OPEN(UNIT=1,FILE="input",STATUS="OLD",FORM="FORMATTED")
  READ(UNIT=1,NML=scalarins)
  READ(UNIT=1,NML=arrayins)
  READ(UNIT=1,NML=mcins)

  CALL gasgrid_init(in_nr,in_ng)
  CALL timestep_init(in_nt)
  CALL particle_init(in_Npartmax,in_Ns)
  
  ! Setting velocity option
  IF (isvelocity.EQV..TRUE.) THEN
     gas_velyes = 1
     gas_velno = 0
  ELSE
     gas_velyes = 0
     gas_velno = 1
  ENDIF
  ! Setting transport option
  puretran = .FALSE.

  CALL globalallocations
  CALL grids
  CALL initials

  CALL CPU_TIME(time_begin)
  tsp_tn = 1
  DO it = 1, tsp_nt 
     WRITE(*,*) 'timestep:',it
     !Calculating opacities (for IMC(transport) and DDMC(diffusion))
     CALL xsections
     !Calculating number of source prt_particles per cell
     CALL sourcenumbers
     !Storing vacant "prt_particles" indexes in ordered array "prt_vacantarr"
     ALLOCATE(prt_vacantarr(prt_nnew))
     CALL vacancies
     !Calculating properties of prt_particles on domain boundary
     !CALL boundary_source
     !Calculating properties of prt_particles emitted in domain interior
     CALL interior_source
     DEALLOCATE(prt_vacantarr)
     !Advancing prt_particles to update radiation field
     CALL advance
     !Updating material state
     CALL material_update
     !Updating elapsed tsp_time and expansion tsp_time
     tsp_time = tsp_time+tsp_dt
     tsp_texp = tsp_texp+tsp_dt
     !Writing data to files
     CALL write_output
     tsp_tn = tsp_tn+1
  ENDDO
  CALL CPU_TIME(time_end)
  WRITE(*,*) 'CPU TIME: ',time_end-time_begin,' seconds'

END PROGRAM driver
