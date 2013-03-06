PROGRAM driver

  USE inputparmod
  USE timestepmod
  USE gasgridmod
  USE particlemod
  IMPLICIT NONE

  INTEGER(iknd) :: it, ipart
  REAL(rknd) :: time_begin, time_end
  
  CALL gasgrid_init(in_nr,in_ng)
  CALL particle_init(in_npartmax,in_ns)
  CALL timestep_init(in_nt)
  
  ! Setting velocity option
  IF (in_isvelocity.EQV..TRUE.) THEN
     velyes = 1
     velno = 0
  ELSE
     velyes = 0
     velno = 1
  ENDIF
  ! Setting transport option
  puretran = .FALSE.

  CALL globalallocations
  CALL grids
  CALL initials

  CALL CPU_TIME(time_begin)
  tn = 1
  DO it = 1, tsp_nt 
     WRITE(*,*) 'timestep:',it
     !Calculating opacities (for IMC(transport) and DDMC(diffusion))
     CALL xsections
     !Calculating number of source particles per cell
     CALL sourcenumbers
     !Storing vacant "particles" indexes in ordered array "vacantarr"
     ALLOCATE(vacantarr(Nnew))
     CALL vacancies
     !Calculating properties of particles on domain boundary
     !CALL boundary_source
     !Calculating properties of particles emitted in domain interior
     CALL interior_source
     DEALLOCATE(vacantarr)
     !Advancing particles to update radiation field
     CALL advance
     !Updating material state
     CALL material_update
     !Updating elapsed time and expansion time
     time = time+dt
     texp = texp+dt
     !Writing data to files
     CALL write_output
     tn = tn+1
  ENDDO
  CALL CPU_TIME(time_end)
  WRITE(*,*) 'CPU TIME: ',time_end-time_begin,' seconds'

END PROGRAM driver
