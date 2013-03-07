PROGRAM supernu

  use mpimod
  USE inputparmod
  USE timestepmod
  USE gasgridmod
  USE particlemod

  use ionsmod, only:ion_read_data,ion_alloc_grndlev
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod

  IMPLICIT NONE
!***********************************************************************
! Main routine
!
! todo:
!  - dummy (drr, 2013/mar/05)
!***********************************************************************
  INTEGER :: it
  REAL*8 :: time_begin, time_end
  integer :: ierr
  LOGICAL :: lmpi0 = .false. !master rank flag
  REAL :: t0,t1  !timing
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
!
!--
!-- SETUP SIMULATION:
!====================
!-- The setup is done by the master task only, and broadcasted to the
!-- other tasks before packet propagation begins.
!--
  if(impi==impi0) then
   lmpi0 = .true. !master rank flag
   call time(t0)
!-- startup message
   call banner
!-- read runtime parameters
   call read_inputpars
!-- parse runtime parameters
   call parse_inputpars(nmpi)
!
!-- time step init
   CALL timestep_init(in_nt,in_alpha)
!
!-- particle init
   CALL particle_init(in_npartmax,in_ns)
!
!-- SETUP GRIDS
!-- setup gas grid (map gstruct to gasgrid)
   CALL gasgrid_init(in_nr,in_ng,in_lr,in_isvelocity)
!! call ggrid_setup(gas_nr,tim_ntim)
!-- read initial temperature structure from file
   !call read_restart_file
!
!-- READ DATA
!-- read ion and level data
   call ion_read_data(gas_nelem)  !ion and level data
   call ion_alloc_grndlev(gas_nr)   !ground state occupation numbers
!-- read bbxs data
   if(.not.in_nobbopac) call read_bbxs_data(gas_nelem)!bound-bound cross section data
!-- read bfxs data
   if(.not.in_nobfopac) call bfxs_read_data          !bound-free cross section data
!-- read bfxs data
   if(.not.in_noffopac) call ffxs_read_data          !free-free cross section data
!
   call time(t1)
   t_setup = t1-t0
  endif !impi
  
  
  ! Setting velocity option
  IF (in_isvelocity.EQV..TRUE.) THEN
     gas_velyes = 1
     gas_velno = 0
  ELSE
     gas_velyes = 0
     gas_velno = 1
  ENDIF
  ! Setting transport option
  in_puretran = .FALSE.

  CALL grids
  CALL initials

  ! Beginning time step loop
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
!    CALL write_restart
     tsp_tn = tsp_tn+1
  ENDDO
  CALL CPU_TIME(time_end)
  WRITE(*,*) 'CPU TIME: ',time_end-time_begin,' seconds'
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output.
  if(lmpi0) then
   call time(t1)
   t_all = dble(t1 - t0)
   call print_timing                 !print timing results
   write(6,*)
   write(6,*) 'SuperNu finished'
   if(in_grab_stdout)write(0,'(a,f8.2,"s")')'SuperNu finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This help to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

END PROGRAM supernu
