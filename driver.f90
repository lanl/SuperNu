PROGRAM driver

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

  INTEGER(iknd) :: it, ipart
  REAL(rknd) :: time_begin, time_end
  logical :: lmpi0 = .false. !master rank flag
  real :: t0,t1  !timing
************************************************************************
* Main routine
*
* todo:
*  - dummy (drr, 2013/mar/05)
************************************************************************
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
  CALL timestep_init(in_nt)
!
!-- particle init
  CALL particle_init(in_npartmax,in_ns)
!
!-- SETUP GRIDS
!-- setup gas grid (map gstruct to ggrid)
   CALL gasgrid_init(in_nr,in_ng)
   call ggrid_setup(ncg,tim_ntim)
!-- read initial temperature structure from file
   !call read_restart_file
!
!-- READ DATA
!-- read ion and level data
   call ion_read_data(gas_nelem)  !ion and level data
   call ion_alloc_grndlev(ncg)   !ground state occupation numbers
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
     velyes = 1
     velno = 0
  ELSE
     velyes = 0
     velno = 1
  ENDIF
  ! Setting transport option
  in_puretran = .FALSE.

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

END PROGRAM driver
