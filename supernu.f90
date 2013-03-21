program supernu

  use mpimod
  use inputparmod
  use timestepmod
  use gasgridmod
  use particlemod
  use physconstmod

  use ionsmod, only:ion_read_data,ion_alloc_grndlev
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod

  implicit none
!***********************************************************************
! Main routine
!
! todo:
!  - dummy (drr, 2013/mar/05)
!***********************************************************************
  real*8 :: help, dt
  real*8 :: t_elapsed
  integer :: ierr
  logical :: lmpi0 = .false. !master rank flag
  real :: t0,t1  !timing
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
!-- parse and verify runtime parameters
   call parse_inputpars(nmpi)
!
!-- init random number generator
   help = rand(in_seed)
!
!-- time step init
!-- constant time step, may be coded to loop if time step is not uniform
   t_elapsed = (in_tlast - in_tfirst) * pc_day  !convert input from days to seconds
   dt = t_elapsed/in_nt
   call timestep_init(in_nt,in_alpha,in_tfirst,dt)
!
!-- particle init
   call particle_init(in_npartmax,in_ns)
!
!-- SETUP GRIDS
   call gasgrid_init(in_nt)
   call gasgrid_setup
   call wlgrid_setup
!-- read initial temperature structure from file
!   call read_restart_file
!-- hard coded temperature structure
!   do ir = 1, gas_nr
!     if (gas_vals2(ir)%tempkev<1.e-6) then
!       gas_vals2(ir)%tempkev = 1.e-6
!     endif
!   enddo

!-- READ DATA
!-- read ion and level data
   call ion_read_data(gas_nelem)  !ion and level data
   call ion_alloc_grndlev(gas_nr)  !ground state occupation numbers
!-- read bbxs data
   if(.not.in_nobbopac) call read_bbxs_data(gas_nelem)!bound-bound cross section data
!-- read bfxs data
   if(.not.in_nobfopac) call bfxs_read_data           !bound-free cross section data
!-- read bfxs data
   if(.not.in_noffopac) call ffxs_read_data           !free-free cross section data
!
   call time(t1)
   t_setup = t1-t0
  endif !impi

  gas_ng = 2
  ! Beginning time step loop
  do tsp_tn = 1, tsp_nt
    write(6,'(a,i5,f8.3,"d")') 'timestep:',tsp_tn,tsp_texp/pc_day
    !
    call gasgrid_update
    !Calculating Fleck factor, leakage opacities
    call xsections
    !Calculating number of source prt_particles per cell
    call sourcenumbers
    !Storing vacant "prt_particles" indexes in ordered array "prt_vacantarr"
    allocate(prt_vacantarr(prt_nnew))
    call vacancies
    !Calculating properties of prt_particles on domain boundary
    call boundary_source
    !Calculating properties of prt_particles emitted in domain interior
    call interior_source
    deallocate(prt_vacantarr)
    !Advancing prt_particles to update radiation field
    call advance

    call temperature_update
    call timestep_update(dt) !Updating elapsed tsp_time and expansion tsp_time

    call write_output
!   call write_restart
  enddo
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

end program supernu
