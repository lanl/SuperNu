program supernu

  use mpimod
  use inputparmod
  use timestepmod
  use gasgridmod
  use particlemod
  use physconstmod

  use inputstrmod

  use ionsmod, only:ion_read_data,ion_alloc_grndlev
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod

  implicit none
!***********************************************************************
! Main routine
!
! todo:
! - clean up manual stratification usage throughout code (into dedicated
!   subroutine?) (13/05/31)
! - verify linear and uniform velocity grid in input.str (13/05/31)
! - check wavelength units for bb and bf data
! - fix gas_wl indexing BUG in physical_opacity
! - improve grid construction: store unit-sphere radii separately (13/06/01)
! - clean up dealloc_all and add all arrays allocated in mpimod_mpi.f
!
!***********************************************************************
  real*8 :: help, dt
  real*8 :: t_elapsed
  integer :: ierr, iimpi, iseed
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
   lmpi0 = .true. !master rank flag!{{{
   call time(t0)
!-- startup message
   write(0,*) 'here'
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
   call particle_init(in_npartmax,in_ns,in_isimcanlog,in_isddmcanlog,in_tauddmc)
!
!-- read input structure
   if(.not.in_noreadstruct) then
    call read_inputstr(in_nr,gas_velout)
   else
    gas_velout = in_velout
   endif
!
!-- SETUP GRIDS
   call wlgrid_setup
   call gasgrid_init(in_nt)
   call gasgrid_setup
!

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
   t_setup = t1-t0!}}}
  endif !impi

  call bcast_permanent !MPI

!-- calculating analytic initial particle distribution (if any)
! not working currently, as prt_particle data structure is only
! available after bcast_nonpermanent
!  call initialnumbers
!
!-- time step loop
!=================
  do iimpi = 1, impi
     iseed = 10*iimpi+3
     help = rand(iseed)
  enddo
  do tsp_it = 1, tsp_nt
    if(impi==impi0) then
      write(6,'(a,i5,f8.3,"d")') 'timestep:',tsp_it,tsp_texp/pc_day
!-- update all non-permanent variables
      call gasgrid_update
!-- number of source prt_particles per cell
      call sourcenumbers
    endif !impi

!-- broadcast to all workers
    call bcast_nonpermanent !MPI

!-- Storing vacant "prt_particles" indexes in ordered array "prt_vacantarr"
    allocate(prt_vacantarr(prt_nnew))
    call vacancies
    !Calculating properties of prt_particles on domain boundary
    call boundary_source
    !Calculating properties of prt_particles emitted in domain interior
write(0,*) 'test1',impi
    call interior_source
write(0,*) 'test2',impi
    deallocate(prt_vacantarr)
    !Advancing prt_particles to update radiation field    

!-- advance particles
    call particle_advance

!-- collect particle results from all workers
    call reduce_tally !MPI
    
    if(impi==impi0) then
       ! averaging reduced results
       !if(nmpi>1) then
!-- dim==0
          gas_erad = gas_erad/real(nmpi)
          gas_eright = gas_eright/real(nmpi)
          gas_eleft = gas_eleft/real(nmpi)
!-- dim==1
          gas_edep = gas_edep/real(nmpi)
!-- dim==2
          gas_eraddens = gas_eraddens/real(nmpi)
       !endif
       !
      call temperature_update
      call timestep_update(dt) !Updating elapsed tsp_time and expansion tsp_time

      call write_output
      call write_restart_file
    endif !impi
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
  write(6,*) 'here 1'
!-- Clean up memory. (This help to locate memory leaks)
  call dealloc_all
  write(6,*) 'here 2'
  call mpi_finalize(ierr) !MPI
  write(6,*) 'here 3'

end program supernu
