program supernu

  use mpimod
  use inputparmod
  use timestepmod
  use gasgridmod
  use particlemod
  use physconstmod
  use profiledatamod
  use miscmod

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
! - check wavelength units for bb and bf data
! - fix gas_wl indexing BUG in physical_opacity
!
!***********************************************************************
  real*8 :: help,dt
  real*8 :: t_elapsed
  integer :: ierr,ng,ns,it
  integer,external :: memusg
  logical :: lmpi0 = .false. !master rank flag
  real*8 :: t0,t1  !timing
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
!
!-- initialize timing module
  call timing_init
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
   call banner
!-- read runtime parameters
   call read_inputpars
!-- parse and verify runtime parameters
   call parse_inputpars(nmpi)
!
!-- time step init
   call timestep_init(in_nt,in_ntres,in_alpha,in_tfirst)
!-- constant time step, may be coded to loop if time step is not uniform
   t_elapsed = (in_tlast - in_tfirst) * pc_day  !convert input from days to seconds
   dt = t_elapsed/in_nt
!
!-- particle init
   ns = in_ns/nmpi
   call particle_init(in_npartmax,ns,in_ns0,in_isimcanlog, &
        in_isddmcanlog,in_tauddmc,in_taulump,in_tauvtime,nmpi, &
        in_norestart)
!
!-- rand() count and prt restarts
   if(tsp_ntres>1.and..not.in_norestart) then
!-- read rand() count
     call read_restart_randcount
!-- read particle properties
     call read_restart_particles
   endif
!
!-- read input structure
   if(.not.in_noreadstruct.and.in_isvelocity) then
     call read_inputstr(in_igeom,in_ndim)
   else
!== generate_inputstr development in progress
     call generate_inputstr(in_igeom)
   endif
!
!-- read gamma deposition profiles
   if(in_isvelocity.and.in_srctype=='none') then
     if(in_igeom>1) stop 'supernu: read_gam_prof: no 2D/3D'
     call read_gamma_profiles(in_ndim)
   endif
!
!-- SETUP GRIDS
   call wlgrid_setup(ng)
   call gasgrid_init(ng)
   call gasgrid_setup
!

!-- READ DATA
!-- read ion and level data
   call ion_read_data(gas_nelem)  !ion and level data
   call ion_alloc_grndlev(gas_nx,gas_ny,gas_nz)  !ground state occupation numbers
!-- read bbxs data
   if(.not.in_nobbopac) call read_bbxs_data(gas_nelem)!bound-bound cross section data
!-- read bfxs data
   if(.not.in_nobfopac) call bfxs_read_data           !bound-free cross section data
!-- read ffxs data
   if(.not.in_noffopac) call ffxs_read_data           !free-free cross section data
!
!-- initial radiation energy
   call initialnumbers
!
!-- memory statistics
   write(6,*) 'memusg: after setup:',memusg()
!
   call time(t1)
   t_setup = t1-t0!}}}
  endif !impi

  call bcast_permanent !MPI
!
!-- initialize random number generator, use different seeds for each rank
  if(in_nomp==0) stop 'supernu: in_nomp == 0'
  help = rand(in_nomp*impi)
!
!-- reading restart rand() count
  if(tsp_ntres>1.and..not.in_norestart) then
    call scatter_restart_data !MPI
!-- mimicking end of tsp reset
    prt_tlyrand = 0
  else
    prt_tlyrand = 1
  endif
!
!-- instantiating initial particles (if any)
  call initial_particles
!
!
!-- time step loop
!=================
  if(impi==impi0) then
     write(6,*)
     write(6,*) "starting time loop:"
     write(6,*) "===================="
  endif
!
  do it=in_ntres,tsp_nt
!-- allow negative and zero it for temperature initialization purposes!{{{
    tsp_it = max(it,1)
!
!-- reset particle clocks
    if(tsp_it<=tsp_ntres) then
      where(.not.prt_isvacant) prt_particles%tsrc = tsp_t
    endif
!
!-- Update tsp_t etc
    if(impi==impi0) call timestep_update(dt)
!
!-- updating prt_tauddmc and prt_taulump
    call tau_update

    if(impi==impi0) then
      write(6,'(1x,a,i5,f8.3,"d",i12)') 'timestep:',it,tsp_t/pc_day, &
         count(.not.prt_isvacant)
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
!-- Calculating properties of prt_particles on domain boundary
    call boundary_source
!-- Calculating properties of prt_particles emitted in domain interior
    call interior_source
    
    deallocate(prt_vacantarr)

!-- advance particles
    call particle_advance

!-- collect particle results from all workers
    call reduce_tally !MPI
!
!-- print packet advance load-balancing info
    !if(impi==impi0) write(6,'(1x,a,3(f9.2,"s"))') 'packets time(min|mean|max):',t_pckt_stat
    if(impi==impi0) call timereg(t_pckt_allrank,t_pckt_stat(3))  !register particle advance timing

!-- collect data necessary for restart (tobe written by impi0)
    if(.not.in_norestart) call collect_restart_data !MPI
!
    if(impi==impi0) then
!-- luminosity statistics
      where(gas_lumnum>0) gas_lumdev = ( &
         (gas_lumdev/gas_lumnum - (gas_luminos/gas_lumnum)**2) * &
         gas_lumnum/dble(gas_lumnum - 1) )**.5d0
!-- update temperature
      call temperature_update
!-- check energy (particle weight) is accounted
      call energy_check
!-- write output
      if(it>0) call write_output
!
!-- restart writers
      if(.not.in_norestart .and. it>0) then
         call write_restart_file !-- temp
         call write_restart_randcount !-- rand() count
         call write_restart_particles !-- particle properties of current time step
      endif
!
    endif !impi
!
!-- reset rand counters
    prt_tlyrand = 0
!
!-- write timestep timing to file
    if(it>0) call timing_timestep(impi)
!!}}}
  enddo !tsp_it
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output.
  if(lmpi0) then
   call time(t1)
   t_all = t1 - t0
   call print_timing                 !print timing results
   write(6,*)
   write(6,*) 'SuperNu finished'
   if(in_grab_stdout)write(0,'(a,f8.2,"s")')'SuperNu finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This help to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

end program supernu
