program supernu

  use randommod
  use sourcemod
  use mpimod
  use transportmod
  use inputparmod
  use timestepmod
  use groupmod
  use gridmod
  use gasmod
  use particlemod
  use physconstmod
  use miscmod
  use totalsmod

  use inputstrmod
  use fluxmod

  use ionsmod, only:ions_read_data,ions_alloc_grndlev
  use bfxsmod, only:bfxs_read_data
  use ffxsmod, only:ffxs_read_data
  use timingmod

  implicit none
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  real*8 :: help
  real*8 :: dt
  integer :: ierr, it
  integer :: icell1, ncell !number of cells per rank (gas_ncell)
  integer :: ns, nsinit, mpart
  real*8 :: tfirst, tlast
  real*8 :: t0, t1 !timing
  character(15) :: msg
!
!-- mpi initialization
  call mpi_init(ierr) !MPI
  call mpi_comm_rank(MPI_COMM_WORLD,impi,ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr) !MPI
  lmpi0 = impi==impi0
!
!-- initialize timing module
  call timingmod_init
!
!--
!-- read and distribut input data
!================================
!-- this is done by the master task only and then broadcasted
  if(lmpi0) then
     t0 = t_time()!{{{
!-- startup message
     call banner
!-- read runtime parameters
     call read_inputpars
!-- parse and verify runtime parameters
     call parse_inputpars(nmpi)
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
       call read_inputstr(in_igeom,in_ndim,in_voidcorners,nmpi)
     else
!== generate_inputstr development in progress
       call generate_inputstr(in_igeom)
     endif
!-- compressed domain, serialize non-void cells
     call inputstr_compress

!-- READ DATA
!-- read ion and level data
     call ions_read_data(gas_nelem)  !ion and level data
!-- read bbxs data
     if(.not.in_nobbopac) call read_bbxs_data(gas_nelem)!bound-bound cross section data
!-- read bfxs data
     if(.not.in_nobfopac) call bfxs_read_data           !bound-free cross section data
!-- read ffxs data
     if(.not.in_noffopac) call ffxs_read_data           !free-free cross section data
!
     t1 = t_time()
     t_setup = t1-t0!}}}
  endif !impi
!-- broadcast init info from impi0 rank to all others
  call bcast_permanent !MPI
!-- domain-decompose input structure
  call scatter_inputstruct(in_ndim,icell1,ncell) !MPI

!--
!-- setup remaining modules
!==========================
!-- time step init
  tfirst = max(in_tsp_tfirst,in_tfirst*pc_day)
  tlast = max(in_tsp_tlast,in_tlast*pc_day)
  call timestepmod_init(in_nt,in_ntres,tfirst)
!-- constant time step, may be coded to loop if time step is not uniform
  dt = (tlast - tfirst)/in_nt

!-- wlgrid (before grid setup)
  call groupmod_init(in_ng,in_wldex,in_wlmin,in_wlmax)
  call fluxgrid_setup(in_flx_ndim,in_flx_wlmin,in_flx_wlmax)

!-- setup spatial grid
  call gridmod_init(lmpi0,grp_ng,in_igeom,in_ndim,str_nc,icell1,ncell, &
          str_lvoid,in_isvelocity)
  call grid_setup
!-- setup gas
  call gasmod_init(lmpi0,icell1,ncell,grp_ng)
  call gas_setup
!-- inputstr no longer needed
  call inputstr_dealloc

!-- particle init
  mpart = int(int(2,8)**in_trn_n2part/nmpi) !max number of particles
  mpart = max(mpart,in_prt_nmax/nmpi)
  call particlemod_init(mpart,in_isimcanlog, &
       in_isddmcanlog,in_tauddmc,in_taulump,in_tauvtime)
!-- create procedure pointers for the selected geometry
  call transportmod_init(in_igeom)
!-- source particles
  ns = int(int(2,8)**in_src_n2s/nmpi)
  ns = max(ns,in_ns/nmpi)
  nsinit = int(int(2,8)**in_src_n2sinit/nmpi)
  nsinit = max(nsinit,in_ns/nmpi)
  call sourcemod_init(ns,nsinit)

!-- allocate arrays of sizes retreived in bcast_permanent
  call ions_alloc_grndlev(gas_nelem,gas_ncell)  !ground state occupation numbers
  call particle_alloc(lmpi0,in_norestart,nmpi)

!-- initialize random number generator, use different seeds for each rank
  call rnd_init(in_nomp,impi)
!-- use extra stream for non-threaded code
  rnd_state = rnd_states(1)

!-- reading restart rand() count
  if(tsp_ntres>1.and..not.in_norestart) then
     call scatter_restart_data !MPI
  endif

!-- initial radiation energy
  call initialnumbers
!-- instantiating initial particles (if any)
  call initial_particles



!-- time step loop
!=================
  if(lmpi0) then
     msg = 'post setup:'
     write(6,*) 'memusg: ',msg,memusg()
!
     write(6,*)
     write(6,*) "starting time loop:"
     write(6,*) "===================="
  endif
!
  do it=in_ntres,tsp_nt
     t_timelin(1) = t_time() !timeline
!-- allow negative and zero it for temperature initialization purposes
     tsp_it = max(it,1)

!-- Update tsp_t etc
     call timestep_update(dt)  !dt is input here, any value can be passed
     call tau_update(tsp_t,tfirst,tlast) !updating prt_tauddmc and prt_taulump

!-- write timestep
     help = merge(tot_eerror,tot_erad,it>1)
     if(lmpi0) write(6,'(1x,a,i5,f8.3,"d",i10,1p,2e10.2)') 'timestep:', &
        it,tsp_t/pc_day,count(.not.prt_isvacant),help

!-- update all non-permanent variables
     call grid_update(tsp_t)
     call gas_update(it)

!-- source energy: gamma and material
     call sourceenergy(nmpi)

!-- grey gamma ray transport
     call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
     t_timelin(2) = t_time() !timeline
     if(in_srctype=='none' .and. .not.in_novolsrc) then
        call allgather_gammacap
        call particle_advance_gamgrey(nmpi)
        call allreduce_gammaenergy !MPI
!       grd_edep = grd_emitex !for testing: local deposition
     else
        grd_edep = 0d0
     endif

!-- gather from gas workers and broadcast to world ranks
     call bcast_nonpermanent !MPI
     t_timelin(3) = t_time() !timeline

!-- source energy
     call sourceenergy_misc(lmpi0)      !add gamma source energy and amplification-factor energy
     call sourceenergy_analytic(lmpi0)  !gas_emitex from analytic distribution

     call leakage_opacity       !IMC-DDMC albedo coefficients and DDMC leakage opacities
     call emission_probability  !emission probabilities for ep-group in each cell
     call allgather_leakage !MPI

     t_timelin(4) = t_time()    !timeline
     call sourcenumbers(nmpi)   !number of source prt_particles per cell

     if(src_nnew>0) then
        allocate(src_ivacant(src_nnew))
        call vacancies          !Storing vacant "prt_particles" indexes in ordered array "src_ivacant"
        call boundary_source    !properties of prt_particles on domain boundary
        call interior_source    !properties of prt_particles emitted in domain interior
        if(in_isvelocity) call source_transformdirection
        deallocate(src_ivacant)
     endif
     if(tsp_it<=tsp_ntres) where(.not.prt_isvacant) prt_particles%t = tsp_t !reset particle clocks

!-- advance particles
     t_timelin(5) = t_time() !timeline
     call particle_advance
     call reduce_tally !MPI  !collect particle results from all workers
     t_timelin(6) = t_time() !timeline

!-- print packet advance load-balancing info
     !if(lmpi0) write(6,'(1x,a,3(f9.2,"s"))') 'packets time(min|mean|max):',t_pckt_stat
     if(lmpi0) then
        call timereg(t_pcktmin,t_pckt_stat(1))
        call timereg(t_pcktmea,t_pckt_stat(2))
        call timereg(t_pcktmax,t_pckt_stat(3))
     endif

!-- collect data necessary for restart (tobe written by impi0)
     if(.not.in_norestart) call collect_restart_data !MPI

!-- update temperature
     call temperature_update
     call reduce_gastemp !MPI

!-- output
     if(lmpi0) then
!-- total energy startup values and energy conservation
        if(it<1) call totals_startup
        call totals_error !check energy (particle weight) is accounted

!-- write output
        if(it>0) call write_output
!-- restart writers
        if(.not.in_norestart .and. it>0) then
           call write_restart_file !temp
           call write_restart_randcount !rand() count
           call write_restart_particles !particle properties of current time step
        endif
     endif !impi

!-- write timestep timing to file
     if(it>0) call timing_timestep(impi)
     t_timelin(7) = t_time() !timeline
     t_timeline(:6) = t_timeline(:6) + (t_timelin(2:) - t_timelin(:6))
  enddo !tsp_it
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output.
  if(lmpi0) then
!
!-- print memory usage
     msg = 'post loop:'
     write(6,*)
     write(6,*) 'memusg: ',msg,memusg()

!-- print cpu timing usage
     t1 = t_time()
     t_all = t1 - t0
     call print_timing  !print timing results
     write(6,*)
     write(6,*) 'SuperNu finished'
     if(in_grabstdout) write(0,'(a,f8.2,"s")')'SuperNu finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This helps to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

end program supernu
