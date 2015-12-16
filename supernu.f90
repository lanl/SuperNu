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
  use countersmod

  implicit none
!***********************************************************************
! TODO and wishlist:
!***********************************************************************
  integer :: ierr, it
  integer :: icell1, ncell !number of cells per rank (gas_ncell)
  integer :: iitflux,itflux=0
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
  call countersmod_init
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
!-- redirect stdout to file if selected
     call open_logfiles
!-- read runtime parameters
     call deprecate_inputpars
!-- parse and verify runtime parameters
     call parse_inputpars(nmpi)
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
  call provide_inputpars(nmpi)
!-- domain-decompose input structure
  call scatter_inputstruct(in_ndim,icell1,ncell) !MPI

!--
!-- setup remaining modules
!==========================
  call timestepmod_init

!-- wlgrid (before grid setup)
  call groupmod_init(in_wldex)
  call fluxgrid_setup

!-- setup spatial grid
  call gridmod_init(lmpi0,grp_ng,str_nc,str_lvoid,icell1,ncell)
  call grid_setup
!-- setup gas
  call gasmod_init(lmpi0,icell1,ncell,grp_ng)
  call gas_setup
!-- inputstr no longer needed
  call inputstr_dealloc

!-- create procedure pointers for the selected geometry
  call transportmod_init(in_igeom)
!-- source particles
  call sourcemod_init(nmpi)

!-- allocate arrays of sizes retreived in bcast_permanent
  call ions_alloc_grndlev(gas_nelem,gas_ncell)  !ground state occupation numbers
  call particle_alloc(lmpi0)

!-- initialize random number generator, use different seeds for each rank
  call rnd_init(in_nomp,impi)
!-- use extra stream for non-threaded code
  rnd_state = rnd_states(1)

!-- initial radiation energy
  if(src_ninit>0) then
     call allgather_initialrad !MPI
     call initialnumbers(nmpi)
!-- instantiating initial particles (if any)
     call initial_particles
  endif

!-- memory usage
  if(lmpi0) then
     msg = 'post setup:'
     write(6,*) 'memusg: ',msg,memusg()
  endif


!-- time step loop
!=================
  if(lmpi0) then
     write(6,*)
     write(6,*) "timestep loop:"
     write(6,*) "===================="
     write(6,'(1x,a5,a9,a10,4(a7,1x),a7)') 'it','t[day]','e_err','nsrc','ncens','nflux','nactiv','usage'
  endif
!
  do it=tsp_itrestart,tsp_nt
     t_timelin(1) = t_time() !timeline
!-- allow negative and zero it for temperature initialization purposes
     tsp_it = max(it,1)

!-- Update tsp_t etc
     call timestep_update
     call tau_update(tsp_t,tsp_tfirst,tsp_tlast) !updating trn_tauddmc and trn_taulump

!-- write timestep
     if(lmpi0) write(6,'(1x,i5,f8.3,"d")',advance='no') it,tsp_t/pc_day

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
!       grd_tally(1,:) = grd_emitex !for testing: local deposition
     else
        grd_tally(1,:) = 0d0 !edep
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
     call sourcenumbers(.false.)!number of source prt_particles

     if(src_nnew>0) then
        allocate(src_ivacant(src_nnew))
        call vacancies          !Storing vacant "prt_particles" indexes in ordered array "src_ivacant"
        call boundary_source    !properties of prt_particles on domain boundary
        call interior_source    !properties of prt_particles emitted in domain interior
        if(in_isvelocity) call source_transformdirection
        deallocate(src_ivacant)
     endif
!-- 0 particle clocks before and in the first time step
     if(tsp_it<=1) where(.not.prt_isvacant) prt_particles%t = tsp_t !reset particle clocks

!-- advance particles
     t_timelin(5) = t_time() !timeline
     call particle_advance
     call reduce_gridtally !MPI  !collect transport results from all workers

!-- print packet advance load-balancing info
     !if(lmpi0) write(6,'(1x,a,3(f9.2,"s"))') 'packets time(min|mean|max):',t_pckt_stat
     if(lmpi0) then
        call timereg(t_pcktmin,t_pckt_stat(1))
        call timereg(t_pcktmea,t_pckt_stat(2))
        call timereg(t_pcktmax,t_pckt_stat(3))
     endif
     t_timelin(6) = t_time() !timeline

!-- flux loop
     do iitflux=itflux+1,tsp_it
        if(it<1) call fluxtally(iitflux)  !drop flux packets in init-iterations
        if(it<1) exit
        if(grd_isvelocity .and. .not.flx_noobservertime .and. &
              tsp_tarr(iitflux+1) > tsp_t*(1-grd_rout/pc_c)) exit
!-- tally flux packets
        call fluxtally(iitflux)
!-- output
        call reduce_fluxtally !MPI  !collect flux results from all workers
        if(lmpi0) call output_flux(iitflux)
        itflux = iitflux
     enddo
     t_timelin(7) = t_time() !timeline

!-- update temperature
     call temperature_update
     call reduce_gastemp !MPI  !for output

!-- output
     if(lmpi0) then
!-- total energy startup values and energy conservation
        if(it<1) call totals_startup
        call totals_error !check energy (particle weight) is accounted

!-- write output
        if(it>0) call output_gamflux
        if(it>0) call output_grid

!-- write stdout
        if(ct_npactive(2)<1000000) then
           write(6,'(1p,e10.2,4(i7,1x),0p,f7.4)') tot_eerror, &
              ct_npcreate(2),(ct_npcensimc(2)+ct_npcensddmc(2)), &
              ct_npflux(2),ct_npactive(2), &
              ct_nnonvacant(2)/dble(prt_npartmax)
        else
           write(6,'(1p,e10.2,4(i7,"k"),0p,f7.4)') tot_eerror, &
              ct_npcreate(2)/1000,(ct_npcensimc(2)+ct_npcensddmc(2))/1000, &
              ct_npflux(2)/1000,ct_npactive(2)/1000, &
              ct_nnonvacant(2)/dble(prt_npartmax)
        endif
     endif !impi

!-- write timestep timing to file
     call timing_timestep(impi,it<=0)
     call counters_timestep(impi,it<=0)
     t_timelin(8) = t_time() !timeline
     t_timeline = t_timeline + (t_timelin(2:) - t_timelin(:7))
  enddo !tsp_it
!
!
!--
!-- FINISH UP:
!=============
  call mpi_barrier(MPI_COMM_WORLD,ierr) !MPI
!-- Print timing output
  if(lmpi0) then
!
!-- print memory usage
     msg = 'post loop:'
     write(6,*)
     write(6,*) 'memusg: ',msg,memusg()
     write(6,*) 'particle array usage:',ct_nnonvacant(4)/dble(prt_npartmax)

!-- print cpu timing usage
     t1 = t_time()
     t_all = t1 - t0
     call print_counters
     call print_timing
     write(6,*)
     write(6,*) 'SuperNu finished'
     if(in_grabstdout) write(0,'(a,f8.2,"s")')'SuperNu finished',t_all!repeat to stderr
  endif
!-- Clean up memory. (This helps to locate memory leaks)
  call dealloc_all
  call mpi_finalize(ierr) !MPI

end program supernu
