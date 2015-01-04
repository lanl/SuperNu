      module mpimod
c     -------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: impi0=0 !the master rank
      integer :: impi !mpi rank
      integer :: nmpi !number of mpi tasks
      integer,private :: ierr
c
      save
c
      contains
c
c
c
      subroutine bcast_permanent
c     --------------------------!{{{
      use inputparmod
      use inputstrmod
      use ionsmod
      use ffxsmod
      use bfxsmod
      use bbxsmod
      use gridmod
      use gasmod
      use groupmod
      use particlemod
      use timestepmod
      use fluxmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
************************************************************************
      integer :: i,n
      integer :: il,ii,ir,ic
      integer :: nx,ny,nz
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      character(4),allocatable :: csndvec(:)
c
c-- inputparmod variables
c========================
c-- create pointer arrays
      call inputpar_create_pointers(il,ii,ir,ic)
c
c-- broadcast logicals
      n = il
      allocate(lsndvec(n))
      forall(i=1:n) lsndvec(i) = in_l(i)%p
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_l(i)%p = lsndvec(i)
      deallocate(lsndvec)
c-- broadcast integers
      n = ii
      allocate(isndvec(n))
      forall(i=1:n) isndvec(i) = in_i(i)%p
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_i(i)%p = isndvec(i)
      deallocate(isndvec)
c-- broadcast real*8
      n = ir
      allocate(sndvec(n))
      forall(i=1:n) sndvec(i) = in_r(i)%p
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      forall(i=1:n) in_r(i)%p = sndvec(i)
      deallocate(sndvec)
c-- broadcast characters
      n = ic
      allocate(csndvec(n))
      forall(i=1:n) csndvec(i) = in_c(i)%p
      call mpi_bcast(csndvec,n*4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
!     forall(i=1:n) in_c(i)%p = csndvec(i) !avoid gfortran 4.6.3 compiler bug
      do i=1,n
       in_c(i)%p = csndvec(i)
      enddo
      deallocate(csndvec)
c
c
c-- everything else
c==================
c-- broadcast constants
c-- logical
      n = 2
      allocate(lsndvec(n))
      if(impi==impi0) lsndvec = (/prt_isimcanlog,prt_isddmcanlog/)
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_isimcanlog = lsndvec(1)
      prt_isddmcanlog = lsndvec(2)
      deallocate(lsndvec)
c
c-- integer
      n = 16
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/
     &  grp_ng,prt_ns,
     &  prt_npartmax,tsp_nt,tsp_ntres,
     &  prt_ninit,prt_ninitnew,
     &  ion_nion,ion_iionmax,bb_nline,
     &  flx_ng,flx_nmu,flx_nom,
     &  str_nc,str_ncp,str_nabund/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      grp_ng       = isndvec(1)
      prt_ns       = isndvec(2)
      prt_npartmax = isndvec(3)
      tsp_nt       = isndvec(4)
      tsp_ntres    = isndvec(5)
      prt_ninit    = isndvec(6)
      prt_ninitnew = isndvec(7)
      ion_nion     = isndvec(8)
      ion_iionmax  = isndvec(9)
      bb_nline     = isndvec(10)
      flx_ng       = isndvec(11)
      flx_nmu      = isndvec(12)
      flx_nom      = isndvec(13)
      str_nc       = isndvec(14)
      str_ncp      = isndvec(15)
      str_nabund   = isndvec(16)
      deallocate(isndvec)
c
c-- real*8
      n = 5
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/prt_tauddmc,prt_taulump,
     &  tsp_t,tsp_dt,tsp_alpha/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_tauddmc  = sndvec(1)
      prt_taulump  = sndvec(2)
      tsp_t        = sndvec(3)
      tsp_dt       = sndvec(4)
      tsp_alpha    = sndvec(5)
      deallocate(sndvec)
c
c-- character
      call mpi_bcast(prt_tauvtime,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c-- dimenstions
      nx = in_ndim(1)
      ny = in_ndim(2)
      nz = in_ndim(3)
c
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       allocate(grp_wl(grp_ng+1))
       allocate(grp_wlinv(grp_ng+1))
       allocate(flx_wl(flx_ng+1))
       allocate(flx_mu(flx_nmu+1))
       allocate(flx_om(flx_nom+1))
       if(bb_nline>0) allocate(bb_xs(bb_nline))
       allocate(str_xleft(nx+1))
       allocate(str_yleft(ny+1))
       allocate(str_zleft(nz+1))
       allocate(str_idcell(str_ncp))
       if(str_nabund>0) allocate(str_iabund(str_nabund))
      endif
c
c-- inputstr
      call mpi_bcast(str_xleft,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_yleft,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_zleft,nz+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_idcell,str_ncp,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      if(str_nabund>0) then
       call mpi_bcast(str_iabund,str_nabund,MPI_INTEGER,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c
c-- broadcast data
      call mpi_bcast(grp_wl,grp_ng+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(grp_wlinv,grp_ng+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(flx_wl,flx_ng+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(flx_mu,flx_nmu+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(flx_om,flx_nom+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- bound-bound
      if(bb_nline>0) then
       call mpi_bcast(bb_xs,sizeof(bb_xs),MPI_BYTE,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c-- bound-free
      call mpi_bcast(bf_ph1,6*7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(bf_ph2,7*30*30,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- free-free
      call mpi_bcast(ff_gff,ff_nu*ff_ngg,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call bcast_ions
c
c
      contains
c
      subroutine bcast_ions
c     ---------------------!{{{
      implicit none
************************************************************************
* broadcast the ions data structure
************************************************************************
      integer :: ii,iz,iion,n
      real*8 :: vec(1000)
      integer :: nion(gas_nelem)
      integer :: nlev(ion_nion)
      real*8 :: e(ion_nion)
c
c-- evaluate shape info
      if(impi==impi0) then
       iion = 0
       do iz=1,gas_nelem
        nion(iz) = ion_el(iz)%ni
        do ii=1,ion_el(iz)%ni
         iion = iion + 1
         nlev(iion) = ion_el(iz)%i(ii)%nlev
         e(iion) = ion_el(iz)%i(ii)%e
        enddo !ii
       enddo !iz
c-- sanity check
       if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
      endif
c
c-- bcast shape info and allocate
      call mpi_bcast(nion,gas_nelem,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nlev,ion_nion,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- allocate structure
      if(impi/=impi0) call ion_alloc_el(gas_nelem,nion,ion_nion,nlev)
c
c-- fill structure
      call mpi_bcast(e,ion_nion,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      iion = 0
      do iz=1,gas_nelem
       do ii=1,ion_el(iz)%ni
        iion = iion + 1
        n = nlev(iion)
c-- eion
        if(impi/=impi0) ion_el(iz)%i(ii)%e = e(iion)
c-- elev
        if(impi==impi0) vec(:n) = ion_el(iz)%i(ii)%elev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%elev = vec(:n)
c-- glev
        if(impi==impi0) vec(:n) = ion_el(iz)%i(ii)%glev
        call mpi_bcast(vec(1),n,MPI_REAL8,
     &    impi0,MPI_COMM_WORLD,ierr)
        if(impi/=impi0) ion_el(iz)%i(ii)%glev = vec(:n)
       enddo !ii
      enddo !iz
c-- sanity check
      if(iion/=ion_nion) stop "bcast_perm: ion_nion problem"
c!}}}
      end subroutine bcast_ions
c!}}}
      end subroutine bcast_permanent
c
c
c
      subroutine scatter_inputstruct(ndim,ncell)
c     ------------------------------------------!{{{
      use inputstrmod
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(out) :: ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      integer :: nx,ny,nz
c
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)
c
      ncell = str_ncp/nmpi
c
c-- allocate domain decomposed and domain compressed
      if(impi/=impi0) allocate(str_massdc(str_ncp))
      allocate(str_massdd(ncell))
      call mpi_scatter(str_massdc,ncell,MPI_REAL8,
     &  str_massdd,ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- mass fractions if available
      if(str_nabund>0) then
       if(impi/=impi0) allocate(str_massfrdc(str_nabund,str_ncp))
       allocate(str_massfrdd(str_nabund,ncell))
       call mpi_scatter(str_massfrdc,str_nabund*ncell,MPI_REAL8,
     &   str_massfrdd,str_nabund*ncell,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
!}}}
      end subroutine scatter_inputstruct
c
c
c
      subroutine allgather_gammacap
c     -----------------------------!{{{
      use gridmod
      use gasmod
      implicit none
************************************************************************
* gather gas_capgam to grd_capgrey
************************************************************************
      call mpi_allgather(gas_emitex,gas_ncell,MPI_REAL8,
     &  grd_emitex,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgather(gas_capgam,gas_ncell,MPI_REAL8,
     &  grd_capgrey,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)!}}}
      end subroutine allgather_gammacap
c
c
      subroutine allreduce_gammaenergy
c     ------------------------------------!{{{
      use gridmod,nx=>grd_nx,ny=>grd_ny,nz=>grd_nz
      use timingmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      real*8 :: snd(grd_ncp)
      real*8 :: t0,t1
c
      t0 = t_time()
c
      snd = grd_edep
      call mpi_allreduce(snd,grd_edep,grd_ncp,MPI_REAL8,MPI_SUM,
     &  MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpigamma, t1-t0)
c!}}}
      end subroutine allreduce_gammaenergy
c
c
      subroutine bcast_nonpermanent
c     -----------------------------!{{{
      use gridmod
      use gasmod
      use groupmod
      use totalsmod
      use particlemod
      use timestepmod
      use timingmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      real*8 :: t0,t1
      real*8 :: snd(grd_ncp)
c
      t0 = t_time()
c
c-- gather
      call mpi_allgather(gas_temp,gas_ncell,MPI_REAL8,
     &  grd_temp,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgather(gas_fcoef,gas_ncell,MPI_REAL8,
     &  grd_fcoef,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgather(gas_capgrey,gas_ncell,MPI_REAL8,
     &  grd_capgrey,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgather(gas_emit,gas_ncell,MPI_REAL8,
     &  grd_emit,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgather(gas_emitex,gas_ncell,MPI_REAL8,
     &  grd_emitex,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
c
      call mpi_allgather(gas_sig,gas_ncell,MPI_REAL8,
     &  grd_sig,gas_ncell,MPI_REAL8,
     &  MPI_COMM_WORLD,ierr)
      call mpi_allgather(gas_cap,grp_ng*gas_ncell,MPI_REAL,
     &  grd_cap,grp_ng*gas_ncell,MPI_REAL,
     &  MPI_COMM_WORLD,ierr)
c
c-- broadcast
      call mpi_bcast(tot_esurf,1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- allreduce
      snd = grd_eamp
      call mpi_allreduce(snd,grd_eamp,grd_ncp,MPI_REAL8,MPI_SUM,
     &  MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpibcast, t1-t0)
c!}}}
      end subroutine bcast_nonpermanent
c
c
c
      subroutine reduce_tally
c     -----------------------!{{{
      use gridmod,nx=>grd_nx,ny=>grd_ny,nz=>grd_nz
      use totalsmod
      use gasmod
      use timingmod
      use fluxmod
      implicit none
************************************************************************
* Reduce the results from particle_advance that are needed for the
* temperature correction.
************************************************************************
      integer :: n
      real*8,allocatable :: sndvec(:),rcvvec(:)
      integer :: isnd(grd_ncp)
      real*8 :: snd(grd_ncp)
      integer :: isnd3f(flx_ng,flx_nmu,flx_nom)
      real*8 :: snd3f(flx_ng,flx_nmu,flx_nom)
      integer :: isnd2f(flx_nmu,flx_nom)
      real*8 :: snd2f(flx_nmu,flx_nom)
      real*8 :: help
      real*8 :: t0,t1
c
      t0 = t_time()
c
c-- dim==0
      n = 6
      allocate(sndvec(n))
      allocate(rcvvec(n))
      sndvec = [tot_erad,tot_eout,tot_eext,tot_evelo,tot_emat,
     &  tot_eext0]
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(impi==impi0) then
         tot_erad = rcvvec(1)
         tot_eout = rcvvec(2)
         tot_eext = rcvvec(3)
         tot_evelo = rcvvec(4)
         tot_emat = rcvvec(5)
         tot_eext0 = rcvvec(6)
      else
c-- zero out cumulative values on all other ranks to avoid double counting.
         tot_eout = 0d0
         tot_eext = 0d0
         tot_evelo = 0d0
      endif
      deallocate(sndvec)
      deallocate(rcvvec)

c
c-- flux dim==2
      n = flx_nmu*flx_nom
      isnd2f = flx_gamlumnum
      call mpi_reduce(isnd2f,flx_gamlumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = flx_gamluminos
      call mpi_reduce(snd2f,flx_gamluminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd2f = flx_gamlumdev
      call mpi_reduce(snd2f,flx_gamlumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- flux dim==3
      n = flx_ng*flx_nmu*flx_nom
      isnd3f = flx_lumnum
      call mpi_reduce(isnd3f,flx_lumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_luminos
      call mpi_reduce(snd3f,flx_luminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_lumdev
      call mpi_reduce(snd3f,flx_lumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- dim==3
      n = grd_ncp
      isnd = grd_numcensus
      call mpi_reduce(isnd,grd_numcensus,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      isnd = grd_methodswap
      call mpi_reduce(isnd,grd_methodswap,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd = grd_edep
      call mpi_reduce(snd,grd_edep,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd = grd_eraddens
      call mpi_reduce(snd,grd_eraddens,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- scatter
      call mpi_scatter(grd_edep,gas_ncell,MPI_REAL8,
     &  gas_edep,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_scatter(grd_eraddens,gas_ncell,MPI_REAL8,
     &  gas_eraddens,gas_ncell,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- timing statistics
      help = t_pckt_stat(1)
      call mpi_reduce(help,t_pckt_stat(1),1,MPI_REAL8,MPI_MIN,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(2)/nmpi
      call mpi_reduce(help,t_pckt_stat(2),1,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      help = t_pckt_stat(3)
      call mpi_reduce(help,t_pckt_stat(3),1,MPI_REAL8,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      t1 = t_time()
      call timereg(t_mpireduc, t1-t0)
c!}}}
      end subroutine reduce_tally
c
c
c
      subroutine reduce_gastemp
c     -------------------------------------!{{{
      use gridmod
      use gasmod
************************************************************************
* for output
************************************************************************
      call mpi_gather(gas_temp,gas_ncell,MPI_REAL8,
     &   grd_temp,gas_ncell,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
c!}}}
      end subroutine reduce_gastemp
c
c
c
      subroutine scatter_restart_data
c     -------------------------------!{{{
      use particlemod
      use randommod
************************************************************************
* scatter restart data from master rank to subordinate ranks.
* allows for restart at some time step, tsp_it.
************************************************************************
c-- helper variables
      integer :: isq
      real*8 :: hlp
c
c-- scattering part vacancy
      call mpi_scatter(prt_tlyvacant,prt_npartmax,MPI_LOGICAL,
     &     prt_isvacant,prt_npartmax,MPI_LOGICAL,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part zone
      call mpi_scatter(prt_tlyzsrc,prt_npartmax,MPI_INTEGER,
     &     prt_particles%ix,prt_npartmax,MPI_INTEGER,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part transport index
      call mpi_scatter(prt_tlyrtsrc,prt_npartmax,MPI_INTEGER,
     &     prt_particles%itype,prt_npartmax,MPI_INTEGER,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part position
      call mpi_scatter(prt_tlyrsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%x,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part direction
      call mpi_scatter(prt_tlymusrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%mu,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part time
      call mpi_scatter(prt_tlytsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%t,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part energy
      call mpi_scatter(prt_tlyesrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%e,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part birth energy
      call mpi_scatter(prt_tlyebirth,prt_npartmax,MPI_REAL8,
     &     prt_particles%e0,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part wavelength
      call mpi_scatter(prt_tlywlsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%wl,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering rand() count
      call mpi_scatter(prt_tlyrandarr,1,MPI_INTEGER,
     &     prt_tlyrand,1,MPI_INTEGER,impi0,MPI_COMM_WORLD,ierr)
c
c-- iterating to correct rand() count
      do isq = 1, prt_tlyrand-1
         hlp = rnd_r(rnd_state)
      enddo
c
c-- deallocations



c!}}}
      end subroutine scatter_restart_data
c
c
c
      subroutine collect_restart_data
c     -------------------------------!{{{
      use particlemod
************************************************************************
* send particle array info and number of rand calls to master rank.
* allows for restart at some time step, tsp_it.
* Files written here to avoid too many allocations of large particle
* arrays.
************************************************************************
c
c-- gathering part vacancy
      call mpi_gather(prt_isvacant,prt_npartmax,MPI_LOGICAL,
     &     prt_tlyvacant,prt_npartmax,MPI_LOGICAL,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part zone
      call mpi_gather(prt_particles%ix,prt_npartmax,MPI_INTEGER,
     &     prt_tlyzsrc,prt_npartmax,MPI_INTEGER,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part transport index
      call mpi_gather(prt_particles%ix,prt_npartmax,MPI_INTEGER,
     &     prt_tlyrtsrc,prt_npartmax,MPI_INTEGER,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part position
      call mpi_gather(prt_particles%x,prt_npartmax,MPI_REAL8,
     &     prt_tlyrsrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part direction
      call mpi_gather(prt_particles%mu,prt_npartmax,MPI_REAL8,
     &     prt_tlymusrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part time
      call mpi_gather(prt_particles%t,prt_npartmax,MPI_REAL8,
     &     prt_tlytsrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part energy
      call mpi_gather(prt_particles%e,prt_npartmax,MPI_REAL8,
     &     prt_tlyesrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part birth energy
      call mpi_gather(prt_particles%e0,prt_npartmax,MPI_REAL8,
     &     prt_tlyebirth,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part wavelength
      call mpi_gather(prt_particles%wl,prt_npartmax,MPI_REAL8,
     &     prt_tlywlsrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c
c====
c
c-- gathering rand() counts
      call mpi_gather(prt_tlyrand,1,MPI_INTEGER,prt_tlyrandarr,1,
     &     MPI_INTEGER,impi0,MPI_COMM_WORLD,ierr)
c
c-- deallocations


c!}}}
      end subroutine collect_restart_data    
c
      end module mpimod
