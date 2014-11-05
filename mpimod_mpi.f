      module mpimod
c     -------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: impi0 = 0 !the master rank
      integer :: impi !mpi rank
      integer :: nmpi !number of mpi tasks
      integer,private :: ierr
c
      contains
c
c
c
      subroutine bcast_permanent
c     --------------------------!{{{
      use inputparmod
      use ionsmod
      use ffxsmod
      use bfxsmod
      use bbxsmod
      use gasgridmod,nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use particlemod
      use timestepmod
      use fluxmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
************************************************************************
      integer :: n
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
c
c-- broadcast constants
c-- logical
      n = 9
      allocate(lsndvec(n))
      if(impi==impi0) lsndvec = (/in_isvelocity,in_puretran,
     &  prt_isimcanlog,prt_isddmcanlog,in_norestart,in_noeos,
     &  in_novolsrc,in_noreadstruct,in_ismodimc/)
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_isvelocity = lsndvec(1)
      in_puretran = lsndvec(2)
      prt_isimcanlog = lsndvec(3)
      prt_isddmcanlog = lsndvec(4)
      in_norestart = lsndvec(5)
      in_noeos = lsndvec(6)
      in_novolsrc = lsndvec(7)
      in_noreadstruct = lsndvec(8)
      in_ismodimc = lsndvec(9)
      deallocate(lsndvec)
c
c-- integer
      n = 21
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/in_igeom,
     &  in_ndim(1),in_ndim(2),in_ndim(3),gas_ng,prt_ns,
     &  prt_npartmax,in_nomp,tsp_nt,in_ntres,tsp_ntres,
     &  prt_ninit,prt_ninitnew,in_ng,in_nheav,
     &  ion_nion,ion_iionmax,bb_nline,
     &  flx_ng,flx_nmu,flx_nom/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_igeom     = isndvec(1) 
      in_ndim(1)   = isndvec(2)
      in_ndim(2)   = isndvec(3)
      in_ndim(3)   = isndvec(4)
      gas_ng       = isndvec(5)
      prt_ns       = isndvec(6)
      prt_npartmax = isndvec(7)
      in_nomp      = isndvec(8)
      tsp_nt       = isndvec(9)
      in_ntres     = isndvec(10)
      tsp_ntres    = isndvec(11)
      prt_ninit    = isndvec(12)
      prt_ninitnew = isndvec(13)
      in_ng        = isndvec(14)
      in_nheav     = isndvec(15)
      ion_nion     = isndvec(16)
      ion_iionmax  = isndvec(17)
      bb_nline     = isndvec(18)
      flx_ng       = isndvec(19)
      flx_nmu      = isndvec(20)
      flx_nom      = isndvec(21)
      deallocate(isndvec)
c
c-- real*8
      n = 21
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/prt_tauddmc,prt_taulump,
     &  tsp_t,tsp_dt,tsp_alpha,
     &  in_sigcoefs,in_sigtpwrs,in_sigrpwrs,
     &  in_sigcoef, in_sigtpwr, in_sigrpwr,
     &  in_suolpick1,in_ldisp1,in_ldisp2,in_theav,in_srcmax,
     &  in_consttemp,in_tempradinit,
     &  in_cvcoef,in_cvtpwr,in_cvrpwr/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_tauddmc  = sndvec(1)
      prt_taulump  = sndvec(2)
      tsp_t        = sndvec(3)
      tsp_dt       = sndvec(4)
      tsp_alpha    = sndvec(5)
      in_sigcoefs  = sndvec(6)
      in_sigtpwrs  = sndvec(7)
      in_sigrpwrs  = sndvec(8)
      in_sigcoef   = sndvec(9)
      in_sigtpwr   = sndvec(10)
      in_sigrpwr   = sndvec(11)
      in_suolpick1 = sndvec(12)
      in_ldisp1    = sndvec(13)
      in_ldisp2    = sndvec(14)
      in_theav     = sndvec(15)
      in_srcmax    = sndvec(16)
      in_consttemp = sndvec(17)
      in_tempradinit=sndvec(18)
      in_cvcoef    = sndvec(19)
      in_cvtpwr    = sndvec(20)
      in_cvrpwr    = sndvec(21)
      deallocate(sndvec)
c
c-- character
      call mpi_bcast(in_srctype,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(in_suol,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(in_opacanaltype,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(prt_tauvtime,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       allocate(gas_wl(gas_ng+1))
       allocate(flx_wl(flx_ng+1))
       allocate(flx_mu(flx_nmu+1))
       allocate(flx_om(flx_nom+1))
       if(bb_nline>0) allocate(bb_xs(bb_nline))
      endif
c
c-- broadcast data
      call mpi_bcast(gas_wl,gas_ng+1,MPI_REAL8,
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
      call mpi_bcast(bf_ph1,7*30*30,MPI_REAL,
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
      subroutine setup_domain_decomposition
c     -------------------------------------!{{{
************************************************************************
* placeholder
************************************************************************
c!}}}
      end subroutine setup_domain_decomposition
c
c
c
      subroutine scatter_inputstruct(ndim)
c     ------------------------------------!{{{
      use inputstrmod
      use gasgridmod
      implicit none
      integer,intent(in) :: ndim(3)
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
* this is doing bcasts only for now.
************************************************************************
      integer :: nx,ny,nz,n
c
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)

      call mpi_bcast(str_nabund,1,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      if(impi/=impi0) then
       allocate(str_xleft(nx+1))
       allocate(str_yleft(ny+1))
       allocate(str_zleft(nz+1))
       allocate(str_mass(nx,ny,nz))
       if(str_nabund>0) allocate(str_massfr(str_nabund,nx,ny,nz))
       if(str_nabund>0) allocate(str_iabund(str_nabund))
      endif
      call mpi_bcast(str_xleft,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_yleft,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(str_zleft,nz+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      n = nx*ny*nz
      call mpi_bcast(str_mass,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      if(str_nabund>0) then
       n = str_nabund * nx*ny*nz
       call mpi_bcast(str_massfr,n,MPI_REAL8,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
c
      if(str_nabund>0) then
       n = str_nabund
       call mpi_bcast(str_iabund,n,MPI_INTEGER,
     &   impi0,MPI_COMM_WORLD,ierr)
      endif
!}}}
      end subroutine scatter_inputstruct
c
c
c
      subroutine bcast_nonpermanent
c     ------------------------!{{{
      use gasgridmod,nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use particlemod
      use timestepmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
!c
!c-- integer
!      n = 3
!      allocate(isndvec(n))
!      if(impi==impi0) isndvec = (/prt_nnew,prt_nsurf,
!     & prt_nexsrc/)
!      call mpi_bcast(isndvec,n,MPI_INTEGER,
!     &  impi0,MPI_COMM_WORLD,ierr)
!c-- copy back
!      prt_nnew = isndvec(1)
!      prt_nsurf = isndvec(2)
!      prt_nexsrc = isndvec(3)
!      deallocate(isndvec)
!c
!c-- real*8
!      n = 4
!      allocate(sndvec(n))
!      if(impi==impi0) sndvec = (/tsp_t,tsp_dt,gas_esurf,
!     & gas_etot/)
!      call mpi_bcast(sndvec,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!c-- copy back
!      tsp_t = sndvec(1)
!      tsp_dt = sndvec(2)
!      gas_esurf = sndvec(3)
!      gas_etot = sndvec(4)
!      deallocate(sndvec)
!c
!c-- initial send of gas_eext
!      if(tsp_it==1) then
!         call mpi_bcast(gas_eext,1,MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
!      endif
!c
!      n = nx*ny*nz
!      call mpi_bcast(gas_nvol,n,MPI_INTEGER,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_nvolex,n,MPI_INTEGER,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_emit,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_emitex,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!c
!      call mpi_bcast(gas_siggrey,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!c
!      call mpi_bcast(gas_fcoef,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_sig,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_emitprob,n*gas_ng,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_opacleak,6*n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_cap,n*gas_ng,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!c
!c
!c-- domain decomposition (becomes mpi_gather at some point)
!      call mpi_bcast(dd_temp,n,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
c-- special case: domain replicated copy
      gas_temp = dd_temp
c!}}}
      end subroutine bcast_nonpermanent
c
c
c
      subroutine reduce_tally
c     -----------------------!{{{
      use gasgridmod,nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use timingmod
      use fluxmod
      implicit none
************************************************************************
* Reduce the results from particle_advance that are needed for the
* temperature correction.
************************************************************************
      integer :: n
      real*8,allocatable :: sndvec(:),rcvvec(:)
      integer :: isnd3(nx,ny,nz),isnd3f(flx_ng,flx_nmu,flx_nom)
      real*8 :: snd3(nx,ny,nz),snd3f(flx_ng,flx_nmu,flx_nom)
      real*8 :: help

c
c-- dim==0
      n = 5
      allocate(sndvec(n))
      allocate(rcvvec(n))
      sndvec = (/gas_erad,gas_eright,gas_eleft,gas_eext,gas_evelo/)
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(impi==impi0) then
       gas_erad = rcvvec(1)/dble(nmpi)
       gas_eright = rcvvec(2)/dble(nmpi)
       gas_eleft = rcvvec(3)/dble(nmpi)
       gas_eextav = rcvvec(4)/dble(nmpi)
       gas_eveloav = rcvvec(5)/dble(nmpi)
      else
       gas_erad = 0d0
       gas_eright = 0d0
       gas_eleft = 0d0
c-- rtw: can't copy back 0 to eext or evelo.
      endif !impi
      deallocate(sndvec)
      deallocate(rcvvec)
c
c-- dim==1
      n = flx_ng*flx_nmu*flx_nom
      isnd3f = flx_lumnum
      call mpi_reduce(isnd3f,flx_lumnum,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3f = flx_luminos
      call mpi_reduce(snd3f,flx_luminos,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      flx_luminos = flx_luminos/dble(nmpi)
c
      snd3f = flx_lumdev
      call mpi_reduce(snd3f,flx_lumdev,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      flx_lumdev = flx_lumdev/dble(nmpi)
c
c-- dim==3
      n = nx*ny*nz
      isnd3 = gas_numcensus
      call mpi_reduce(isnd3,gas_numcensus,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3 = gas_edep
      call mpi_allreduce(snd3,gas_edep,n,MPI_REAL8,MPI_SUM,
     &  MPI_COMM_WORLD,ierr)
      gas_edep = gas_edep/dble(nmpi)
c
      isnd3 = gas_methodswap
      call mpi_reduce(isnd3,gas_methodswap,n,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      snd3 = gas_eraddens
      call mpi_reduce(snd3,gas_eraddens,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      gas_eraddens = gas_eraddens/dble(nmpi)
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
c!}}}
      end subroutine reduce_tally
c
c
c
      subroutine reduce_gastemp
c     -------------------------------------!{{{
      use gasgridmod
************************************************************************
* placeholder
************************************************************************
      gas_temp = dd_temp
c!}}}
      end subroutine reduce_gastemp
c
c
c
      subroutine scatter_restart_data
c     -------------------------------!{{{
      use particlemod
************************************************************************
* scatter restart data from master rank to subordinate ranks.
* allows for restart at some time step, tsp_it.
************************************************************************
c-- helper variables
      integer :: isq
      real :: hlp
c
c-- scattering part vacancy
      call mpi_scatter(prt_tlyvacant,prt_npartmax,MPI_LOGICAL,
     &     prt_isvacant,prt_npartmax,MPI_LOGICAL,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part zone
      call mpi_scatter(prt_tlyzsrc,prt_npartmax,MPI_INTEGER,
     &     prt_particles%zsrc,prt_npartmax,MPI_INTEGER,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part transport index
      call mpi_scatter(prt_tlyrtsrc,prt_npartmax,MPI_INTEGER,
     &     prt_particles%rtsrc,prt_npartmax,MPI_INTEGER,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part position
      call mpi_scatter(prt_tlyrsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%rsrc,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part direction
      call mpi_scatter(prt_tlymusrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%musrc,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part time
      call mpi_scatter(prt_tlytsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%tsrc,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part energy
      call mpi_scatter(prt_tlyesrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%esrc,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part birth energy
      call mpi_scatter(prt_tlyebirth,prt_npartmax,MPI_REAL8,
     &     prt_particles%ebirth,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering part wavelength
      call mpi_scatter(prt_tlywlsrc,prt_npartmax,MPI_REAL8,
     &     prt_particles%wlsrc,prt_npartmax,MPI_REAL8,impi0,
     &     MPI_COMM_WORLD,ierr)
c
c-- scattering rand() count
      call mpi_scatter(prt_tlyrandarr,1,MPI_INTEGER,
     &     prt_tlyrand,1,MPI_INTEGER,impi0,MPI_COMM_WORLD,ierr)
c
c-- iterating to correct rand() count
      do isq = 1, prt_tlyrand-1
         hlp = rand()
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
      call mpi_gather(prt_particles%zsrc,prt_npartmax,MPI_INTEGER,
     &     prt_tlyzsrc,prt_npartmax,MPI_INTEGER,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part transport index
      call mpi_gather(prt_particles%zsrc,prt_npartmax,MPI_INTEGER,
     &     prt_tlyrtsrc,prt_npartmax,MPI_INTEGER,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part position
      call mpi_gather(prt_particles%rsrc,prt_npartmax,MPI_REAL8,
     &     prt_tlyrsrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part direction
      call mpi_gather(prt_particles%musrc,prt_npartmax,MPI_REAL8,
     &     prt_tlymusrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part time
      call mpi_gather(prt_particles%tsrc,prt_npartmax,MPI_REAL8,
     &     prt_tlytsrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part energy
      call mpi_gather(prt_particles%esrc,prt_npartmax,MPI_REAL8,
     &     prt_tlyesrc,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part birth energy
      call mpi_gather(prt_particles%ebirth,prt_npartmax,MPI_REAL8,
     &     prt_tlyebirth,prt_npartmax,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &     ierr)
c
c-- gathering part wavelength
      call mpi_gather(prt_particles%wlsrc,prt_npartmax,MPI_REAL8,
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
