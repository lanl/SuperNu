      module mpimod
c     -------------
      implicit none
      INCLUDE 'mpif.h'
c
      integer,parameter :: impi0 = 0 !the master rank
      integer :: impi !mpi rank
      integer :: nmpi !number of mpi tasks
      integer,private :: ierr
      integer :: status(MPI_STATUS_SIZE)
c
      contains
c
c
c
      subroutine bcast_permanent
c     --------------------------!{{{
      use inputparmod
      use gasgridmod,nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use particlemod
      use timestepmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
*
* arrays:
* -------
* integer :: gas_nvolinit(nx,ny,nz)
* real*8 :: gas_xarr(nx+1)
* real*8 :: gas_yarr(ny+1)
* real*8 :: gas_zarr(nz+1)
* real*8 :: gas_dxarr(nx)
* real*8 :: gas_dyarr(ny)
* real*8 :: gas_dzarr(nz)
* real*8 :: gas_evolinit(nx,ny,nz)
* real*8 :: gas_wl(gas_ng+1)
*
* scalars:
*----------
*-- logical
* logical :: in_isvelocity
* logical :: in_puretran
* logical :: prt_isimcanlog
* logical :: prt_isddmcanlog
* logical :: in_norestart
*-- integer
* integer :: gas_nx
* integer :: gas_ny
* integer :: gas_nz
* integer :: gas_ng
* integer :: prt_npartmax
* integer :: in_nomp
* integer :: tsp_nt
* integer :: tsp_ntres
* integer :: in_seed
* integer :: prt_ninitnew
*-- real*8
* real*8 :: prt_tauddmc
* real*8 :: prt_taulump
* real*8 :: tsp_t
*-- character
* character(4) :: gas_srctype
* character(4) :: prt_tauvtime
*
************************************************************************
      integer :: n
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      character*4,allocatable :: csndvec(:)
      
c
c-- broadcast constants
c-- logical
      n = 5
      allocate(lsndvec(n))
      if(impi==impi0) lsndvec = (/gas_isvelocity,in_puretran,
     &  prt_isimcanlog,prt_isddmcanlog,in_norestart/)
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      gas_isvelocity = lsndvec(1)
      in_puretran = lsndvec(2)
      prt_isimcanlog = lsndvec(3)
      prt_isddmcanlog = lsndvec(4)
      in_norestart = lsndvec(5)
      deallocate(lsndvec)
c
c-- integer
      n = 11
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/gas_nx,gas_ny,gas_nz,gas_ng,
     &  prt_npartmax,in_nomp,tsp_nt,tsp_ntres,in_seed,prt_ninitnew,
     &     in_ng/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      gas_nx = isndvec(1)
      gas_ny = isndvec(2)
      gas_nz = isndvec(3)
      gas_ng = isndvec(4)
      prt_npartmax = isndvec(5)
      in_nomp = isndvec(6)
      tsp_nt = isndvec(7)
      tsp_ntres = isndvec(8)
      in_seed = isndvec(9)
      prt_ninitnew = isndvec(10)
      in_ng = isndvec(11)
      deallocate(isndvec)
c
c-- real*8
      n = 3
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/prt_tauddmc,prt_taulump,tsp_t/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_tauddmc = sndvec(1)
      prt_taulump = sndvec(2)
      tsp_t = sndvec(3)
      deallocate(sndvec)
c
c-- character
      call mpi_bcast(gas_srctype,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacanaltype,4,MPI_CHARACTER,
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
       allocate(gas_nvolinit(nx,ny,nz))
       allocate(gas_xarr(nx+1))
       allocate(gas_yarr(ny+1))
       allocate(gas_zarr(nz+1))
       allocate(gas_dxarr(nx))
       allocate(gas_dyarr(ny))
       allocate(gas_dxarr(nz))
       allocate(gas_evolinit(nx,ny,nz))
       allocate(gas_wl(gas_ng+1))
       !prt_done = .false.
c-- allocating particle array for helper ranks
       allocate(prt_particles(prt_npartmax))
       prt_particles%isvacant = .true.
c--
      endif
      !if(impi==impi0) then
      !   prt_particles%isvacant=.true.
      !endif
c
c-- broadcast data
      call mpi_bcast(gas_nvolinit,nx*ny*nz,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_xarr,nx+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_yarr,ny+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_dzarr,nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_evolinit,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_wl,gas_ng+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c!}}}
      end subroutine bcast_permanent
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
*-- scalars:
*-- real
* real*8 :: tsp_t
* real*8 :: tsp_dt
* real*8 :: gas_esurf
* real*8 :: gas_etot
* real*8 :: gas_eext
*-- integer
* integer :: prt_nnew
* integer :: prt_nsurf
* integer :: prt_nexsrc
*--
*
*-- arrays:
*-- real
* real*8 :: gas_temp(nx,ny,nz)
* real*8 :: gas_nvol(nx,ny,nz)
* real*8 :: gas_nvolex(nx,ny,nz)
* real*8 :: gas_emit(nx,ny,nz)
* real*8 :: gas_emitex(nx,ny,nz)
* real*8 :: gas_fcoef(nx,ny,nz)
* real*8 :: gas_sig(nx,ny,nz)
* real*8 :: gas_emitprob(gas_ng,nx,ny,nz)
* real*8 :: gas_opacleak(2,nx,ny,nz)
* real*8 :: gas_cap(gas_ng,nx,ny,nz)
* real*8 :: gas_wl(gas_ng+1)
*-- integer
* integer :: gas_methodswap(nx,ny,nz)
*
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      real*8,allocatable :: sndmat(:,:)

c-- variables to be reduced -----------------------------------
c-- dim==1,2
      if(impi/=impi0 .and. tsp_it==tsp_ntres) then
         allocate(gas_numcensus(nx,ny,nz))
         allocate(gas_edep(nx,ny,nz))
         allocate(gas_eraddens(nx,ny,nz))
         allocate(gas_luminos(gas_ng))
         allocate(gas_lumdev(gas_ng))
         allocate(gas_lumnum(gas_ng))
         allocate(gas_methodswap(nx,ny,nz))
      endif
!      call mpi_bcast(gas_edep,nx*ny*nz,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_numcensus,nx*ny*nz,MPI_INTEGER,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_eraddens,nx*ny*nz,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
c--------------------------------------------------------------
c
c-- integer
      n = 3
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/prt_nnew,prt_nsurf,
     & prt_nexsrc/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_nnew = isndvec(1)
      prt_nsurf = isndvec(2)
      prt_nexsrc = isndvec(3)
      deallocate(isndvec)
c
c-- real*8
      n = 4
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/tsp_t,tsp_dt,gas_esurf,
     & gas_etot/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      tsp_t = sndvec(1)
      tsp_dt = sndvec(2)
      gas_esurf = sndvec(3)
      gas_etot = sndvec(4)
      deallocate(sndvec)
c
c-- initial send of gas_eext
      if(tsp_it==1) then
         n=1
         allocate(sndvec(n))
         sndvec=(/gas_eext/)
         call mpi_bcast(sndvec,n,MPI_REAL8,
     &        impi0,MPI_COMM_WORLD,ierr)
c-- copy back
         gas_eext = sndvec(1)
         deallocate(sndvec)
      endif
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0 .and. tsp_it==tsp_ntres) then
       allocate(gas_temp(nx,ny,nz))
       allocate(gas_nvol(nx,ny,nz))
       allocate(gas_nvolex(nx,ny,nz))
       allocate(gas_emit(nx,ny,nz))
       allocate(gas_emitex(nx,ny,nz))
c
       allocate(gas_fcoef(nx,ny,nz))
       allocate(gas_sig(nx,ny,nz))
       allocate(gas_siggrey(nx,ny,nz))
       allocate(gas_emitprob(gas_ng,nx,ny,nz))
       allocate(gas_opacleak(2,nx,ny,nz))
       allocate(gas_cap(gas_ng,nx,ny,nz))
!       allocate(gas_wl(gas_ng+1))
c
      endif
c
c-- broadcasting particle array
C$$$      call mpi_bcast(prt_particles%isvacant,prt_npartmax,
C$$$     &  MPI_LOGICAL,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%zsrc,prt_npartmax,
C$$$     &  MPI_INTEGER,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%rtsrc,prt_npartmax,
C$$$     &  MPI_INTEGER,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%rsrc,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%musrc,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%tsrc,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%esrc,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%ebirth,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
C$$$      call mpi_bcast(prt_particles%wlsrc,prt_npartmax,
C$$$     &  MPI_REAL8,impi0,MPI_COMM_WORLD,ierr)
c--
c
      call mpi_bcast(gas_temp,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_nvol,nx*ny*nz,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_nvolex,nx*ny*nz,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emit,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emitex,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_bcast(gas_siggrey,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_bcast(gas_fcoef,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_sig,nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emitprob,nx*ny*nz*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacleak,2*nx*ny*nz,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_cap,nx*ny*nz*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_wl,gas_ng+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c!}}}
      end subroutine bcast_nonpermanent
c
c
c
      subroutine reduce_tally
c     -----------------------!{{{
      use gasgridmod,nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use timingmod
      implicit none
************************************************************************
* Reduce the results from particle_advance that are needed for the
* temperature correction.
* - t_pckt_stat !min,mean,max
* 
*-- dim==0
* real*8 :: gas_erad
* real*8 :: gas_eright
* real*8 :: gas_eleft
*-- dim==1
* integer :: gas_lumnum(ng)
* real*8 :: gas_luminos(ng)
* real*8 :: gas_lumdev(ng)
*-- dim==3
* integer :: gas_numcensus(nx,ny,nz)
* real*8 :: gas_edep(nx,ny,nz)
* real*8 :: gas_eraddens(nx,ny,nz)
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:),rcvvec(:)
      integer :: isndvec3(nx,ny,nz)
      real*8 :: sndvec3(nx,ny,nz)
      real*8 :: help

c
c-- dim==0
      n = 5
      allocate(sndvec(n))
      allocate(rcvvec(n))
      !if(impi==impi0) 
      sndvec = (/gas_erad,gas_eright,gas_eleft,gas_eext,gas_evelo/)
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(impi==0) then
       gas_erad = rcvvec(1)
       gas_eright = rcvvec(2)
       gas_eleft = rcvvec(3)
       gas_eextav = rcvvec(4)
       gas_eveloav = rcvvec(5)
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
      allocate(isndvec(gas_ng))
      isndvec = gas_lumnum
      call mpi_reduce(isndvec,gas_lumnum,gas_ng,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(isndvec)
c
      allocate(sndvec(gas_ng))
      sndvec = gas_luminos
      call mpi_reduce(sndvec,gas_luminos,gas_ng,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(sndvec)
c
      allocate(sndvec(gas_ng))
      sndvec = gas_lumdev
      call mpi_reduce(sndvec,gas_lumdev,gas_ng,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(sndvec)
c
c-- dim==3
      isndvec3 = gas_numcensus
      call mpi_reduce(isndved3d,gas_numcensus,nx*ny*nz,
     &  MPI_INTEGER,MPI_SUM,impi0,MPI_COMM_WORLD,ierr)
c
      sndvec3 = gas_edep
      call mpi_reduce(sndvec3,gas_edep,nx*ny*nz,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      isndvec3 = gas_methodswap
      call mpi_reduce(isndvec3,gas_methodswap,nx*ny*nz,
     &  MPI_INTEGER,MPI_SUM,impi0,MPI_COMM_WORLD,ierr)
c
      n = nx*ny*nz
      sndvec3 = gas_eraddens
      call mpi_reduce(sndvec3,gas_eraddens,n,MPI_REAL8,MPI_SUM,
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
c!}}}
      end subroutine reduce_tally
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
     &     prt_particles%isvacant,prt_npartmax,MPI_LOGICAL,impi0,
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
      call mpi_gather(prt_particles%isvacant,prt_npartmax,MPI_LOGICAL,
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
