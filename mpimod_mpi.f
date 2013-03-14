      module mpimod
c     -------------
      implicit none
      include 'mpif.h'
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
c     --------------------------
      use inputparmod
      use timestepmod, only:tsp_nt
      use gasgridmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
*
* data:
* - rgrid
* - gas_wl
* constants:
* - in_ncostf,in_nphif,in_nwlf,in_wlmin,in_wlmax
* - gas_ng,in_epsline,in_niwlem
* - in_nr,rg_ncr,gas_nr,gas_npacket
* - in_nomp,tsp_nt,in_ntc
* - gas_xi2beta,gas_dxwin
* - flx_wlhelp,flx_wlminlg
*
* Also scatter the rng_offset_all to the mpi ranks.
************************************************************************
      integer :: i
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
c
c-- broadcast integer constants
      allocate(isndvec(12))
      if(impi==impi0) isndvec = (/in_ncostf,in_nphif,in_nwlf,gas_ng,
     &  in_niwlem,in_nr,rg_ncr,gas_nr,gas_npacket,in_nomp,tsp_nt,
     &  in_ntc/)
      call mpi_bcast(isndvec,12,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_ncostf = isndvec(1)
      in_nphif = isndvec(2)
      in_nwlf = isndvec(3)
      gas_ng = isndvec(4)
      in_niwlem = isndvec(5)
      in_nr = isndvec(6)
      rg_ncr = isndvec(7)
      gas_nr = isndvec(8)
      gas_npacket = isndvec(9)
      in_nomp = isndvec(10)
      tsp_nt = isndvec(11)
      in_ntc = isndvec(12)
      deallocate(isndvec)
c
c-- broadcast real*8 constants
      allocate(sndvec(7))
      if(impi==impi0) then
       sndvec(1:3) = (/in_wlmin,in_wlmax,in_epsline/)
       sndvec(4:7) = (/gas_xi2beta,gas_dxwin,flx_wlhelp,flx_wlminlg/)
      endif
      call mpi_bcast(sndvec,7,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_wlmin = sndvec(1)
      in_wlmax = sndvec(2)
      in_epsline = sndvec(3)
      gas_xi2beta = sndvec(4)
      gas_dxwin = sndvec(5)
      flx_wlhelp = sndvec(6)
      flx_wlminlg = sndvec(7)
      deallocate(sndvec)
c
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       allocate(rgrid(rg_ncr))
       allocate(gas_vals(gas_nr))
       allocate(gas_wl(gas_ng))
       allocate(gas_cap(gas_nr,gas_ng))
       allocate(flx_flux(in_nwlf,in_ncostf,in_nphif))
       allocate(flx_iflux(in_nwlf,in_ncostf,in_nphif))
      endif
      allocate(rng_offset(0:in_nomp-1))
c
c-- broadcast data
      call mpi_bcast(gas_wl,gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c WARNING, sizeof may not work on hetrogeneous clusters!!!
      call mpi_bcast(rgrid,sizeof(rgrid),MPI_BYTE,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- scatter rng_offset
      i = sizeof(rng_offset)
c-- verify sizes
      if(impi==impi0) then
       if(sizeof(rng_offset_all) /= i*nmpi) stop
     &   'bcast_perm_mpi: rng_offset_all size wrong'
      endif
c
c-- allocate dummy source on non-master ranks
!     if(impi/=impi0) allocate(rng_offset_all(nmpi*in_nomp)) !dummy, needed in debug mode
      if(impi/=impi0) allocate(rng_offset_all(1)) !dummy, needed for mpi+ifort in debug mode
c
      call mpi_scatter(rng_offset_all,i,MPI_BYTE,
     &  rng_offset,i,MPI_BYTE,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(rng_offset_all)
!     write(*,'(5i13)') (impi*in_nomp+i,rng_offset(i),i=0,in_nomp-1)
c
      end subroutine bcast_permanent
c
c
c
      subroutine bcast_mutable
c     ------------------------
      use inputparmod
      use gasgridmod
      use timestepmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
* constants:
* - tim_itc
* data:
* - gas_vals
* - gas_cap
************************************************************************
      call mpi_bcast(tim_itc,1,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_bcast(gas_cap,gas_ng*gas_nr,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c WARNING, this may not work on heterogeneous clusters!!!
      call mpi_bcast(gas_vals,sizeof(gas_vals),MPI_BYTE,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine bcast_mutable
c
c
c
      subroutine reduce_enabs
c     -----------------------
      use gasgridmod
      use timestepmod
      use timingmod
      implicit none
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - gas_vals%enabs_e
* - gas_vals%enabs_c
* - t_pckt_stat !min,mean,max
************************************************************************
      real*8 :: sndvec(gas_nr),rcvvec(gas_nr)
      real*8 :: help
c
      sndvec = gas_vals(:)%enabs_e
      call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      gas_vals(:)%enabs_e = rcvvec
c
      sndvec = gas_vals(:)%enabs_c
      call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      gas_vals(:)%enabs_c = rcvvec
c
      if(tim_itc>0) then
       sndvec = gas_vals%enocoll
       call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
     &   impi0,MPI_COMM_WORLD,ierr)
       gas_vals%enocoll = rcvvec
      endif
c
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
      end subroutine reduce_enabs
c
c
      end module mpimod
