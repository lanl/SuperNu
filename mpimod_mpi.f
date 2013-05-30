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
      subroutine bcast_persistent
cc     --------------------------!{{{
c      use inputparmod
c      use timestepmod, only:tsp_nt
c      use gasgridmod
c      implicit none
c************************************************************************
c* Broadcast the data that does not evolve over time (or temperature).
c* Also once the constants are broadcasted, all allocatable arrays are
c* allocated.
c*
c* data:
c* - rgrid
c* - gas_wl
c* constants:
c* - in_ncostf,in_nphif,in_nwlf,in_wlmin,in_wlmax
c* - gas_ng,in_epsline,in_niwlem
c* - in_nr,rg_ncr,gas_nr,gas_npacket
c* - in_nomp,tsp_nt,in_ntc
c* - gas_xi2beta,gas_dxwin
c* - flx_wlhelp,flx_wlminlg
c*
c* Also scatter the rng_offset_all to the mpi ranks.
c************************************************************************
c      integer :: i
c      integer,allocatable :: isndvec(:)
c      real*8,allocatable :: sndvec(:)
cc
cc-- broadcast integer constants
c      allocate(isndvec(12))
c      if(impi==impi0) isndvec = (/in_ncostf,in_nphif,in_nwlf,gas_ng,
c     &  in_niwlem,in_nr,rg_ncr,gas_nr,gas_npacket,in_nomp,tsp_nt,
c     &  in_ntc/)
c      call mpi_bcast(isndvec,12,MPI_INTEGER,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc-- copy back
c      in_ncostf = isndvec(1)
c      in_nphif = isndvec(2)
c      in_nwlf = isndvec(3)
c      gas_ng = isndvec(4)
c      in_niwlem = isndvec(5)
c      in_nr = isndvec(6)
c      rg_ncr = isndvec(7)
c      gas_nr = isndvec(8)
c      gas_npacket = isndvec(9)
c      in_nomp = isndvec(10)
c      tsp_nt = isndvec(11)
c      in_ntc = isndvec(12)
c      deallocate(isndvec)
cc
cc-- broadcast real*8 constants
c      allocate(sndvec(7))
c      if(impi==impi0) then
c       sndvec(1:3) = (/in_wlmin,in_wlmax,in_epsline/)
c       sndvec(4:7) = (/gas_xi2beta,gas_dxwin,flx_wlhelp,flx_wlminlg/)
c      endif
c      call mpi_bcast(sndvec,7,MPI_REAL8,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc-- copy back
c      in_wlmin = sndvec(1)
c      in_wlmax = sndvec(2)
c      in_epsline = sndvec(3)
c      gas_xi2beta = sndvec(4)
c      gas_dxwin = sndvec(5)
c      flx_wlhelp = sndvec(6)
c      flx_wlminlg = sndvec(7)
c      deallocate(sndvec)
cc
cc$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
cc
cc-- allocate all arrays. These are deallocated in dealloc_all.f
c      if(impi/=impi0) then
c       allocate(rgrid(rg_ncr))
c       allocate(gas_vals(gas_nr))
c       allocate(gas_wl(gas_ng))
c       allocate(gas_cap(gas_nr,gas_ng))
c       allocate(flx_flux(in_nwlf,in_ncostf,in_nphif))
c       allocate(flx_iflux(in_nwlf,in_ncostf,in_nphif))
c      endif
c      allocate(rng_offset(0:in_nomp-1))
cc
cc-- broadcast data
c      call mpi_bcast(gas_wl,gas_ng,MPI_REAL8,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc
cc WARNING, sizeof may not work on hetrogeneous clusters!!!
c      call mpi_bcast(rgrid,sizeof(rgrid),MPI_BYTE,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc
cc-- scatter rng_offset
c      i = sizeof(rng_offset)
cc-- verify sizes
c      if(impi==impi0) then
c       if(sizeof(rng_offset_all) /= i*nmpi) stop
c     &   'bcast_perm_mpi: rng_offset_all size wrong'
c      endif
cc
cc-- allocate dummy source on non-master ranks
c!     if(impi/=impi0) allocate(rng_offset_all(nmpi*in_nomp)) !dummy, needed in debug mode
c      if(impi/=impi0) allocate(rng_offset_all(1)) !dummy, needed for mpi+ifort in debug mode
cc
c      call mpi_scatter(rng_offset_all,i,MPI_BYTE,
c     &  rng_offset,i,MPI_BYTE,
c     &  impi0,MPI_COMM_WORLD,ierr)
c      deallocate(rng_offset_all)
c!     write(6,'(5i13)') (impi*in_nomp+i,rng_offset(i),i=0,in_nomp-1)
cc!}}}
      end subroutine bcast_persistent
c
c
c
      subroutine bcast_nonpersistent
cc     ------------------------!{{{
c      use inputparmod
c      use gasgridmod
c      use timestepmod
c      implicit none
c************************************************************************
c* Broadcast the data that changes with time/temperature.
c* constants:
c* - tim_itc
c* data:
c* - gas_vals
c* - gas_cap
c************************************************************************
c      call mpi_bcast(tim_itc,1,MPI_INTEGER,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc
c      call mpi_bcast(gas_cap,gas_ng*gas_nr,MPI_REAL,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc
cc WARNING, this may not work on heterogeneous clusters!!!
c      call mpi_bcast(gas_vals,sizeof(gas_vals),MPI_BYTE,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc!}}}
      end subroutine bcast_nonpersistent
c
c
c
      subroutine reduce_tally
cc     -----------------------!{{{
c      use gasgridmod
c      use timestepmod
c      use timingmod
c      implicit none
c************************************************************************
c* Reduce the results from the packet transport that are needed for the
c* temperature correction.
c* - gas_vals%enabs_e
c* - gas_vals%enabs_c
c* - t_pckt_stat !min,mean,max
c************************************************************************
c      real*8 :: sndvec(gas_nr),rcvvec(gas_nr)
c      real*8 :: help
cc
c      sndvec = gas_vals(:)%enabs_e
c      call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
c     &  impi0,MPI_COMM_WORLD,ierr)
c      gas_vals(:)%enabs_e = rcvvec
cc
c      sndvec = gas_vals(:)%enabs_c
c      call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
c     &  impi0,MPI_COMM_WORLD,ierr)
c      gas_vals(:)%enabs_c = rcvvec
cc
c      if(tim_itc>0) then
c       sndvec = gas_vals%enocoll
c       call mpi_reduce(sndvec,rcvvec,gas_nr,MPI_REAL8,MPI_SUM,
c     &   impi0,MPI_COMM_WORLD,ierr)
c       gas_vals%enocoll = rcvvec
c      endif
cc
c      help = t_pckt_stat(1)
c      call mpi_reduce(help,t_pckt_stat(1),1,MPI_REAL8,MPI_MIN,
c     &  impi0,MPI_COMM_WORLD,ierr)
c      help = t_pckt_stat(2)/nmpi
c      call mpi_reduce(help,t_pckt_stat(2),1,MPI_REAL8,MPI_SUM,
c     &  impi0,MPI_COMM_WORLD,ierr)
c      help = t_pckt_stat(3)
c      call mpi_reduce(help,t_pckt_stat(3),1,MPI_REAL8,MPI_MAX,
c     &  impi0,MPI_COMM_WORLD,ierr)
cc!}}}
      end subroutine reduce_tally
c
c
      end module mpimod
