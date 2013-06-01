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
c     --------------------------!{{{
      use inputparmod
      use gasgridmod
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
* real*8 :: gas_rarr(gas_nr+1)
* real*8 :: gas_drarr(gas_nr)
* real*8 :: gas_curvcent(gas_nr)
*
* scalars:
*----------
*-- logical
* logical :: in_isvelocity
* logical :: in_puretran
* logical :: gas_isshell
* logical :: prt_isimcanlog
* logical :: prt_isddmcanlog
*-- integer
* integer :: gas_nr
* integer :: gas_ng
* integer :: gas_velno
* integer :: gas_velyes
* integer :: prt_npartmax
* integer :: in_nomp
* integer :: tsp_nt
*-- real*8
* real*8 :: prt_tauddmc
*
************************************************************************
      integer :: n
      logical,allocatable :: lsndvec(:)
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
c
c-- broadcast constants
c-- logical
      n = 5
      allocate(lsndvec(n))
      if(impi==impi0) lsndvec = (/in_isvelocity,in_puretran,gas_isshell,
     &  prt_isimcanlog,prt_isddmcanlog/)
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_isvelocity = lsndvec(1)
      in_puretran = lsndvec(2)
      gas_isshell = lsndvec(3)
      prt_isimcanlog = lsndvec(4)
      prt_isddmcanlog = lsndvec(5)
      deallocate(lsndvec)
c
c-- integer
      n = 7
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/gas_nr,gas_ng,gas_velno,gas_velyes,
     &  prt_npartmax,in_nomp,tsp_nt/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      gas_nr = isndvec(1)
      gas_ng = isndvec(2)
      gas_velno = isndvec(3)
      gas_velyes = isndvec(4)
      prt_npartmax = isndvec(5)
      in_nomp = isndvec(6)
      tsp_nt = isndvec(7)
      deallocate(isndvec)
c
c-- real*8
      n = 1
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/prt_tauddmc/)
      call mpi_bcast(sndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_tauddmc = sndvec(1)
      deallocate(sndvec)
c
c
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       allocate(gas_rarr(gas_nr+1))
       allocate(gas_drarr(gas_nr))
       allocate(gas_curvcent(gas_nr))
       allocate(prt_particles(prt_npartmax))
      endif
c
c-- broadcast data
      call mpi_bcast(gas_rarr,gas_nr+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_drarr,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_curvcent,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c!}}}
      end subroutine bcast_permanent
c
c
c
      subroutine bcast_nonpermanent
c     ------------------------!{{{
      use gasgridmod
      use particlemod
      use timestepmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
*-- scalars:
* real*8 :: tsp_time
* real*8 :: tsp_texp
* real*8 :: tsp_dt
* integer :: prt_nnew
* integer :: prt_nsurf
*
*-- arrays:
* real*8 :: gas_tempkev(gas_nr)
* real*8 :: gas_fcoef(gas_nr)
* real*8 :: gas_sig(gas_nr)
* real*8 :: gas_emitprob(gas_ng,gas_nr)
* real*8 :: gas_ppl(gas_ng,gas_nr)
* real*8 :: gas_ppr(gas_ng,gas_nr)
* real*8 :: gas_opacleakl(gas_ng,gas_nr)
* real*8 :: gas_opacleakr(gas_ng,gas_nr)
* real*8 :: gas_cap(gas_ng,gas_nr)
* real*8 :: gas_wl(gas_ng)
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
c
c-- integer
      n = 2
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/prt_nnew,prt_nsurf/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_nnew = isndvec(1)
      prt_nsurf = isndvec(2)
      deallocate(isndvec)
c
c-- real*8
      n = 3
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/tsp_time,tsp_texp,tsp_dt/)
      call mpi_bcast(sndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      tsp_time = sndvec(1)
      tsp_texp = sndvec(2)
      tsp_dt = sndvec(3)
      deallocate(sndvec)
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0 .and. tsp_it==1) then
       allocate(gas_tempkev(gas_nr))
       allocate(gas_fcoef(gas_nr))
       allocate(gas_sig(gas_nr))
       allocate(gas_emitprob(gas_ng,gas_nr))
       allocate(gas_ppl(gas_ng,gas_nr))
       allocate(gas_ppr(gas_ng,gas_nr))
       allocate(gas_opacleakl(gas_ng,gas_nr))
       allocate(gas_opacleakr(gas_ng,gas_nr))
       allocate(gas_cap(gas_ng,gas_nr))
       allocate(gas_wl(gas_ng))
      endif
c
      call mpi_bcast(gas_tempkev,gas_nr,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_fcoef,gas_nr,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_sig,gas_nr,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emitprob,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_ppl,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_ppr,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacleakl,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacleakr,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_cap,gas_nr*gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_wl,gas_ng,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c!}}}
      end subroutine bcast_nonpermanent
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
