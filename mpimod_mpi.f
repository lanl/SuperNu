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
* integer :: prt_npartmax
* integer :: in_nomp
* integer :: tsp_nt
* integer :: in_seed
*-- real*8
* real*8 :: prt_tauddmc
*-- character
* character(4) :: gas_srctype
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
     &  gas_isshell,prt_isimcanlog,prt_isddmcanlog/)
      call mpi_bcast(lsndvec,n,MPI_LOGICAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      gas_isvelocity = lsndvec(1)
      in_puretran = lsndvec(2)
      gas_isshell = lsndvec(3)
      prt_isimcanlog = lsndvec(4)
      prt_isddmcanlog = lsndvec(5)
      deallocate(lsndvec)
c
c-- integer
      n = 6
      allocate(isndvec(n))
      if(impi==impi0) isndvec = (/gas_nr,gas_ng,
     &  prt_npartmax,in_nomp,tsp_nt,in_seed/)
      call mpi_bcast(isndvec,n,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      gas_nr = isndvec(1)
      gas_ng = isndvec(2)
      prt_npartmax = isndvec(3)
      in_nomp = isndvec(4)
      tsp_nt = isndvec(5)
      in_seed = isndvec(6)
      deallocate(isndvec)
c
c-- real*8
      n = 1
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/prt_tauddmc/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      prt_tauddmc = sndvec(1)
      deallocate(sndvec)
c
c-- character
      call mpi_bcast(gas_srctype,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacanaltype,4,MPI_CHARACTER,
     &  impi0,MPI_COMM_WORLD,ierr)
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
       !prt_done = .false.
      endif
      !if(impi==impi0) then
      !   prt_particles%isvacant=.true.
      !endif
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
*-- real
* real*8 :: tsp_time
* real*8 :: tsp_texp
* real*8 :: tsp_dt
* real*8 :: gas_esurf
* real*8 :: gas_etot
*-- integer
* integer :: prt_nnew
* integer :: prt_nsurf
* integer :: prt_nexsrc
*--
*
*-- arrays:
* real*8 :: gas_temp(gas_nr)
* real*8 :: gas_nvol(gas_nr)
* real*8 :: gas_nvolex(gas_nr)
* real*8 :: gas_emit(gas_nr)
* real*8 :: gas_emitex(gas_nr)
* real*8 :: gas_tempb(gas_nr+1)
* real*8 :: gas_fcoef(gas_nr)
* real*8 :: gas_sig(gas_nr)
* real*8 :: gas_emitprob(gas_ng,gas_nr)
* real*8 :: gas_ppl(gas_ng,gas_nr)
* real*8 :: gas_ppr(gas_ng,gas_nr)
* real*8 :: gas_opacleakl(gas_ng,gas_nr)
* real*8 :: gas_opacleakr(gas_ng,gas_nr)
* real*8 :: gas_cap(gas_ng,gas_nr)
* real*8 :: gas_wl(gas_ng+1)
*
* Variables to be reduced
*-- dim==0
* real*8 :: gas_erad
* real*8 :: gas_eright
* real*8 :: gas_eleft
*-- dim==1
* integer :: gas_numcensus(gas_nr)
* real*8 :: gas_edep(gas_nr)
*-- dim==2
* real*8 :: gas_eraddens(gas_ng,gas_nr)
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:)
      real*8,allocatable :: sndvec(:)
      real*8,allocatable :: sndmat(:,:)

c-- variables to be reduced -----------------------------------
c-- dim==0
      n = 3
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/gas_erad,gas_eright,gas_eleft/)
      call mpi_bcast(sndvec,n,MPI_REAL8,impi0,MPI_COMM_WORLD,
     &  ierr)
c-- copy back
      gas_erad = sndvec(1)
      gas_eright = sndvec(2)
      gas_eleft = sndvec(3)
      deallocate(sndvec)
c-- dim==1,2
      if(impi/=impi0 .and. tsp_it==1) then
         allocate(gas_numcensus(gas_nr))
         allocate(gas_edep(gas_nr))
         allocate(gas_eraddens(gas_ng,gas_nr))
         allocate(gas_luminos(gas_ng))
      endif
!      call mpi_bcast(gas_edep,gas_nr,MPI_REAL8,
!     &  impi0,MPI_COMM_WORLD,ierr)
!      call mpi_bcast(gas_numcensus,gas_nr,MPI_INTEGER,
!     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_eraddens,gas_ng*gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
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
      n = 5
      allocate(sndvec(n))
      if(impi==impi0) sndvec = (/tsp_time,tsp_texp,tsp_dt,gas_esurf,
     & gas_etot/)
      call mpi_bcast(sndvec,n,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      tsp_time = sndvec(1)
      tsp_texp = sndvec(2)
      tsp_dt = sndvec(3)
      gas_esurf = sndvec(4)
      gas_etot = sndvec(5)
      deallocate(sndvec)
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0 .and. tsp_it==1) then
       allocate(gas_temp(gas_nr))
       allocate(gas_nvol(gas_nr))
       allocate(gas_nvolex(gas_nr))
       allocate(gas_emit(gas_nr))
       allocate(gas_emitex(gas_nr))
c
       allocate(gas_tempb(gas_nr+1))
       allocate(gas_fcoef(gas_nr))
       allocate(gas_sig(gas_nr))
       allocate(gas_emitprob(gas_ng,gas_nr))
       allocate(gas_ppl(gas_ng,gas_nr))
       allocate(gas_ppr(gas_ng,gas_nr))
       allocate(gas_opacleakl(gas_ng,gas_nr))
       allocate(gas_opacleakr(gas_ng,gas_nr))
       allocate(gas_cap(gas_ng,gas_nr))
       allocate(gas_wl(gas_ng+1))
c
c-- allocating particle array for helper ranks
       allocate(prt_particles(prt_npartmax))
       prt_particles%isvacant = .true.
c--
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
      call mpi_bcast(gas_temp,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_nvol,gas_nr,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_nvolex,gas_nr,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emit,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emitex,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_bcast(gas_tempb,gas_nr+1,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_fcoef,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_sig,gas_nr,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_emitprob,gas_nr*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_ppl,gas_nr*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_ppr,gas_nr*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacleakl,gas_nr*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_opacleakr,gas_nr*gas_ng,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(gas_cap,gas_nr*gas_ng,MPI_REAL8,
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
      use gasgridmod
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
* integer :: gas_numcensus(gas_nr)
* real*8 :: gas_edep(gas_nr)
*-- dim==2
* real*8 :: gas_eraddens(gas_ng,gas_nr)
************************************************************************
      integer :: n
      integer,allocatable :: isndvec(:),ircvvec(:)
      real*8,allocatable :: sndvec(:),rcvvec(:)
      real*8,allocatable :: sndmat(:,:),rcvmat(:,:)
      real*8 :: help

c
c-- dim==0
      n = 3
      allocate(sndvec(n))
      allocate(rcvvec(n))
      !if(impi==impi0) 
      sndvec = (/gas_erad,gas_eright,gas_eleft/)
      call mpi_reduce(sndvec,rcvvec,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      if(impi==0) then
       gas_erad = rcvvec(1)
       gas_eright = rcvvec(2)
       gas_eleft = rcvvec(3)
      else
       gas_erad = 0d0
       gas_eright = 0d0
       gas_eleft = 0d0
      endif !impi
      deallocate(sndvec)
      deallocate(rcvvec)
c
c-- dim==1
      allocate(isndvec(gas_nr))
      isndvec = gas_numcensus
      call mpi_reduce(isndvec,gas_numcensus,gas_nr,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(isndvec)
c--
      allocate(sndvec(gas_nr))
      sndvec = gas_edep
      call mpi_reduce(sndvec,gas_edep,gas_nr,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(sndvec)
c
      allocate(sndvec(gas_ng))
      sndvec = gas_luminos
      call mpi_reduce(sndvec,gas_luminos,gas_ng,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(sndvec)
c
c-- dim==2
      allocate(sndmat(gas_ng,gas_nr))
      n = gas_ng*gas_nr
      sndmat = gas_eraddens
      call mpi_reduce(sndmat,gas_eraddens,n,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      deallocate(sndmat)
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
      end module mpimod
