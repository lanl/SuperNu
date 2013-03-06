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
      use timestepmod, only:tim_ntim
      use rgridmod
      use rngenmod, only:rng_offset_all,rng_offset
      use ggridmod
      use fluxmod
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* Also once the constants are broadcasted, all allocatable arrays are
* allocated.
*
* data:
* - rgrid
* - ggrid_wl
* constants:
* - in_ncostf,in_nphif,in_nwlf,in_wlmin,in_wlmax
* - in_nwlg,in_epsline,in_niwlem
* - in_nr,rg_ncr,gg_ncg,gg_npacket
* - in_nomp,tim_ntim,in_ntc
* - gg_xi2beta,gg_dxwin
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
      if(impi==impi0) isndvec = (/in_ncostf,in_nphif,in_nwlf,in_nwlg,
     &  in_niwlem,in_nr,rg_ncr,gg_ncg,gg_npacket,in_nomp,tim_ntim,
     &  in_ntc/)
      call mpi_bcast(isndvec,12,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_ncostf = isndvec(1)
      in_nphif = isndvec(2)
      in_nwlf = isndvec(3)
      in_nwlg = isndvec(4)
      in_niwlem = isndvec(5)
      in_nr = isndvec(6)
      rg_ncr = isndvec(7)
      gg_ncg = isndvec(8)
      gg_npacket = isndvec(9)
      in_nomp = isndvec(10)
      tim_ntim = isndvec(11)
      in_ntc = isndvec(12)
      deallocate(isndvec)
c
c-- broadcast real*8 constants
      allocate(sndvec(7))
      if(impi==impi0) then
       sndvec(1:3) = (/in_wlmin,in_wlmax,in_epsline/)
       sndvec(4:7) = (/gg_xi2beta,gg_dxwin,flx_wlhelp,flx_wlminlg/)
      endif
      call mpi_bcast(sndvec,7,MPI_REAL8,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- copy back
      in_wlmin = sndvec(1)
      in_wlmax = sndvec(2)
      in_epsline = sndvec(3)
      gg_xi2beta = sndvec(4)
      gg_dxwin = sndvec(5)
      flx_wlhelp = sndvec(6)
      flx_wlminlg = sndvec(7)
      deallocate(sndvec)
c
c$    if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c
c-- allocate all arrays. These are deallocated in dealloc_all.f
      if(impi/=impi0) then
       allocate(rgrid(rg_ncr))
       allocate(ggrid(gg_ncg))
       allocate(ggrid_wl(in_nwlg))
       allocate(ggrid_cap(gg_ncg,in_nwlg))
       allocate(ggrid_icapbb(0:in_niwlem,gg_ncg))
       allocate(flx_flux(in_nwlf,in_ncostf,in_nphif))
       allocate(flx_iflux(in_nwlf,in_ncostf,in_nphif))
      endif
      allocate(rng_offset(0:in_nomp-1))
c
c-- broadcast data
      call mpi_bcast(ggrid_wl,in_nwlg,MPI_REAL8,
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
      use ggridmod
      use timestepmod
      implicit none
************************************************************************
* Broadcast the data that changes with time/temperature.
* constants:
* - tim_itc
* data:
* - ggrid
* - ggrid_cap
* - ggrid_icapbb
************************************************************************
      call mpi_bcast(tim_itc,1,MPI_INTEGER,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      call mpi_bcast(ggrid_cap,in_nwlg*gg_ncg,MPI_REAL,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c WARNING, this may not work on heterogeneous clusters!!!
      call mpi_bcast(ggrid_icapbb,sizeof(ggrid_icapbb),MPI_BYTE,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c WARNING, this may not work on heterogeneous clusters!!!
      call mpi_bcast(ggrid,sizeof(ggrid),MPI_BYTE,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine bcast_mutable
c
c
c
      subroutine reduce_enabs
c     -----------------------
      use ggridmod
      use timestepmod
      use timingmod
      implicit none
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - ggrid%enabs_e
* - ggrid%enabs_c
* - t_pckt_stat !min,mean,max
************************************************************************
      real*8 :: sndvec(gg_ncg),rcvvec(gg_ncg)
      real*8 :: help
c
      sndvec = ggrid(:)%enabs_e
      call mpi_reduce(sndvec,rcvvec,gg_ncg,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      ggrid(:)%enabs_e = rcvvec
c
      sndvec = ggrid(:)%enabs_c
      call mpi_reduce(sndvec,rcvvec,gg_ncg,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
      ggrid(:)%enabs_c = rcvvec
c
      if(tim_itc>0) then
       sndvec = ggrid%enocoll
       call mpi_reduce(sndvec,rcvvec,gg_ncg,MPI_REAL8,MPI_SUM,
     &   impi0,MPI_COMM_WORLD,ierr)
       ggrid%enocoll = rcvvec
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
c
      subroutine reduce_fluxes
c     ------------------------
      use inputparmod
      use fluxmod
      use ggridmod
      implicit none
************************************************************************
* Reduce the results from the last packet transport.
* - flx_flux
* - flx_iflux
* - ggrid%enostor
*
* - flx_npckt
* - flx_dxtotmax
* - flx_nabsmax
* - flx_nscatmax
************************************************************************
      integer :: i
      real*8 :: help
      integer :: isndmat(in_nwlf,in_ncostf,in_nphif)
      real*8 :: sndmat(in_nwlf,in_ncostf,in_nphif)
c
c-- reduce fluxes
c-- flx_flux
      sndmat = flx_flux
      i = in_nwlf*in_ncostf*in_nphif
      call mpi_reduce(sndmat,flx_flux,i,MPI_REAL8,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- flx_iflux
      isndmat = flx_iflux
      i = in_nwlf*in_ncostf*in_nphif
      call mpi_reduce(isndmat,flx_iflux,i,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- reduce total variables
      i = flx_npckt
      call mpi_reduce(i,flx_npckt,1,MPI_INTEGER,MPI_SUM,
     &  impi0,MPI_COMM_WORLD,ierr)
c
c-- reduce real*8 max variables
      help = flx_dxtotmax
      call mpi_reduce(help,flx_dxtotmax,1,MPI_REAL8,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
c-- reduce integer max variables
      i = flx_nabsmax
      call mpi_reduce(i,flx_nabsmax,1,MPI_INTEGER,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
      i = flx_nscatmax
      call mpi_reduce(i,flx_nscatmax,1,MPI_INTEGER,MPI_MAX,
     &  impi0,MPI_COMM_WORLD,ierr)
c
      end subroutine reduce_fluxes
c
      end module mpimod
