      module mpimod
c     -------------
      implicit none
      integer :: MPI_COMM_WORLD=0
      integer :: MPI_MAX_PROCESSOR_NAME=13
      integer,private :: ierr=0
      integer :: impi=0  !mpi rank
      integer :: impi0=0 !master mpi rank
      integer :: nmpi=1  !number of mpi tasks
c
      contains
c
      subroutine bcast_permanent
c     --------------------------
      use inputparmod, only:in_nomp
      implicit none
************************************************************************
* Broadcast the data that does not evolve over time (or temperature).
* - stub
************************************************************************
      end subroutine bcast_permanent
c
c
      subroutine bcast_nonpermanent
c     ------------------------
************************************************************************
* Broadcast the data that changes with time/temperature.
* - stub
************************************************************************
      end subroutine bcast_nonpermanent
c
c
      subroutine reduce_tally
c     --------------------------
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      end subroutine reduce_tally
c
c
      subroutine reduce_fluxes
c     -------------------------
************************************************************************
* Reduce the results from the packet transport.
* - stub
************************************************************************
      end subroutine reduce_fluxes
c
c
      subroutine collect_restart_data
************************************************************************
* Collect data required for restarting a simulation at some time
*
************************************************************************
      use particlemod
      prt_tlyrandarr=prt_tlyrand
      end subroutine collect_restart_data
c
c
      subroutine mpi_init(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_init
c
      subroutine mpi_comm_rank(mpi_comm,impi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: impi_,ierr_
      ierr_ = ierr
      impi_ = impi
      mpi_comm = MPI_COMM_WORLD
      impi = impi
      end subroutine mpi_comm_rank
c
      subroutine mpi_comm_size(mpi_comm,nmpi_,ierr_)
      implicit none
      integer :: mpi_comm
      integer :: nmpi_,ierr_
      ierr_ = ierr
      nmpi_ = nmpi
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_comm_size
c
      subroutine mpi_get_processor_name(pname,ilen_,ierr_)
      implicit none
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen_,ierr_
      pname = 'NOT AVAILABLE'
      ierr_ = ierr
      ilen_ = 1
      end subroutine mpi_get_processor_name
c
      subroutine mpi_barrier(mpi_comm,ierr_)
      implicit none
      integer :: mpi_comm,ierr_
      ierr_ = ierr
      mpi_comm = MPI_COMM_WORLD
      end subroutine mpi_barrier
c
      subroutine mpi_finalize(ierr_)
      implicit none
      integer :: ierr_
      ierr_ = ierr
      end subroutine mpi_finalize
c
      end module mpimod
