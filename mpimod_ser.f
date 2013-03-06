      module mpimod
c     -------------
      implicit none
      integer :: MPI_COMM_WORLD=0
      integer :: MPI_MAX_PROCESSOR_NAME=1
      integer,private :: ierr
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
      subroutine bcast_mutable
c     ------------------------
************************************************************************
* Broadcast the data that changes with time/temperature.
* - stub
************************************************************************
      end subroutine bcast_mutable
c
c
      subroutine reduce_enabs
c     --------------------------
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      end subroutine reduce_enabs
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
      subroutine mpi_init(ierr)
      implicit none
      integer :: ierr
      end subroutine mpi_init
c
      subroutine mpi_comm_rank(MPI_COMM_WORLD,impi,ierr)
      implicit none
      integer :: MPI_COMM_WORLD
      integer :: impi,ierr
      end subroutine mpi_comm_rank
c
      subroutine mpi_comm_size(MPI_COMM_WORLD,nmpi,ierr)
      implicit none
      integer :: MPI_COMM_WORLD
      integer :: nmpi,ierr
      end subroutine mpi_comm_size
c
      subroutine mpi_get_processor_name(pname,ilen,ierr)
      implicit none
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen,ierr
      end subroutine mpi_get_processor_name
c
      subroutine mpi_barrier(MPI_COMM_WORLD,ierr)
      implicit none
      integer :: MPI_COMM_WORLD,ierr
      end subroutine mpi_barrier
c
      subroutine mpi_finalize(ierr)
      implicit none
      integer :: ierr
      end subroutine mpi_finalize
c
      end module mpimod
