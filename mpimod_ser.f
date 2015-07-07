      module mpimod
c     -------------
      implicit none
      integer :: MPI_COMM_WORLD=0
      integer,parameter :: MPI_MAX_PROCESSOR_NAME=13
      integer,private :: ierr=0
      integer :: impi=0  !mpi rank
      integer,parameter :: impi0=0 !master mpi rank
      logical :: lmpi0   !true on the master mpi rank
      integer :: nmpi=1  !number of mpi tasks
c
      save
c
      contains
c
      subroutine bcast_permanent
      end subroutine bcast_permanent
c
c
      subroutine scatter_inputstruct(ndim,icell1,ncell)
      use inputstrmod!{{{
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(out) :: icell1,ncell
************************************************************************
* mpi_scatter the input structure to all ranks in the worker comm.
************************************************************************
      integer :: dmy
c
      icell1 = 1
      ncell = str_nc
c
      dmy = ndim(1) !use the intent(in) variable
      allocate(str_massdd(ncell))
      if(str_nabund>0) then
       allocate(str_massfrdd(str_nabund,ncell))
       str_massfrdd = reshape(str_massfrdc,[str_nabund,ncell])
      endif
      str_massdd = reshape(str_massdc,[ncell]) !}}}
      end subroutine scatter_inputstruct
c
c
      subroutine allgather_gammacap
c     ---------------------------
      use gridmod
      use gasmod
      grd_capgrey = reshape(gas_capgam,[grd_ncell])
      grd_emitex = reshape(gas_emitex,[grd_ncell])
      end subroutine allgather_gammacap
c
c
      subroutine allreduce_gammaenergy
      end subroutine allreduce_gammaenergy
c
c
      subroutine bcast_nonpermanent
      use gridmod
      use gasmod
      use groupmod
      use sourcemod
      use particlemod
************************************************************************
* Broadcast the data that changes with time.
* - stub
************************************************************************
c-- domain decomposition
      grd_tempinv = reshape(1d0/gas_temp,[grd_ncell])
      grd_emit = reshape(gas_emit,[grd_ncell])
      grd_emitex = reshape(gas_emitex,[grd_ncell])
      grd_evolinit = reshape(gas_evolinit,[grd_ncell])
c
      grd_cap = reshape(gas_cap,[grp_ng,grd_ncell])
      grd_sig = reshape(gas_sig,[grd_ncell])
      grd_capgrey = reshape(gas_capgrey,[grd_ncell])
      grd_fcoef = reshape(gas_fcoef,[grd_ncell])
c
      src_nvacantall(1) = count(prt_isvacant)
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
      end subroutine allgather_leakage
c
c
      subroutine reduce_tally
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      use totalsmod
      use gridmod
      use gasmod
      gas_edep = reshape(grd_tally(1,:),[grd_ncell])
      gas_eraddens = reshape(grd_tally(2,:),[grd_ncell])
      end subroutine reduce_tally
c
c
      subroutine reduce_fluxes
      end subroutine reduce_fluxes
c
      subroutine reduce_gastemp
c     -------------------------
      use gridmod
      use gasmod
      grd_tempinv = reshape(1d0/gas_temp,[grd_ncell])
      end subroutine reduce_gastemp
c
      subroutine scatter_restart_data
      end subroutine scatter_restart_data
c
c
      subroutine collect_restart_data
      end subroutine collect_restart_data
c
c
      subroutine mpimod_dealloc
      end subroutine mpimod_dealloc
c
c-- MPI intrinsics
c-----------------
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
