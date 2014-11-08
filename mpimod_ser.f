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
* Broadcast the data that does not evolve over time.
* - stub
************************************************************************
      end subroutine bcast_permanent
c
c
      subroutine setup_domain_decomposition
      end subroutine setup_domain_decomposition
c
c
      subroutine scatter_inputstruct(ndim,ncell)
      use inputstrmod
      use gasmod
      implicit none
      integer,intent(in) :: ndim(3)
      integer,intent(in) :: ncell
c
      allocate(str_massdd(ncell))
      if(str_nabund>0) allocate(str_massfrdd(str_nabund,ncell))
      str_massdd = reshape(str_mass,[ncell])
      str_massfrdd = reshape(str_massfr,[str_nabund,ncell])
      end subroutine scatter_inputstruct
c
c
      subroutine bcast_nonpermanent
      use gridmod
************************************************************************
* Broadcast the data that changes with time.
* - stub
************************************************************************
c-- domain decomposition
      grd_temp = reshape(gas_temp,[grd_nx,grd_ny,grd_nz])
      grd_emit = reshape(gas_emit,[grd_nx,grd_ny,grd_nz])
      grd_emitex = reshape(gas_emitex,[grd_nx,grd_ny,grd_nz])
      grd_evolinit = reshape(gas_evolinit,[grd_nx,grd_ny,grd_nz])

      grd_emitprob = reshape(gas_emitprob,[grd_ng,grd_nx,grd_ny,grd_nz])
      grd_cap = reshape(gas_cap,[grd_ng,grd_nx,grd_ny,grd_nz])
      grd_sig = reshape(gas_sig,[grd_nx,grd_ny,grd_nz])
      grd_capgam = reshape(gas_capgam,[grd_nx,grd_ny,grd_nz])
      grd_siggrey = reshape(gas_siggrey,[grd_nx,grd_ny,grd_nz])
      grd_fcoef = reshape(gas_fcoef,[grd_nx,grd_ny,grd_nz])

      end subroutine bcast_nonpermanent
c
c
      subroutine reduce_tally
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      use gridmod
      use gasmod
c
      tot_eextav=tot_eext
      tot_eveloav = tot_evelo
c
      gas_edep = reshape(grd_edep,[grd_nx*grd_ny*grd_nz])
      gas_eraddens = reshape(grd_eraddens,[grd_nx*grd_ny*grd_nz])
      end subroutine reduce_tally
c
c
      subroutine reduce_fluxes
************************************************************************
* Reduce the results from the packet transport.
* - stub
************************************************************************
      end subroutine reduce_fluxes
c
      subroutine reduce_gastemp
c     -------------------------
      use gridmod
      use gasmod
      grd_temp = reshape(gas_temp,[grd_nx,grd_ny,grd_nz])
      end subroutine reduce_gastemp
c
      subroutine scatter_restart_data
************************************************************************
* Distribute data required for restarting a simulation at some time
* - stub
************************************************************************
      end subroutine scatter_restart_data
c
c
      subroutine collect_restart_data
      use particlemod
************************************************************************
* Collect data required for restarting a simulation at some time
*
************************************************************************
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
