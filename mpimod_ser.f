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
      use gasgridmod
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
      use gasgridmod
************************************************************************
* Broadcast the data that changes with time.
* - stub
************************************************************************
c-- domain decomposition
      grd_temp = reshape(dd_temp,[grd_nx,grd_ny,grd_nz])
      grd_emit = reshape(dd_emit,[grd_nx,grd_ny,grd_nz])
      grd_emitex = reshape(dd_emitex,[grd_nx,grd_ny,grd_nz])
      grd_evolinit = reshape(dd_evolinit,[grd_nx,grd_ny,grd_nz])

      grd_emitprob = reshape(dd_emitprob,[gas_ng,grd_nx,grd_ny,grd_nz])
      grd_cap = reshape(dd_cap,[gas_ng,grd_nx,grd_ny,grd_nz])
      grd_sig = reshape(dd_sig,[grd_nx,grd_ny,grd_nz])
      grd_capgam = reshape(dd_capgam,[grd_nx,grd_ny,grd_nz])
      grd_siggrey = reshape(dd_siggrey,[grd_nx,grd_ny,grd_nz])
      grd_fcoef = reshape(dd_fcoef,[grd_nx,grd_ny,grd_nz])

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
      use gasgridmod
c
      gas_eextav=gas_eext
      gas_eveloav = gas_evelo
c
      dd_edep = reshape(grd_edep,[grd_nx*grd_ny*grd_nz])
      dd_eraddens = reshape(grd_eraddens,[grd_nx*grd_ny*grd_nz])
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
      use gasgridmod
      grd_temp = reshape(dd_temp,[grd_nx,grd_ny,grd_nz])
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
