* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
      str_massdd = reshape(str_massdc,[ncell]) !}}}
      if(str_nabund>0) then
       allocate(str_massfrdd(str_nabund,ncell))
       str_massfrdd = reshape(str_massfrdc,[str_nabund,ncell])
      endif
      if(str_ltemp) then
       allocate(str_tempdd(ncell))
       str_tempdd = reshape(str_tempdc,[ncell])
      endif
      if(str_lye) then
       allocate(str_yedd(ncell))
       str_yedd = reshape(str_yedc,[ncell])
      endif
      if(str_lcap) then
       allocate(str_capdd(ncell))
       str_capdd = reshape(str_capdc,[ncell])
      endif
      if(str_ldynfr) then
       allocate(str_dynfrdd(ncell))
       str_dynfrdd = reshape(str_dynfrdc,[ncell])
      endif
      end subroutine scatter_inputstruct
c
c
      subroutine bcast_tbxs(ng)
c     -------------------------
      implicit none
      integer,intent(in) :: ng
************************************************************************
      integer :: dmy
      dmy = ng !use the intent(in) variable
      end subroutine bcast_tbxs
c
c
      subroutine bcast_tbsrc
      end subroutine bcast_tbsrc
c
c
      subroutine allgather_initialrad
c     -------------------------------
      use gridmod
      use gasmod
      grd_evolinit = reshape(gas_eraddens,[grd_ncell])
      end subroutine allgather_initialrad
c
c
      subroutine allgather_gammacap
c     ---------------------------
      use gridmod
      use gasmod
      grd_capgam = reshape(gas_capgam,[grd_ncell])
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
      use inputparmod, only:in_doemiss !TODO: consider making argument
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
      if (in_doemiss) then
        grd_em_cap = reshape(gas_em_cap,[grp_ng,grd_ncell])
        grd_em_capgrey = reshape(gas_em_capgrey,[grd_ncell])
      endif
c
      src_nvacantall(1) = count(prt_isvacant)
      end subroutine bcast_nonpermanent
c
c
      subroutine allgather_leakage
      end subroutine allgather_leakage
c
c
      subroutine reduce_gridtally
************************************************************************
* Reduce the results from the packet transport that are needed for the
* temperature correction.
* - stub
************************************************************************
      use gridmod
      use gasmod
      gas_edep = reshape(grd_tally(1,:),[grd_ncell])
      gas_eraddens = reshape(grd_tally(2,:),[grd_ncell])
      end subroutine reduce_gridtally
c
c
      subroutine reduce_fluxtally
      end subroutine reduce_fluxtally
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
c vim: fdm=marker
