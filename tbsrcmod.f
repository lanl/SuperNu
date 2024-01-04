* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
      module tbsrcmod
c     ---------------
************************************************************************
* Heating rates from runtime files
************************************************************************
      character(13),private :: fname_dyn='hrate_dyn.dat'
      character(13),private :: fname_wnd='hrate_wnd.dat'
      integer, parameter :: tbs_ncol = 8
      integer :: tbs_ntdyn = 0
      integer :: tbs_ntwnd = 0
c-- thermalization formula options (first 3 rate columns time+totals)
      integer :: tbs_dynopt(tbs_ncol-3)
      integer :: tbs_wndopt(tbs_ncol-3)
c-- thermalization parameters (use is based on formula options)
      real*8 :: tbs_dynpars(tbs_ncol-3)
      real*8 :: tbs_wndpars(tbs_ncol-3)
c-- coarsened heating data
      real*8,allocatable :: tbs_dynhr(:,:) !(ncol,nt)
      real*8,allocatable :: tbs_wndhr(:,:) !(ncol,nt)
c
      save
      public
c
      contains
c
c
c
      subroutine tbsrc_dealloc
c     ------------------------
      implicit none
      if(allocated(tbs_dynhr)) deallocate(tbs_dynhr)
      if(allocated(tbs_wndhr)) deallocate(tbs_wndhr)
      end subroutine tbsrc_dealloc
c
c
c
      subroutine read_tbsrc
c     ---------------------
      implicit none
************************************************************************
* Read heating rate tables to heating rate arrays.
************************************************************************
      logical :: lexist_dyn, lexist_wnd
      integer :: ierr
      character(2) :: dmy
c
c-- open file for dynamical ejecta rates
      inquire(file='Source/'//fname_dyn,exist=lexist_dyn)
      if(lexist_dyn) then
        open(4,file='Source/'//fname_dyn,status='old',iostat=ierr)
c-- read past first headers
        read(4,*)
        read(4,*)
        read(4,*)
        read(4,*)
        read(4,*)
c-- read thermalization parameters
        read(4,*,iostat=ierr) tbs_dynopt
        read(4,*,iostat=ierr) tbs_dynpars
        read(4,*)
c-- sanity check
        if(any(tbs_dynopt>1)) stop 'read_tbsrc: dynopt>1'
        if(any(tbs_dynopt<0)) stop 'read_tbsrc: dynopt<0'
c-- read number of time steps
        read(4,*,iostat=ierr) dmy, tbs_ntdyn
c-- sanity check
        if(tbs_ntdyn <= 0) stop 'read_tbsrc: ntdyn <= 0'
c-- allocate dynamical ejecta source table
        allocate(tbs_dynhr(tbs_ncol,tbs_ntdyn))
c-- read past headers
        read(4,*)
        read(4,*)
c-- read table
        read(4,*,iostat=ierr) tbs_dynhr
        close(4)
      endif
c
c-- open file for wind rates
      inquire(file='Source/'//fname_wnd,exist=lexist_wnd)
      if(lexist_wnd) then
        open(4,file='Source/'//fname_wnd,status='old',iostat=ierr)
c-- read past first headers
        read(4,*)
        read(4,*)
        read(4,*)
        read(4,*)
        read(4,*)
c-- read thermalization parameters
        read(4,*,iostat=ierr) tbs_wndopt
        read(4,*,iostat=ierr) tbs_wndpars
        read(4,*)
c-- sanity check
        if(any(tbs_wndopt>1)) stop 'read_tbsrc: wndopt>1'
        if(any(tbs_wndopt<0)) stop 'read_tbsrc: wndopt<0'
c-- read number of time steps
        read(4,*,iostat=ierr) dmy, tbs_ntwnd
c-- sanity check
        if(tbs_ntwnd <= 0) stop 'read_tbsrc: ntwnd <= 0'
c-- allocate dynamical ejecta source table
        allocate(tbs_wndhr(tbs_ncol,tbs_ntwnd))
c-- read past headers
        read(4,*)
        read(4,*)
c-- read table
        read(4,*,iostat=ierr) tbs_wndhr
        close(4)
      endif
c
c-- sanity check
      if(.not.lexist_dyn .and. .not.lexist_wnd) stop
     &   'read_tbsrc: no source files'
c
      end subroutine read_tbsrc
c
      end module tbsrcmod
