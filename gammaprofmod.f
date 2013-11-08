      module gammaprofmod
c     -------------------
      implicit none
************************************************************************
* Gamma energy deposition profiles exported from Phoenix for the w7 test
* problem.  The profile integrates to 1 in early epochs and to <1 later.
* The profile describes the fraction of the total decay energy that is
* deposited in each of the bins.
************************************************************************
      integer :: gamprf_nt = 0
      real*8,allocatable :: gamprf_time(:) !(gamprf_nt)  !time in seconds
      real*8,allocatable :: gamprf_prof(:,:) !(gas_nr,gamprf_nt)  !local deposition fraction
c
      contains
c
      subroutine read_gamma_profile(nrin)
c     ---------------------------------
      implicit none
      integer,intent(in) :: nrin
c
      integer :: nr,istat
c
c-- open input file
      open(4,file='data.gamma_profile',status='old',iostat=istat)
      if(istat/=0) stop 'read_gamprf: file io error'
c-- read header
      read(4,*) !header description
      read(4,*,iostat=istat) gamprf_nt,nr
      if(istat/=0) stop 'read_gamprf: header error'
c-- allocate arrays
      allocate(gamprf_time(gamprf_nt))
      allocate(gamprf_prof(nr,gamprf_nt))
c-- read data
      read(4,*,iostat=istat) gamprf_time
      if(istat/=0) stop 'read_gamprf: body error 1'
      read(4,*,iostat=istat) gamprf_prof
      if(istat/=0) stop 'read_gamprf: body error 2'
      close(4)
c
c-- verify number of zones
      if(nr/=nrin) stop 'read_gamprf: incompatible dimensions: nr/=nr'
c
      end subroutine read_gamma_profile
c
      end module gammaprofmod
