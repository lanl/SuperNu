      module gammaprofmod
c     -------------------
      implicit none
************************************************************************
* Gamma energy deposition profiles exported from Phoenix for the w7 test
* problem.  The profile integrates to 1 in early epochs and to <1 later.
* The profile describes the fraction of the total decay energy that is
* deposited in each of the bins.
************************************************************************
      integer :: gamprf_nt=0
      integer,private :: nr=0
      real*8,allocatable,private :: timevec(:) !(gamprf_nt)  !time in seconds
      real*8,allocatable,private :: profvec(:,:) !(gas_nr,gamprf_nt)  !local deposition fraction
c
      contains
c
      function interp_gamma_profile(t) result(prof)
c     ---------------------------------------------!{{{
      implicit none
      real*8,intent(in) :: t
      real*8 :: prof(nr)
************************************************************************
* interpolate the gamma deposition data for the time requested
************************************************************************
      integer,external :: locate
c
      real*8 :: help
      integer :: it
c
c-- previous time slice in input data
      it = locate(timevec,gamprf_nt,t)
c
c-- linear interpolation
      if(it == gamprf_nt) then
c-- use latest profile in dataset
       prof = profvec(:,it)
      else
       help = timevec(it)
       help = (t - help)/(timevec(it+1) - help)
       prof = (1d0 - help)*profvec(:,it) + help*profvec(:,it+1)
!      write(6,*) 'interp_gam_prof:',timevec(it),timevec(it+1),t,help
      endif
c-- debug output
c!}}}
      end function interp_gamma_profile
c
      subroutine read_gamma_profile(nrin)
c     -----------------------------------!{{{
      implicit none
      integer,intent(in) :: nrin
c
      integer :: istat
c
c-- open input file
      open(4,file='data.gamma_profile',status='old',iostat=istat)
      if(istat/=0) stop 'read_gamprf: file io error'
c-- read header
      read(4,*) !header description
      read(4,*,iostat=istat) gamprf_nt,nr
      if(istat/=0) stop 'read_gamprf: header error'
c-- allocate arrays
      allocate(timevec(gamprf_nt))
      allocate(profvec(nr,gamprf_nt))
c-- read data
      read(4,*,iostat=istat) timevec
      if(istat/=0) stop 'read_gamprf: body error 1'
      read(4,*,iostat=istat) profvec
      if(istat/=0) stop 'read_gamprf: body error 2'
      close(4)
c
c-- verify number of zones
      if(nr/=nrin) stop 'read_gamprf: incompatible dims: nr/=nr'
c!}}}
      end subroutine read_gamma_profile
c
      end module gammaprofmod
