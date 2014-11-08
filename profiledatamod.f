      module profiledatamod
c     -------------------
      implicit none
************************************************************************
* Gamma energy deposition profiles exported from Phoenix for the w7 test
* problem.  The profile integrates to 1 in early epochs and to <1 later.
* The profile describes the fraction of the total decay energy that is
* deposited in each of the bins.
************************************************************************
      integer :: prof_nr=0
c-- gamma profiles
      integer :: prof_ntgam=0
      real*8,allocatable :: prof_timegamvec(:) !(prof_ntgam)  !time in seconds
      real*8,allocatable :: prof_profgamvec(:,:) !(grd_nx,prof_ntgam)  !local deposition fraction
c
      contains
c
c
c
      function gamma_profile(t) result(prof)
c     ---------------------------------------------!{{{
      implicit none
      real*8,intent(in) :: t
      real*8 :: prof(prof_nr)
************************************************************************
* interpolate the gamma deposition data for the time requested
************************************************************************
      integer,external :: locate
c
      real*8 :: help
      integer :: it
c
c-- data sanity test
      if(prof_ntgam==0) stop 'tgam_profile: data not loaded'
c
c-- previous time slice in input data
      it = locate(prof_timegamvec,prof_ntgam,t)
c
c-- linear interpolation
      if(it == prof_ntgam) then
c-- use latest profile in dataset
       prof = prof_profgamvec(:,it)
      else
       help = prof_timegamvec(it)
       help = (t - help)/(prof_timegamvec(it+1) - help)
       prof = (1d0 - help)*prof_profgamvec(:,it) +
     &   help*prof_profgamvec(:,it+1)
!      write(6,*) 'interp_gam_prof:',prof_timegamvec(it),prof_timegamvec(it+1),t,help
      endif
c-- debug output
c!}}}
      end function gamma_profile
c
      subroutine read_gamma_profiles(ndim)
c     -----------------------------------!{{{
      implicit none
      integer,intent(in) :: ndim(3)
c
      integer :: istat
c
c-- open input file
      open(4,file='data.gamma_profiles',status='old',iostat=istat)
      if(istat/=0) stop 'read_gamprf: file io error'
c-- read header
      read(4,*) !header description
      read(4,*) !header description
      read(4,*,iostat=istat) prof_ntgam,prof_nr
      if(istat/=0) stop 'read_gamprf: header error'
c-- allocate arrays
      allocate(prof_timegamvec(prof_ntgam))
      allocate(prof_profgamvec(prof_nr,prof_ntgam))
c-- read data
      read(4,*,iostat=istat) prof_timegamvec
      if(istat/=0) stop 'read_gamprf: body error 1'
      read(4,*,iostat=istat) prof_profgamvec
      if(istat/=0) stop 'read_gamprf: body error 2'
      close(4)
c
c-- verify number of zones
      if(prof_nr/=ndim(1)) stop 'read_gamprf: incompatible dims: nr/=nr'
c!}}}
      end subroutine read_gamma_profiles
c
      end module profiledatamod
