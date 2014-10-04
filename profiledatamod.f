      module profiledatamod
c     -------------------
      implicit none
************************************************************************
* Gamma energy deposition profiles exported from Phoenix for the w7 test
* problem.  The profile integrates to 1 in early epochs and to <1 later.
* The profile describes the fraction of the total decay energy that is
* deposited in each of the bins.
************************************************************************
      integer,private :: nr=0
c-- gamma profiles
      integer :: prof_ntgam=0
      real*8,allocatable,private :: timegamvec(:) !(prof_ntgam)  !time in seconds
      real*8,allocatable,private :: profgamvec(:,:) !(gas_nx,prof_ntgam)  !local deposition fraction
c-- trad profiles
      integer :: prof_nttrad=0
      real*8,allocatable,private :: timetradvec(:) !(prof_nttrad)  !time in seconds
      real*8,allocatable,private :: proftradvec(:,:) !(gas_nx,prof_nttrad)  !local deposition fraction
c
      contains
c
c
c
      function gamma_profile(t) result(prof)
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
c-- data sanity test
      if(prof_ntgam==0) stop 'tgam_profile: data not loaded'
c
c-- previous time slice in input data
      it = locate(timegamvec,prof_ntgam,t)
c
c-- linear interpolation
      if(it == prof_ntgam) then
c-- use latest profile in dataset
       prof = profgamvec(:,it)
      else
       help = timegamvec(it)
       help = (t - help)/(timegamvec(it+1) - help)
       prof = (1d0 - help)*profgamvec(:,it) + help*profgamvec(:,it+1)
!      write(6,*) 'interp_gam_prof:',timegamvec(it),timegamvec(it+1),t,help
      endif
c-- debug output
c!}}}
      end function gamma_profile
c
      subroutine read_gamma_profiles(nrin)
c     -----------------------------------!{{{
      implicit none
      integer,intent(in) :: nrin
c
      integer :: istat
c
c-- open input file
      open(4,file='data.gamma_profiles',status='old',iostat=istat)
      if(istat/=0) stop 'read_gamprf: file io error'
c-- read header
      read(4,*) !header description
      read(4,*) !header description
      read(4,*,iostat=istat) prof_ntgam,nr
      if(istat/=0) stop 'read_gamprf: header error'
c-- allocate arrays
      allocate(timegamvec(prof_ntgam))
      allocate(profgamvec(nr,prof_ntgam))
c-- read data
      read(4,*,iostat=istat) timegamvec
      if(istat/=0) stop 'read_gamprf: body error 1'
      read(4,*,iostat=istat) profgamvec
      if(istat/=0) stop 'read_gamprf: body error 2'
      close(4)
c
c-- verify number of zones
      if(nr/=nrin) stop 'read_gamprf: incompatible dims: nr/=nr'
c!}}}
      end subroutine read_gamma_profiles
c
c
c
      function trad_profile(t) result(prof)
c     ---------------------------------------------!{{{
      implicit none
      real*8,intent(in) :: t
      real*8 :: prof(nr)
************************************************************************
* interpolate the trad deposition data for the time requested
************************************************************************
      integer,external :: locate
c
      real*8 :: help
      integer :: it
c
c-- data sanity test
      if(prof_nttrad==0) stop 'trad_profile: data not loaded'
c
c-- previous time slice in input data
      it = locate(timetradvec,prof_nttrad,t)
c
c-- linear interpolation
      if(it == prof_nttrad) then
c-- use latest profile in dataset
       prof = proftradvec(:,it)
      else
       help = timetradvec(it)
       help = (t - help)/(timetradvec(it+1) - help)
       prof = (1d0 - help)*proftradvec(:,it) + help*proftradvec(:,it+1)
       write(6,*) 'interp_trad_prof:',timetradvec(it),timetradvec(it+1),
     &   t,help
      endif
c-- debug output
c!}}}
      end function trad_profile
c
      subroutine read_trad_profiles(nrin)
c     -----------------------------------!{{{
      implicit none
      integer,intent(in) :: nrin
c
      integer :: istat
c
c-- open input file
      open(4,file='data.trad_profiles',status='old',iostat=istat)
      if(istat/=0) stop 'read_tradprf: file io error'
c-- read header
      read(4,*) !header description
      read(4,*) !header description
      read(4,*,iostat=istat) prof_nttrad,nr
      if(istat/=0) stop 'read_tradprf: header error'
c-- allocate arrays
      allocate(timetradvec(prof_nttrad))
      allocate(proftradvec(nr,prof_nttrad))
c-- read data
      read(4,*,iostat=istat) timetradvec
      if(istat/=0) stop 'read_tradprf: body error 1'
      read(4,*,iostat=istat) proftradvec
      if(istat/=0) stop 'read_tradprf: body error 2'
      close(4)
c
c-- verify number of zones
      if(nr/=nrin) stop 'read_tradprf: incompatible dims: nr/=nr'
c!}}}
      end subroutine read_trad_profiles
c
      end module profiledatamod
