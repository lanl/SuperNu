      module inputparmod
c     ------------------
      use physconstmod
      implicit none
************************************************************************
* input parameters
************************************************************************
c-- write stdout to file
      logical :: in_grab_stdout = .false. !write stdout to file
c-- gas grid
      integer :: in_nr = 0 !# spatial grid in spherical geom
      integer :: in_ng = 0 !# groups
      logical :: in_isvelocity = .true.  !switch underlying grid between spatial+static to velocity+expanding
      logical :: in_isshell = .false.
      real*8 :: in_l0 = 0d0  !innermost radius of the domain
      real*8 :: in_lr = 0d0  !spatial length of the domain
      real*8 :: in_velout = 0d0  !cm/s
      real*8 :: in_totmass = 0d0  !g
      real*8 :: in_templ0 = 0d0 !inner bound temperature in keV
      real*8 :: in_sigcoef = 0d0 !power law opacity coefficient
      real*8 :: in_sigtpwr = 0d0 !power law opacity temperature exponent
      real*8 :: in_sigrpwr = 0d0 !power law opacity density exponent
c-- flat-structure parameters
      real*8 :: in_consttempkev = 0d0  !keV
      logical :: in_solidni56 = .false.  !pure nickel56 atmosphere
c
c-- random number generator
      integer :: in_seed = 1984117 !starting point of random number generator
c
c-- particles
      integer :: in_ns = 0   !# source particles generated per time step
      integer :: in_npartmax = 0 !total # particles allowed
      logical :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      real*8 :: in_alpha = 1d0 !time centering control parameter [0,1]
c
c-- time step
      real*8 :: in_tfirst = 0d0 !first point in time evolution
      real*8 :: in_tlast = 0d0  !last point in time evolution
      integer :: in_nt = -1   !# time steps
c
c-- parallelization
      integer :: in_nomp = 1       !# openmp threads
c
c-- group structure
      real*8 :: in_wlmin = 1000d0     !lower wavelength boundary in output spectrum
      real*8 :: in_wlmax = 30000d0    !upper wavelength boundary in output spectrum
c
c-- opacity (cm^2/gram)
      real*8 :: in_opcapgam = .06d0   ![cm^2/g] extinction coefficient for gamma radiation
      real*8 :: in_opcap = 0d0        !additional gray opacity (for testing with in_nobbopac only!)
      real*8 :: in_epsline = 1d0      !line absorption fraction (the rest is scattering)
      logical :: in_nobbopac = .false.    !turn off bound-bound opacity
      logical :: in_nobfopac = .false.    !turn off bound-bound opacity
      logical :: in_noffopac = .false.    !turn off bound-bound opacity
c
c-- misc
      character(4) :: in_opacdump = 'off'    !off|one|each|all: write opacity data to file
      character(4) :: in_pdensdump = 'off'   !off|one|each: write partial densities to file
c
c-- runtime parameter namelist
      namelist /inputpars/
     & in_nr,in_ng,in_isvelocity,in_isshell,in_lr,in_l0,
     & in_totmass,in_templ0,in_velout,in_consttempkev,in_solidni56,
     & in_sigcoef,in_sigtpwr,in_sigrpwr,
     & in_seed,in_ns,in_npartmax,in_puretran,in_alpha,
     & in_tfirst,in_tlast,in_nt,
     & in_grab_stdout,in_nomp,
     & in_wlmin,in_wlmax,
     & in_opcapgam,in_opcap,in_epsline,in_nobbopac,in_nobfopac,
     & in_noffopac,
     & in_opacdump,in_pdensdump
c
      public
      private inputpars
      save
c
      contains
c
      subroutine read_inputpars
c     -------------------------
      implicit none
************************************************************************
* read the input parameter namelist
************************************************************************
      character(15),parameter :: fname='input.par'
c
c-- read namelist
      open(4,file=fname,status='old',err=66)
      read(4,nml=inputpars,end=67,err=68)
      close(4)
c
      return
66    stop 'read_inputpars: namelist input file missing: input.par'
67    stop 'read_inputpars: namelist missing or bad in input.par'
68    stop 'read_inputpars: ivalid parameters or values in namelist'
      end subroutine read_inputpars
c
c
c
      subroutine parse_inputpars(nmpi)
c     --------------------------------
      use miscmod, only:warn
c$    use omp_lib
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* parse the input parameter namelist
************************************************************************
      integer :: istat
c$    integer :: i
c
c-- redirect stdout to file if selected
      if(in_grab_stdout) then!{{{
       write(6,*) 'write stdout to fort.6'
       open(6,file='fort.6',action='write',status='replace',recl=2000,
     &   iostat=istat) !write stdout to file
       if(istat/=0) stop 'parse_inputpars: open fort.6 error'
       call banner
      endif!}}}
c
c-- dump namelist to stdout
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'!{{{
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(in_nr<=0) stop 'in_nr invalid'
      if(in_ng<=0) stop 'in_ng invalid'
c
      if(in_isvelocity) then
       if(in_lr>0d0) stop 'vel grid: use in_velout, not in_lr'
       if(in_velout<=0d0) stop 'vel grid: use in_velout, not in_lr'
      else
       if(in_lr<=0d0) stop 'static grid: use in_lr, not in_velout'
       if(in_l0<0d0) stop 'static grid: in_l0 invalid'
       if(in_velout>0d0) stop 'static grid: use in_lr, not in_velout'
      endif
c
      if(in_ns<=0) stop 'in_ns invalid'
      if(in_npartmax<=0) stop 'in_npartmax invalid'
      if(in_alpha>1d0 .or. in_alpha<0d0) stop 'in_alpha invalid'
c
      if(in_totmass<=0d0) stop 'in_totmass <= 0'
      if(in_consttempkev<=0d0) stop 'in_consttempkev <= 0'
c
      if(in_nt<1) stop 'in_nt invalid'
      if(in_tfirst<=0d0) stop 'in_tfirst invalid'
      if(in_tlast<in_tfirst) stop 'in_tlast invalid'
c
      if(in_nr<=0) stop 'in_nr invalid'
c
      if(in_wlmin<0d0) stop 'in_wlmin invalid'
      if(in_wlmax<=0d0 .or. in_wlmax<in_wlmin) stop 'in_wlmax invalid'
c
      if(in_opcapgam<0d0) stop 'in_opcapgam invalid'
      if(in_opcap<0d0) then
       stop 'in_opcap invalid'
      elseif(in_opcap>0d0) then
       call warn('read_inputpars',
     &   'gray opacity added! For testing uses only!')
      endif
      if(in_epsline<0d0 .or. in_epsline>1d0) stop 'in_epsline invalid'
c
      if(in_nobbopac) call warn('read_inputpars','bb opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','ff opacity disabled!')
c
      if(trim(in_opacdump)=='off') then
      elseif(trim(in_opacdump)=='one') then
      elseif(trim(in_opacdump)=='each') then
      elseif(trim(in_opacdump)=='all') then
       call warn('read_inputpars',
     &   "in_opacdump=='all' will generate a big data file!")
      else
       stop 'in_opacdump invalid'
      endif
c
      if(trim(in_pdensdump)=='off') then
      elseif(trim(in_pdensdump)=='one') then
      elseif(trim(in_pdensdump)=='each') then
      else
       stop 'in_pdensdump invalid'
      endif
c
c-- set the number of threads
c$    if(.false.) then!{{{
c-- serial run
       in_nomp = 1
c$    else
c-- openmp run
c$     if(in_nomp/=0) call omp_set_num_threads(in_nomp)
c$omp parallel shared(in_nomp) private(i)
c$     i = omp_get_num_threads()
c$     if(in_nomp/=0 .and. i/=in_nomp)
c$   &   stop 'read_inputpars: in_nomp error'
c$     in_nomp = i
c$omp end parallel
c$    endif
      write(6,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &  nmpi,in_nomp,nmpi*in_nomp
      if(in_grab_stdout) then
       write(0,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &   nmpi,in_nomp,nmpi*in_nomp
      endif!}}}!}}}
c
      end subroutine parse_inputpars
c
      end module inputparmod
