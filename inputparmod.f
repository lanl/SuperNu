      module inputparmod
c     ------------------
      USE kindmod
      implicit none
************************************************************************
* input parameters
************************************************************************
c-- gas grid
      INTEGER(iknd) :: in_nr = 0 !# spatial grid in spherical geom
      INTEGER(iknd) :: in_ng = 0 !# groups
      LOGICAL :: in_isvelocity = .true.  !switch underlying grid between spatial+static to velocity+expanding
      REAL(rknd) :: in_lr = 0d0  !spatial length of the domain
c
c-- particles
      INTEGER(iknd) :: in_seed = 0 !starting point of random number generator
      INTEGER(iknd) :: in_ns = 0   !# source particles generated per time step
      INTEGER(iknd) :: in_npartmax = 0 !total # particles allowed
      LOGICAL :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      REAL(rknd) :: in_alpha = 1d0 !time centering control parameter [0,1]
c
c-- time step
      REAL(rknd) :: in_tfirst = 0d0 !first point in time evolution
      REAL(rknd) :: in_tlast = 0d0  !last point in time evolution
      INTEGER(iknd) :: in_nt = -1   !# time steps

c
c-- parallelization
      LOGICAL :: in_grab_stdout = .false. !write stdout to file
      INTEGER(ikind) :: in_nomp = 1       !# openmp threads
c
c-- group structure
      REAL(rknd) :: in_wlmin = 1000d0     !lower wavelength boundary in output spectrum
      REAL(rknd) :: in_wlmax = 30000d0    !upper wavelength boundary in output spectrum
c
c-- opacity (cm^2/gram)
      REAL(rknd) :: in_opcapgam = .06d0    ![cm^2/g] extinction coefficient for gamma radiation
      REAL(rknd) :: in_opcap = .0d0       !additional gray opacity (for testing with in_nobbopac only!)
      REAL(rknd) :: in_epsline = 1d0      !line absorption fraction (the rest is scattering)
      LOGICAL :: in_nobbopac = .false.    !turn off bound-bound opacity
      LOGICAL :: in_nobfopac = .false.    !turn off bound-bound opacity
      LOGICAL :: in_noffopac = .false.    !turn off bound-bound opacity
c
c-- misc
      CHARACTER(4) :: in_opacdump = 'off'    !off|one|each|all: write opacity data to file
      CHARACTER(4) :: in_pdensdump = 'off'   !off|one|each: write partial densities to file
c
c-- runtime parameter namelist
      namelist /inputpars/
     & in_nr,in_ng,in_isvelocity,in_lr,
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
      character(9),parameter :: fname='input.par'
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
      real :: rhelp
c$    integer :: i
c
c-- redirect stdout to file if chosen so
      if(in_grab_stdout) then!{{{
       write(6,*) 'write stdout to fort.6'
       open(6,file='fort.6',action='write',status='replace',recl=2000,
     &   iostat=istat) !write stdout to file
       if(istat/=0) stop 'parse_inputpars: open fort.6 error'
       call banner
      endif!}}}
c
c-- check read values
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)

      if(in_nr<=0) stop 'invalid'
      if(in_ng<=0) stop 'invalid'
      if(in_lr<=0d0) stop 'invalid'
      
      if(in_ns<=0) stop 'in_ns invalid'
      if(in_npartmax<=0) stop 'in_npartmax invalid'
      if(in_alpha>1d0 .or. in_alpha<0d0) stop 'in_alpha invalid'
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'!{{{
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(in_nt<1) stop 'in_nt invalid'
      if(in_tfirst<=0d0) stop 'in_tfirst invalid'
      if(in_tlast<in_tfirst) stop 'in_tlast invalid'
      if(trim(in_enoresrv_init)=='locdep') then
      elseif(trim(in_enoresrv_init)=='restart') then
      elseif(trim(in_enoresrv_init)=='off') then
       call warn('read_inputpars','decay energy before tfirst dropped!')
      else
       stop 'in_enoresrv_init invalid'
      endif
c
      if(in_ntc<0) stop 'in_ntc invalid'
c
      if(in_nr<=0) stop 'in_nr invalid'
      if(in_l2packet<0) stop 'in_l2packet < 0'
      if(in_l2packet>40) stop 'in_l2packet > 40'
c
      if(in_ndim<=0 .or. in_ndim>3) stop 'in_ndim invalid'
      if(in_nwlg<=0) stop 'in_nwlg invalid'
      if(in_niwlem<=0) stop 'in_niwlem invalid'
c
      if(in_nwlf<=0) stop 'in_nwlf invalid'
      if(in_ncostf<=0) stop 'in_ncostf invalid'
      if(in_nphif<=0) stop 'in_nphif invalid'
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
      if(in_nobbopac) call warn('read_inputpars','ff opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','bb opacity disabled!')
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
      if(trim(in_tempcorrdump)=='off') then
      elseif(trim(in_tempcorrdump)=='all') then
      else
       stop 'in_tempcorrdump invalid'
      endif!}}}
c
c-- override namelist values
      if(.not.in_force3dflux) then!{{{
       if(in_ndim==1 .and. (in_ncostf/=1 .or. in_nphif/=1)) then
        write(6,*)
     &    'override namelist: in_ndim==1 => in_ncostf=1, in_nphif=1'
        in_ncostf = 1
        in_nphif = 1
       elseif(in_ndim==2 .and. in_nphif/=1) then
        write(6,*) 'override namelist: in_ndim==2 => in_nphif=1'
        in_nphif = 1
       endif
      endif
c
      if(in_notimedep) then
       call warn('read_inputpars',
     &   'in_notimedep selected (generate packets thermally)')
       if(trim(in_enoresrv_init)/='off') write(6,*)
     &   "override namelist: in_notimedep => in_enoresrv='off'"
      endif!}}}
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
      endif!}}}
c
      end subroutine parse_inputpars
c
      end module inputparmod
