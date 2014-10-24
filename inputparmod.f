      module inputparmod
c     ------------------
      use physconstmod
      implicit none
************************************************************************
* input parameters
************************************************************************
c-- write stdout to file
      character(80) :: in_comment = "" !why did I run this simulation?
      logical :: in_grab_stdout = .false. !write stdout to file
c-- parallelization
      integer :: in_nomp = 1       !number of openmp threads
c
c-- grid geometry and dimensions
      integer :: in_igeom = 0 !geometry: 1=1Dsph, 2=2Dcyl, 3=3Dcar
      integer :: in_ndim(3) = [1, 1, 1] !number of x-direction cells

      real*8 :: in_lx = 0d0  !spatial length of x-direction
      real*8 :: in_ly = 0d0  !spatial length of y-direction
      real*8 :: in_lz = 0d0  !spatial length of z-direction
c
c-- do read input structure file instead of specifying the stucture with input parameters
c==================
      logical :: in_noreadstruct = .false.
c-- special grid
      logical :: in_isvelocity = .true.  !switch underlying grid between spatial+static to velocity+expanding
c-- specify the atmospheric stratification
      real*8 :: in_velout = 0d0  !cm/s, velocity of outer bound
      real*8 :: in_totmass = 0d0  !g
      character(4) :: in_dentype = 'unif' ! unif|mass: 'unif' for uniform density, 'mass' for equal mass accross cells
c============
c
c-- temperature parameters
      real*8 :: in_consttemp = 0d0 !non-zero will not read temp from file. units: K
      real*8 :: in_tempradinit = 0d0 !initial radiation temperature.  Use gas_temp by default
c
c-- analytic heat capacity terms
      real*8 :: in_cvcoef = 1d0 !power law heat capacity coefficient
      real*8 :: in_cvtpwr = 0d0 !power law heat capacity temperature exponent
      real*8 :: in_cvrpwr = 0d0 !power law heat capacity density exponent
c
c-- random number generator
      integer :: in_seed = 1984117 !starting point of random number generator
c
c-- particles
      integer :: in_ns = 0    !number of source particles generated per time step (total over all ranks)
      integer :: in_ns0 = 0   !number of initial particles at in_tfirst
      integer :: in_npartmax = 0 !total number of particles allowed PER MPI RANK
      logical :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      logical :: in_isimcanlog = .false. !use analog IMC tally if true
      logical :: in_isddmcanlog = .true. !use analog DDMC tally if true
      real*8 :: in_tauddmc = 5d0 !number of mean free paths per cell required for DDMC
      real*8 :: in_taulump = 10d0 !number of of mean free paths needed to lump DDMC groups
c-- time dependence of in_tauddmc and in_taulump
      character(4) :: in_tauvtime = 'unif' ! unif|incr = constant or limiting (s-curve) to more conservative constant
c
      real*8 :: in_alpha = 1d0 !time centering control parameter [0,1]
c
c-- time step
      real*8 :: in_tfirst = 0d0 !first point in time evolution
      real*8 :: in_tlast = 0d0  !last point in time evolution
      integer :: in_nt = 0      !number of time steps.  <0 means read timeline from input.tsp_time
      integer :: in_ntres = -1   !restart time step number
      logical :: in_norestart = .true.
      logical :: in_ismodimc=.true. !Gentile-Fleck factor switch
c
c
c-- group structure
      integer :: in_ng = -1      !number of groups: 0 uses in_wldex
      integer :: in_ngs = 0      !>0 number of subgroups per opacity group
                                 ! 0 non-subgridded physical_opacities
                                 !<0 target subgrid resolution for automatic subgridding -> lambda/(Delta lambda) = abs(in_ngs)
      integer :: in_wldex = 1    !if in_iswlread = t, selects group grid from formatted group grid file
      real*8  :: in_wlmin =   100e-8 !lower wavelength boundary [cm]
      real*8  :: in_wlmax = 32000e-8 !upper wavelength boundary [cm]
c
c
c-- physical opacities
      real*8 :: in_opcapgam = .06d0   ![cm^2/g] extinction coefficient for gamma radiation
      real*8 :: in_epsline = 1d0      !line absorption fraction (the rest is scattering)
      logical :: in_noplanckweighting = .false. !disable planck weighting of rosseland opacities within group
      real*8 :: in_opacmixrossel = 1d0 !mix rosseland with planck average, 1=pure rosseland
c-- test switches
      logical :: in_nobbopac = .false.    !turn off bound-bound opacity
      logical :: in_nobfopac = .false.    !turn off bound-bound opacity
      logical :: in_noffopac = .false.    !turn off bound-bound opacity
      logical :: in_nothmson = .false.    !turn off thomson scattering
c
c
c-- analytic opacities
      character(4) :: in_opacanaltype = 'none'    !none|grey|mono|pick|line: group opacity structure type
c-- picket fence specific group structure
      character(4) :: in_suol = 'tsta'    !tsta|tstb|tstc: Su&Olson picket fence (pick) test cases 
      real*8 :: in_suolpick1 = 1d0  !in [0,1]: probability of being at first picket
c-- line specific group structure
      real*8 :: in_ldisp1 = 1d0  !loosely speaking, the analytic odd group line strength
      real*8 :: in_ldisp2 = 1d0  !loosely speaking, the analytic even group line strength
c-- scattering terms:
      real*8 :: in_sigcoefs = 0d0 !power law absorption opacity coefficient
      real*8 :: in_sigtpwrs = 0d0 !power law absorption opacity temperature exponent
      real*8 :: in_sigrpwrs = 0d0 !power law absorption opacity density exponent
c-- absorption terms:
      real*8 :: in_sigcoef = 0d0 !power law absorption opacity coefficient
      real*8 :: in_sigtpwr = 0d0 !power law absorption opacity temperature exponent
      real*8 :: in_sigrpwr = 0d0 !power law absorption opacity density exponent
c
c
c-- external source structure
      character(4) :: in_srctype='none'   !none|heav|strt|manu: external source structure type
      integer :: in_nheav = 0   !outer cell bound if heaviside ('heav') source
      real*8 :: in_theav = 0d0 !duration of heaviside source
      real*8 :: in_srcmax = 0d0 !peak source strength (ergs/cm^3/s)
c-- disable all sources
      logical :: in_novolsrc = .true.  !switch to turn off any volume source (could be useful for debugs)
c
c-- misc
      character(4) :: in_opacdump = 'off'    !off|one|each|all: write opacity data to file
      character(4) :: in_pdensdump = 'off'   !off|one|each: write partial densities to file
c     
c-- runtime parameter namelist
      namelist /inputpars/
     & in_igeom,in_ndim,
     & in_isvelocity,in_novolsrc,
     & in_lx,in_ly,in_lz,
     & in_ng,in_ngs,in_wldex,in_wlmin,in_wlmax,
     & in_totmass,in_velout,
     & in_consttemp,
     & in_seed,in_ns,in_ns0,in_npartmax,in_puretran,in_alpha,
     & in_tfirst,in_tlast,in_nt,in_ntres,
     & in_grab_stdout,in_nomp,
     & in_opcapgam,in_epsline,in_nobbopac,in_nobfopac,
     & in_noffopac,in_nothmson,in_noplanckweighting,in_opacmixrossel,
     & in_opacdump,in_pdensdump,
     & in_sigcoefs,in_sigtpwrs,in_sigrpwrs,
     & in_sigcoef,in_sigtpwr,in_sigrpwr,
     & in_cvcoef,in_cvtpwr,in_cvrpwr,
     & in_opacanaltype,in_suol,
     & in_suolpick1, in_ldisp1, in_ldisp2,
     & in_srctype, in_theav, in_nheav, in_srcmax,
     & in_isimcanlog, in_isddmcanlog,
     & in_tauddmc, in_dentype, in_noreadstruct,
     & in_norestart, in_taulump, in_tauvtime,
     & in_tempradinit, in_ismodimc,
     & in_comment
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
       open(6,file='fort.6',action='write',status='replace',recl=3000,
     &   iostat=istat) !write stdout to file
       if(istat/=0) stop 'parse_inputpars: open fort.6 error'
       call banner
      endif!}}}
c
c-- dump namelist to stdout
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)
      write(6,*)
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'!{{{
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(any(in_ndim<1)) stop 'in_ndim invalid'
      select case(in_igeom)
      case(:0)
       stop 'in_igeom invalid'
      case(1)
       if(in_ndim(2)>1 .or. in_ndim(3)>1) stop 'in_ndim invalid'
      case(2)
       if(in_ndim(3)>1) stop 'in_ndim invalid'
      case(4:)
       stop 'in_igeom invalid'
      endselect
c
      if(in_isvelocity) then
       if(in_lx>0d0) stop 'vel grid: use in_velout, not in_lx'
       if(in_ly>0d0) stop 'vel grid: use in_velout, not in_ly'
       if(in_lz>0d0) stop 'vel grid: use in_velout, not in_lz'
       if(in_velout<=0d0 .and. in_noreadstruct) stop
     &   'vel grid: use in_velout, not in_lx,in_ly,in_lz'
      else
       if(in_lx<=0d0) stop 'static grid: use in_lx, not in_velout'
       if(in_ly<=0d0) stop 'static grid: use in_ly, not in_velout'
       if(in_lz<=0d0) stop 'static grid: use in_lz, not in_velout'
       if(in_velout>0d0) 
     &      stop 'static grid: use in_lx,in_ly,in_lz, not in_velout'
      endif
c
      if(in_ng==0) then
       if(in_wldex<1) stop 'in_wldex invalid'
      elseif(in_ng<0) then
       stop 'in_ng invalid'
      endif
      if(in_wlmin<0) stop 'in_wlmin invalid'
      if(in_wlmax<=in_wlmin) stop 'in_wlmax invalid'
c
      if(in_ns<=0) stop 'in_ns <= 0'
      if(in_ns0<0) stop 'in_ns0 < 0'
      if(in_npartmax<=in_ns+in_ns0) stop 'in_npartmax invalid'
      if(in_alpha>1d0 .or. in_alpha<0d0) stop 'in_alpha invalid'
      if(in_taulump<in_tauddmc) stop 'in_taulump<in_tauddmc'
c
      if(in_totmass<=0d0 .and. in_noreadstruct) stop 'in_totmass <= 0'
c
c-- temp init
      if(in_consttemp<0d0) stop 'in_consttemp < 0'
      if(in_tempradinit<0d0) stop 'in_tempradinit < 0'
c
      if(in_nt==0) stop 'in_nt invalid'
      if(in_tfirst<0d0) stop 'in_tfirst invalid'
      if(in_tlast<in_tfirst) stop 'in_tlast invalid'
c
c-- special grid
      if(.not.in_noreadstruct) then
       if(in_velout/=0d0) stop 'velout incomp. with struct'
       if(in_totmass/=0d0) stop 'totmass incomp. with struct'
      endif
c
      select case(in_srctype)
      case('none')
      case('heav')
      case('strt')
      case('manu')
      case default
       stop 'in_srctype unknown'
      end select
c
      select case(in_opacanaltype)
      case('none')
c--R.W.: condition under case(pick) supposed to be here? (rev 243)
         if(in_nobbopac.and.in_nobfopac.and.in_noffopac)
     &        stop 'no phys opac + in_opacanaltype==none'
      case('grey')
      case('mono')
      case('pick')
C$$$       if(.not.in_nobbopac) stop 'no phys opac + in_grptyp==none'
C$$$       if(.not.in_nobfopac) stop 'no phys opac + in_grptyp==none'
C$$$       if(.not.in_noffopac) stop 'no phys opac + in_grptyp==none'
      case('line')
      case default
       stop 'in_opacanaltype unknown'
      end select
c-- disallow physical opacities when analytic opacities are selected
      if(in_opacanaltype/='none' .and. .not.(in_nobbopac .and.
     &  in_nobfopac .and. in_noffopac)) then
       stop 'in_no??opac: no physical opacities allowed with anal opac'
      endif
c
      if(in_opcapgam<0d0) stop 'in_opcapgam invalid'
      if(in_epsline<0d0 .or. in_epsline>1d0) stop 'in_epsline invalid'
c
      if(in_nobbopac) call warn('read_inputpars','bb opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','ff opacity disabled!')
      if(in_nothmson) call warn('read_inputpars','Thomson disabled')
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
      if(in_opacmixrossel<0d0 .or. in_opacmixrossel>1d0) then
       stop 'in_opacmixrossel invalid'
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
