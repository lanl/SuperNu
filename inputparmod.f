      module inputparmod
c     ------------------!{{{
      implicit none
************************************************************************
* input parameters
* New parameters need to be added in this routine in four places:
*  1) variable declaration
*  2) namelist
*  3) pointer array
*  4) sanity conditions
************************************************************************
c-- write stdout to file
      character(40) :: in_name = "spn"  !simulation name/title, for post-processing identification
      character(80) :: in_comment = ""  !why did I run this simulation?
c-- parallelization
      integer :: in_nomp = 1  !number of openmp threads
c
c
c-- grid geometry and dimensions
      integer :: in_igeom = 0 !geometry: 1=sph, 2=cyl, 3=car, 11=1Dsph
      integer :: in_ndim(3) = [1, 1, 1]  !number of x-direction cells
      logical :: in_isvelocity = .true.  !switch underlying grid between spatial+static to velocity+expanding
c
c
c-- read input structure file instead of specifying the stucture with input parameters
      logical :: in_voidcorners = .false.  !zero mass in cells outside central sphere in domain
      logical :: in_noreadstruct = .false.
c-- static grid size
      real*8 :: in_lx = 0d0  !spatial length of x-direction
      real*8 :: in_ly = 0d0  !spatial length of y-direction
      real*8 :: in_lz = 0d0  !spatial length of z-direction
c-- parameterized input structure
      real*8 :: in_velout = 0d0  !cm/s, velocity of outer bound
      real*8 :: in_totmass = 0d0  !g
      character(4) :: in_dentype = 'none' ! unif|mass: 'unif' for uniform density, 'mass' for equal mass accross cells
c
c
c-- time step
      real*8 :: in_tsp_tfirst = 0d0  !first point in time evolution [sec]
      real*8 :: in_tsp_tlast = 0d0   !last point in time evolution [sec]
      integer :: in_tsp_nt = 0       !number of time steps
      integer :: in_tsp_itrestart = -1   !restart time step number
        real*8 :: in_tfirst = 0d0
        real*8 :: in_tlast = 0d0
        integer :: in_nt = 0
        integer :: in_ntres = -1
      character(4) :: in_tsp_gridtype = 'lin ' ! line|expo|read: linear, exponential or read-in timestep grid
c
c
c-- group structure
      integer :: in_grp_ng = -1      !number of groups: 0 read wl-grid from file
      integer :: in_grp_ngs = 1      !number of subgroups per opacity group
      real*8 :: in_grp_wlmin =   100d-8  !lower wavelength boundary [cm]
      real*8 :: in_grp_wlmax = 32000d-8  !upper wavelength boundary [cm]
      integer :: in_wldex = 0        !selects group grid from formatted group grid file
        integer :: in_ng = -1
        integer :: in_ngs = 1
        real*8 :: in_wlmin =   100d-8
        real*8 :: in_wlmax = 32000d-8
c
c
c-- outbound flux group and direction bins
      integer :: in_flx_ndim(3) = [0, 1, 1]
      real*8 :: in_flx_wlmin =  1000d-8  !lower wavelength flux boundary [cm]
      real*8 :: in_flx_wlmax = 32000d-8  !upper wavelength flux boundary [cm]
      logical :: in_flx_noobservertime = .false. ! record flux in escape time instead of observer time
c
c
c-- particles
      integer :: in_prt_nmax = 0    !length of particle array
      integer :: in_prt_n2max = -1  !2^n length of particle array
        integer :: in_trn_n2part = -1  !2^n length of particle array
c
c
c-- source
      integer :: in_src_ns = 0     !number of source particles generated per time step (total over all ranks)
      integer :: in_src_n2s = -1   !2^n source particles generated per time step (total over all ranks)
        integer :: in_ns = 0
      integer :: in_src_nsinit = 0   !number of initial particles at in_tsp_tfirst
      integer :: in_src_n2sinit = -1 !2^n number of initial particles at in_tsp_tfirst
        integer :: in_ns0 = 0
      logical :: in_novolsrc = .false.  !switch to turn off any volume source (could be useful for debugs)
      real*8 :: in_srcepwr = 1d0  !source particle number-energy slope, 1 is linear, equal number of packets per erg.
c-- analytic power-law source terms
      real*8 :: in_gas_srccoef = 0d0
      real*8 :: in_gas_srcrpwr = 0d0
      real*8 :: in_gas_srctpwr = 0d0
      real*8 :: in_gas_srctimepwr = 0d0
c-- external source structure
      character(4) :: in_srctype = 'none'   !none|heav|strt|manu|surf: external source structure type
      character(4) :: in_surfsrcloc = 'out ' !in|out|up|down|top|botm: surface source location
      character(4) :: in_surfsrcmu = 'isot' !isot|beam: surface source direction distribution
      integer :: in_nheav = 0   !outer cell bound if heaviside ('heav') source
      real*8 :: in_theav = 0d0 !duration of heaviside source
      real*8 :: in_srcmax = 0d0 !peak source strength
c
c
c-- transport
      logical :: in_trn_errorfatal = .true. !stop on transport error, disable for production runs
      real*8 :: in_tauddmc = 5d0  !number of mean free paths per cell required for DDMC
      real*8 :: in_taulump = 10d0 !number of of mean free paths needed to lump DDMC groups
      logical :: in_puretran = .false. !use IMC only instead of IMC+DDMC hybrid
      logical :: in_isimcanlog = .false. !use analog IMC tally if true
      logical :: in_isddmcanlog = .true. !use analog DDMC tally if true
      logical :: in_trn_noamp = .true.  !disable amplification factor
      real*8 :: in_alpha = 1d0 !time centering control parameter [0,1]
      logical :: in_trn_nolumpshortcut = .false. !disable approximation for large emitlump that sampling outside the lump collapses to the single most likely group
c-- time dependence of in_tauddmc and in_taulump
      character(4) :: in_tauvtime = 'unif' ! unif|incr = constant or limiting (s-curve) to more conservative constant
      logical :: in_ismodimc = .true. !Gentile-Fleck factor switch
c
c
c-- gas grid parameters
      real*8 :: in_gas_gastempinit = 0d0 !non-zero will not read temp from file. units: K
      real*8 :: in_gas_radtempinit = 0d0 !initial radiation temperature.  Use grd_temp by default
        real*8 :: in_consttemp = 0d0
        real*8 :: in_tempradinit = 0d0
c-- analytic heat capacity terms
      real*8 :: in_gas_cvcoef = 1d7 !power law heat capacity coefficient
      real*8 :: in_gas_cvtpwr = 0d0 !power law heat capacity temperature exponent
      real*8 :: in_gas_cvrpwr = 1d0 !power law heat capacity density exponent
c-- debugging
      logical :: in_noeos = .false.     !don't use the EOS
c
c-- physical opacities
      real*8 :: in_opcapgam = .06d0   ![cm^2/g] extinction coefficient for gamma radiation
      logical :: in_noplanckweighting = .false. !disable planck weighting of rosseland opacities within group
      real*8 :: in_opacmixrossel = 0d0 !mix rosseland with planck average, 1=pure rosseland
c-- physical opacities debuging
      logical :: in_nobbopac = .false.  !turn off bound-bound opacity
      logical :: in_nobfopac = .false.  !turn off bound-bound opacity
      logical :: in_noffopac = .false.  !turn off bound-bound opacity
      logical :: in_nothmson = .false.  !turn off thomson scattering
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
      real*8 :: in_gas_sigcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_sigtpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_sigrpwr = 0d0  !power law absorption opacity density exponent
c-- absorption terms:
      real*8 :: in_gas_capcoef = 0d0  !power law absorption opacity coefficient
      real*8 :: in_gas_captpwr = 0d0  !power law absorption opacity temperature exponent
      real*8 :: in_gas_caprpwr = 0d0  !power law absorption opacity density exponent
c
c
c-- output
      logical :: in_grabstdout = .false.  !write stdout to file
      logical :: in_io_dogrdtally = .false. !write transport tallies per grid cell
      logical :: in_nogriddump = .false.  !don't write grid cell variables
      character(4) :: in_opacdump = 'off '    !off|one|each|all: write opacity data to file
      character(4) :: in_pdensdump = 'off '   !off|one|each: write partial densities to file
c     
c-- runtime parameter namelist
      namelist /inputpars/
     & in_name,in_comment,
     & in_nomp,
!grd
     & in_igeom,in_ndim,
     & in_isvelocity,
!str
     & in_voidcorners,in_noreadstruct,
     & in_lx,in_ly,in_lz,
     & in_velout,in_totmass,in_dentype,
!tsp
     & in_tsp_tfirst,in_tsp_tlast,
     & in_tsp_gridtype,in_tsp_nt,in_tsp_itrestart,
     &   in_tfirst,in_tlast,in_nt,in_ntres, !deprec
!grp
     & in_grp_ng,in_grp_ngs,in_grp_wlmin,in_grp_wlmax,
     &   in_ng,in_ngs,in_wlmin,in_wlmax, !deprec
     & in_wldex,
!flx
     & in_flx_ndim,in_flx_wlmin,in_flx_wlmax,in_flx_noobservertime,
!prt
     & in_prt_nmax,in_prt_n2max,
     &   in_trn_n2part, !deprec
!src
     & in_src_ns,in_src_nsinit,
     & in_src_n2s,in_src_n2sinit,
     &   in_ns,in_ns0, !deprec
     & in_srcepwr,
     & in_srctype,in_theav,in_nheav,in_srcmax,
     & in_surfsrcloc,in_surfsrcmu,
     & in_novolsrc,
     & in_gas_srccoef,in_gas_srcrpwr,in_gas_srctpwr,in_gas_srctimepwr,
!trn
     & in_trn_nolumpshortcut,in_trn_errorfatal,in_puretran,in_alpha,
     & in_isimcanlog,in_isddmcanlog,in_trn_noamp,in_ismodimc,
     & in_tauddmc,in_taulump,in_tauvtime,
!gas
     & in_gas_gastempinit,in_gas_radtempinit,
     &   in_consttemp,in_tempradinit, !deprec
     & in_gas_cvcoef,in_gas_cvtpwr,in_gas_cvrpwr,
     & in_noeos,
!phys opac
     & in_opcapgam,
     & in_noplanckweighting,in_opacmixrossel,
     & in_nobbopac,in_nobfopac,in_noffopac,in_nothmson,
!anal opac
     & in_opacanaltype,in_suol,
     & in_suolpick1,in_ldisp1,in_ldisp2,
     & in_gas_sigcoef,in_gas_sigtpwr,in_gas_sigrpwr,
     & in_gas_capcoef,in_gas_captpwr,in_gas_caprpwr,
!io
     & in_grabstdout,
     & in_nogriddump,in_io_dogrdtally,
     & in_opacdump,in_pdensdump
c
c-- pointers
c
      integer,parameter,private :: npointers = 100
c
      type lptr
       logical,pointer :: p
      endtype lptr
      type(lptr) :: in_l(npointers)
c
      type iptr
       integer,pointer :: p
      endtype iptr
      type(iptr) :: in_i(npointers)
c
      type rptr
       real*8,pointer :: p
      endtype rptr
      type(rptr) :: in_r(npointers)
c
      type cptr
       character(4),pointer :: p
      endtype cptr
      type(cptr) :: in_c(npointers)
c
      public
      private inputpars
      save
c!}}}
c
      contains
c
      subroutine inputpar_create_pointers(il,ii,ir,ic)
c     ------------------------------------------------!{{{
      implicit none
************************************************************************
* create pointer arrays (not to be confused with a array pointers) to
* all input parameters.  These arrays are used in mpimod to broadcast
* all input parameters at once.
************************************************************************
      integer,intent(out) :: il,ii,ir,ic
c
c-- init
      il=0
      ii=0
      ir=0
      ic=0
c
      call insertl(in_grabstdout,in_l,il)
      call inserti(in_nomp,in_i,ii)
      call inserti(in_igeom,in_i,ii)
      call inserti(in_ndim(1),in_i,ii)
      call inserti(in_ndim(2),in_i,ii)
      call inserti(in_ndim(3),in_i,ii)
      call insertr(in_lx,in_r,ir)
      call insertr(in_ly,in_r,ir)
      call insertr(in_lz,in_r,ir)
      call inserti(in_flx_ndim(1),in_i,ii)
      call inserti(in_flx_ndim(2),in_i,ii)
      call inserti(in_flx_ndim(3),in_i,ii)
      call insertr(in_flx_wlmin,in_r,ir)
      call insertr(in_flx_wlmax,in_r,ir)
      call insertl(in_flx_noobservertime,in_l,il)
      call insertl(in_nogriddump,in_l,il)
      call insertl(in_io_dogrdtally,in_l,il)
      call insertl(in_noreadstruct,in_l,il)
      call insertl(in_isvelocity,in_l,il)
      call insertl(in_voidcorners,in_l,il)
      call insertr(in_velout,in_r,ir)
      call insertr(in_totmass,in_r,ir)
      call insertc(in_dentype,in_c,ic)
      call insertr(in_gas_gastempinit,in_r,ir)
      call insertr(in_gas_radtempinit,in_r,ir)
      call insertr(in_gas_cvcoef,in_r,ir)
      call insertr(in_gas_cvtpwr,in_r,ir)
      call insertr(in_gas_cvrpwr,in_r,ir)
      call inserti(in_src_ns,in_i,ii)
      call inserti(in_src_nsinit,in_i,ii)
      call inserti(in_src_n2s,in_i,ii)
      call inserti(in_src_n2sinit,in_i,ii)
      call inserti(in_prt_nmax,in_i,ii)
      call inserti(in_prt_n2max,in_i,ii)
      call insertl(in_puretran,in_l,il)
      call insertl(in_trn_nolumpshortcut,in_l,il)
      call insertl(in_trn_errorfatal,in_l,il)
      call insertl(in_isimcanlog,in_l,il)
      call insertl(in_isddmcanlog,in_l,il)
      call insertr(in_tauddmc,in_r,ir)
      call insertr(in_taulump,in_r,ir)
      call insertc(in_tauvtime,in_c,ic)
      call insertr(in_alpha,in_r,ir)
      call insertr(in_tsp_tfirst,in_r,ir)
      call insertr(in_tsp_tlast,in_r,ir)
      call insertc(in_tsp_gridtype,in_c,ic)
      call inserti(in_tsp_nt,in_i,ii)
      call inserti(in_tsp_itrestart,in_i,ii)
      call insertl(in_ismodimc,in_l,il)
      call inserti(in_grp_ng,in_i,ii)
      call inserti(in_grp_ngs,in_i,ii)
      call insertr(in_grp_wlmin,in_r,ir)
      call insertr(in_grp_wlmax,in_r,ir)
      call inserti(in_wldex,in_i,ii)
      call insertr(in_opcapgam,in_r,ir)
      call insertl(in_noplanckweighting,in_l,il)
      call insertr(in_opacmixrossel,in_r,ir)
      call insertc(in_opacanaltype,in_c,ic)
      call insertc(in_suol,in_c,ic)
      call insertr(in_suolpick1,in_r,ir)
      call insertr(in_ldisp1,in_r,ir)
      call insertr(in_ldisp2,in_r,ir)
      call insertr(in_gas_sigcoef,in_r,ir)
      call insertr(in_gas_sigtpwr,in_r,ir)
      call insertr(in_gas_sigrpwr,in_r,ir)
      call insertr(in_gas_capcoef,in_r,ir)
      call insertr(in_gas_captpwr,in_r,ir)
      call insertr(in_gas_caprpwr,in_r,ir)
      call insertr(in_gas_srccoef,in_r,ir)
      call insertr(in_gas_srcrpwr,in_r,ir)
      call insertr(in_gas_srctpwr,in_r,ir)
      call insertr(in_gas_srctimepwr,in_r,ir)
      call insertc(in_srctype,in_c,ic)
      call insertc(in_surfsrcloc,in_c,ic)
      call insertc(in_surfsrcmu,in_c,ic)
      call inserti(in_nheav,in_i,ii)
      call insertr(in_theav,in_r,ir)
      call insertr(in_srcmax,in_r,ir)
      call insertr(in_srcepwr,in_r,ir)
      call insertc(in_opacdump,in_c,ic)
      call insertc(in_pdensdump,in_c,ic)
      call insertl(in_noeos,in_l,il)
      call insertl(in_novolsrc,in_l,il)
      call insertl(in_nobbopac,in_l,il)
      call insertl(in_nobfopac,in_l,il)
      call insertl(in_noffopac,in_l,il)
      call insertl(in_nothmson,in_l,il)
c
      contains
c
      subroutine insertl(par,arr,i)
      implicit none!{{{
      logical,intent(in),target :: par
      type(lptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertl: i>npointers'
      arr(i)%p => par!}}}
      end subroutine insertl
c
      subroutine inserti(par,arr,i)
      implicit none!{{{
      integer,intent(in),target :: par
      type(iptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par!}}}
      end subroutine inserti
c
      subroutine insertr(par,arr,i)
      implicit none!{{{
      real*8,intent(in),target :: par
      type(rptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'insertr: i>npointers'
      arr(i)%p => par!}}}
      end subroutine insertr
c
      subroutine insertc(par,arr,i)
      implicit none!{{{
      character(4),intent(in),target :: par
      type(cptr),intent(inout) :: arr(npointers)
      integer,intent(inout) :: i
      i = i + 1
      if(i>npointers) stop 'inserti: i>npointers'
      arr(i)%p => par!}}}
      end subroutine insertc
!}}}
      end subroutine inputpar_create_pointers
c
c
      subroutine read_inputpars
c     -------------------------!{{{
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
!}}}
      end subroutine read_inputpars
c
c
c
      subroutine deprecate_inputpars
c     -------------------------------!{{{
***********************************************************************
* Warn about deprecated input variables and error stop if new and old
* variable are both set.

* The purpose is to help the user transition away from deprecated variables
***********************************************************************
      if(in_ns/=0) stop 'DEPRECATED: in_ns => in_src_ns'
      if(in_ns0/=0) stop 'DEPRECATED: in_ns0 => in_src_nsinit'
      if(in_trn_n2part/=-1) stop
     &  'DEPRECATED: in_trn_n2part => in_prt_n2max'
c-- tsp
      if(in_tfirst/=0d0) stop 'DEPRECATED: in_tfirst => in_tsp_tfirst'
      if(in_tlast/=0d0) stop 'DEPRECATED: in_tlast => in_tsp_tlast'
      if(in_nt/=0) then
       if(in_tsp_nt/=0) then
        stop 'DEPRECATED: in_nt => in_tsp_nt'
       else
        in_tsp_nt = in_nt
        write(6,*) 'DEPRECATED: in_nt => in_tsp_nt'
       endif
      endif
      if(in_ntres/=-1) then
       if(in_tsp_itrestart/=-1) then
        stop 'DEPRECATED: in_ntres => in_tsp_itrestart'
       else
        in_tsp_itrestart = in_ntres
        write(6,*) 'DEPRECATED: in_ntres => in_tsp_itrestart'
       endif
      endif
c-- gas
      if(in_consttemp/=0d0) then
       if(in_gas_gastempinit/=0d0) then
        stop 'DEPRECATED: in_consttemp => in_gas_gastempinit'
       else
        in_gas_gastempinit = in_consttemp
        write(6,*) 'DEPRECATED: in_consttemp => in_gas_gastempinit'
       endif
      endif
      if(in_tempradinit/=0d0) then
       if(in_gas_radtempinit/=0d0) then
        stop 'DEPRECATED: in_tempradinit => in_gas_radtempinit'
       else
        in_gas_radtempinit = in_tempradinit
        write(6,*) 'DEPRECATED: in_tempradinit => in_gas_radtempinit'
       endif
      endif
c-- grp
      if(in_ng/=-1) then
       if(in_grp_ng/=-1) then
        stop 'DEPRECATED: in_ng => in_grp_ng'
       else
        in_grp_ng = in_ng
        write(6,*) 'DEPRECATED: in_ng => in_grp_ng'
       endif
      endif
      if(in_ngs/=1) then
       if(in_grp_ngs/=1) then
        stop 'DEPRECATED: in_ngs => in_grp_ngs'
       else
        in_grp_ngs = in_ngs
        write(6,*) 'DEPRECATED: in_ngs => in_grp_ngs'
       endif
      endif
      if(in_wlmin/=100d-8) then
       if(in_grp_wlmin/=100d-8) then
        stop 'DEPRECATED: in_wlmin => in_grp_wlmin'
       else
        in_grp_wlmin = in_wlmin
        write(6,*) 'DEPRECATED: in_wlmin => in_grp_wlmin'
       endif
      endif
      if(in_wlmax/=32000d-8) then
       if(in_grp_wlmax/=32000d-8) then
        stop 'DEPRECATED: in_wlmax => in_grp_wlmax'
       else
        in_grp_wlmax = in_wlmax
        write(6,*) 'DEPRECATED: in_wlmax => in_grp_wlmax'
       endif
      endif
      write(6,*)
!}}}
      end subroutine deprecate_inputpars
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
c-- dump namelist to stdout
      write(6,*) 'namelist read:'
      write(6,nml=inputpars)
      write(6,*)
c
c-- write simulation name to file
      open(4,file='output.name',iostat=istat)
      if(istat/=0) stop 'parse_inputpars: open output.name error'
      write(4,'(a)') trim(in_name)
      close(4)
c
c-- check input parameter validity
      if(in_nomp<0) stop 'in_nomp invalid'!{{{
      if(in_nomp==0 .and. nmpi>1) stop 'no in_nomp==0 in mpi mode'
c
      if(any(in_ndim<1)) stop 'in_ndim invalid'
c
      select case(in_igeom)
      case(1)
       if(in_srctype=='surf'.and.in_surfsrcloc/='out') stop
     &   'in_srctype and in_surfsrcloc invalid'
      case(2)
       if(in_srctype=='surf' .and.
     &      any((/'in  ','top ','botm'/)==in_surfsrcloc))
     &      stop 'in_srctype and in_surfsrcloc invalid'
      case(3)
      case(11)
       if(in_ndim(2)>1 .or. in_ndim(3)>1) stop 'in_ndim invalid'
       if(in_srctype=='surf'.and.in_surfsrcloc/='out') stop
     &   'in_srctype and in_surfsrcloc invalid'
       if(in_flx_ndim(2)/=1) stop 'in_flx_ndim(2) inval'
       if(in_flx_ndim(3)/=1) stop 'in_flx_ndim(3) inval'
      case default
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
      if(in_nogriddump .and. in_io_dogrdtally) stop
     &   'dogridtally and !griddump'
c
      if(in_voidcorners.and.in_igeom==1) stop 'voidcorners && igeom=1'
      if(in_voidcorners.and.in_igeom==11) stop 'voidcorners && igeom=11'
c
      if(in_grp_ng<0) stop 'in_grp_ng invalid'
      if(in_grp_ng==0 .and. in_wldex<1) stop 'in_wldex invalid'
      if(in_grp_ngs<=0) stop 'in_grp_ngs invalid'
      if(in_grp_wlmin<0) stop 'in_grp_wlmin invalid'
      if(in_grp_wlmax<=in_grp_wlmin) stop 'in_grp_wlmax invalid'
c
      if(in_src_ns<=0 .eqv. in_src_n2s<0) stop
     &  'use in_src_ns or in_src_n2s'
      if(in_src_nsinit>0 .and. in_src_n2sinit>=0) stop
     &  'use in_src_nsinit or in_src_n2sinit'
      if(in_prt_nmax>0) then
       if(in_prt_n2max>-1) stop 'use in_prt_nmax OR in_prt_n2max'
       if(in_prt_nmax<max(int(in_src_ns,8),int(2,8)**in_src_n2s)) stop
     &   'in_prt_nmax too small'
      else
       if(in_prt_nmax>0) stop 'use in_prt_nmax OR in_prt_n2max'
       if(int(2,8)**in_prt_n2max < max(int(in_src_ns,8),
     &   int(2,8)**in_src_n2s)) stop 'in_prt_nmax too small'
      endif
c
      if(in_alpha>1d0 .or. in_alpha<0d0) stop 'in_alpha invalid'
      if(in_taulump<=.05d0*in_tauddmc) stop 'in_taulump too small' !don't let scattering dominate
c
      if(in_noreadstruct) then
        if(in_totmass<=0d0) stop 'in_totmass <= 0'
        if(in_noreadstruct.and.in_novolsrc.and.in_gas_cvcoef<=0d0)
     &     stop 'in_noreadstruct && in_novolsrc && in_gas_cvcoef<=0'
        if(in_dentype=='none') stop 'noreadstruct & in_dentype==none'
      else
        if(in_dentype/='none') stop '!noreadstruct & in_dentype/=none'
      endif
c
c-- temp init
      if(in_gas_gastempinit<0d0) stop 'in_gas_gastempinit < 0'
      if(in_gas_radtempinit<0d0) stop 'in_gas_radtempinit < 0'
c
c-- timestepping
      select case(in_tsp_gridtype)
      case('lin ')
       if(in_tsp_nt<=0) stop 'in_tsp_gridtype==line & in_tsp_nt<=0'
      case('expo')
       if(in_tsp_nt<=0) stop 'in_tsp_gridtype==expo & in_tsp_nt<=0'
       if(in_tsp_tfirst<=0d0) stop
     &   'in_tsp_gridtype==expo & in_tsp_tfirst<=0'
      case('read')
       if(in_tsp_nt/=0) stop 'in_tsp_gridtype==read & in_tsp_nt/=0'
      case default
       stop 'in_tsp_gridtype invalid'
      end select

      if(in_tsp_tlast<in_tsp_tfirst) stop 'in_tsp_tlast invalid'
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
      case('surf')
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
c
      if(in_nobbopac) call warn('read_inputpars','bb opacity disabled!')
      if(in_nobfopac) call warn('read_inputpars','bf opacity disabled!')
      if(in_noffopac) call warn('read_inputpars','ff opacity disabled!')
      if(in_nothmson) call warn('read_inputpars','Thomson disabled')
c
      select case(in_opacdump)
      case('off ')
      case('one ')
      case('each')
      case('all ')
       call warn('read_inputpars',
     &   "in_opacdump=='all' will generate a big data file!")
      case default
       stop 'in_opacdump invalid'
      end select
c
      if(in_opacmixrossel<0d0 .or. in_opacmixrossel>1d0) then
       stop 'in_opacmixrossel invalid'
      endif
c
      select case(in_pdensdump)
      case('off ')
      case('one ')
      case('each')
      case default
       stop 'in_pdensdump invalid'
      end select
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
      if(in_grabstdout) then
       write(0,'(1x,a,2i5,i7)') 'nmpi,in_nomp,#threads        :',
     &   nmpi,in_nomp,nmpi*in_nomp
      endif!}}}!}}}
c
      end subroutine parse_inputpars
c
c
c
      subroutine provide_inputpars(nmpi)
c     -------------------------------------!{{{
      use physconstmod
      use particlemod
      use sourcemod
      use transportmod
      use timestepmod
      use fluxmod
      use groupmod
      use gasmod
      use gridmod
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* Distribute the input parameter values to the respective modules.
* This needs to be called AFTER the values are bcast.
************************************************************************
      integer :: mpart
      integer :: ns,nsinit
c
      mpart = int(int(2,8)**in_prt_n2max/nmpi) !max number of particles
      mpart = max(mpart,in_prt_nmax/nmpi)
      prt_npartmax = mpart
c
      ns = int(int(2,8)**in_src_n2s/nmpi)
      ns = max(ns,in_src_ns/nmpi)
      nsinit = int(int(2,8)**in_src_n2sinit/nmpi)
      nsinit = max(nsinit,in_src_nsinit/nmpi)
      src_ns = ns
      src_ninit = nsinit
c
      tsp_gridtype = in_tsp_gridtype
      tsp_nt     = in_tsp_nt
      tsp_itrestart  = in_tsp_itrestart
      tsp_tfirst = in_tsp_tfirst
      tsp_tlast  = in_tsp_tlast
c
      gas_gastempinit = in_gas_gastempinit
      gas_radtempinit = in_gas_radtempinit
      !gas_sigcoef = in_gas_sigcoef
      !gas_sigtpwr = in_gas_sigtpwr
      !gas_sigrpwr = in_gas_sigrpwr
      !gas_capcoef = in_gas_capcoef
      !gas_captpwr = in_gas_captpwr
      !gas_caprpwr = in_gas_caprpwr
      !gas_cvcoef  = in_gas_cvcoef
      !gas_cvtpwr  = in_gas_cvtpwr
      !gas_cvrpwr  = in_gas_cvrpwr
      !gas_srccoef = in_gas_srccoef      
      !gas_srcrpwr = in_gas_srcrpwr      
      !gas_srctpwr = in_gas_srctpwr      
      !gas_srctimepwr = in_gas_srctimepwr
c
      trn_isimcanlog = in_isimcanlog
      trn_isddmcanlog = in_isddmcanlog
      trn_tauddmc = in_tauddmc
      trn_taulump = in_taulump
      trn_tauvtime = in_tauvtime
      trn_nolumpshortcut = in_trn_nolumpshortcut
      trn_errorfatal = in_trn_errorfatal
      trn_noampfact = in_trn_noamp
c
      !io_dogrdtally = in_io_dogrdtally
c
      flx_ndim  = in_flx_ndim
      flx_wlmin = in_flx_wlmin
      flx_wlmax = in_flx_wlmax
      flx_noobservertime = in_flx_noobservertime
c
      grp_ng    = in_grp_ng
      grp_ngs   = in_grp_ngs
      grp_wlmin = in_grp_wlmin
      grp_wlmax = in_grp_wlmax
c
      grd_igeom = in_igeom
      grd_nx    = in_ndim(1)
      grd_ny    = in_ndim(2)
      grd_nz    = in_ndim(3)
      grd_isvelocity = in_isvelocity
c!}}}
      end subroutine provide_inputpars
c
      end module inputparmod
