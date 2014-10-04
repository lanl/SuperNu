module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: gas_nx = 0
  integer :: gas_ny = 0
  integer :: gas_nz = 0
  integer :: gas_ng = 0

  real*8 :: gas_lx = 0
  real*8 :: gas_ly = 0
  real*8 :: gas_lz = 0

!
  real*8 :: gas_velout = 0d0 !outer boundary velocity

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
  real*8,allocatable :: gas_cap(:,:,:,:) !(gas_ng,gas_nx,gas_ny,gas_nz) Line+Cont extinction coeff
  real*8,allocatable :: gas_sig(:,:,:) !(gas_nx,gas_ny,gas_nz) scattering coefficient
!---------------------------------------------------------------
  real*8,allocatable :: gas_capgam(:,:,:) !(gas_nx,gas_ny,gas_nz) Gamma ray gray opacity
!
!-- temperature structure history
  real*8,allocatable :: gas_temphist(:,:,:,:) !(gas_nx,gas_ny,gas_nz,tim_nt)
  real*8,allocatable :: gas_temppreset(:,:,:,:) !(gas_nx,gas_ny,gas_nz,tim_nt)
  logical :: gas_isvelocity = .false.
  logical :: gas_novolsrc = .false. !no external volume source (e.g. radioactivity)
!-(rev. 121)
  real*8 :: gas_sigcoefs=0  !analytic scattering opacity power law coefficient
  real*8 :: gas_sigtpwrs=0  !analytic scattering opacity power law temperature exponent
  real*8 :: gas_sigrpwrs=0  !analytic scattering opacity power law density exponent
!-
  real*8 :: gas_sigcoef=0  !analytic absorption opacity power law coefficient
  real*8 :: gas_sigtpwr=0  !analytic absorption opacity power law temperature exponent
  real*8 :: gas_sigrpwr=0  !analytic absorption opacity power law density exponent
  real*8 :: gas_cvcoef=0  !analytic heat capacity power law coefficient
  real*8 :: gas_cvtpwr=0  !analytic heat capacity power law temperature exponent
  real*8 :: gas_cvrpwr=0  !analytic heat capacity power law density exponent

  character(4) :: gas_opacanaltype = 'grey' !analytic opacity dependence on group.
  !is used with power law to create group opacities each timestep (see 
  !inputparmod for possible values).
  character(4) :: gas_suol = 'tsta' !if gas_opacanaltype='pick', sets picket
  !magnitudes with values from cases in literature (Su&Olson 1999).
  
  real*8 :: gas_ldisp1 = 0d0 !if gas_opacanaltype='line',
  real*8 :: gas_ldisp2 = 0d0 !if gas_opacanaltype-'line'
  !ratio of strong line group opacity strength to weak group

  character(4) :: gas_srctype = 'none' !analytic external source dependence
  !on space (and group if gas_srctype='manu')
  integer :: gas_nheav = 0 !outer cell bound of external heaviside ('heav') source
  real*8 :: gas_theav = 0d0 !duration of heaviside source
  real*8 :: gas_srcmax = 0d0 !peak strength (ergs*s^2/cm^3 if isvelocity, else ergs/cm^3/s) of external source
  real*8 :: gas_tempradinit = 0d0 !initial radiation temperature
!
!-- energy conservation check quantities
  real*8 :: gas_eext = 0d0 !time-integrated input energy from external source
  real*8 :: gas_emat = 0d0 !material energy
  real*8 :: gas_erad = 0d0 !census radiation energy
  real*8 :: gas_eleft= 0d0 !left (inward) leaked energy from domain
  real*8 :: gas_eright=0d0 !right (outward) leaked energy from domain
  real*8 :: gas_etot = 0d0 !total source energy added per time step
  real*8 :: gas_evelo= 0d0 !total energy change to rad field from fluid
  real*8 :: gas_eerror= 0d0 !error in integral problem energy
!-- average quantities used for energy check with MPI
  real*8 :: gas_eextav = 0d0
  real*8 :: gas_eveloav = 0d0
!--
  real*8 :: gas_esurf

  real*8,allocatable :: gas_xarr(:)   !(gas_nx+1), left cell edge values
  real*8,allocatable :: gas_yarr(:)   !(gas_ny+1), left cell edge values
  real*8,allocatable :: gas_zarr(:)   !(gas_nz+1), left cell edge values
  real*8,allocatable :: gas_dxarr(:)  !(gas_nx)
  real*8,allocatable :: gas_dyarr(:)  !(gas_ny)
  real*8,allocatable :: gas_dzarr(:)  !(gas_nz)
  real*8,allocatable :: gas_edep(:,:,:)   !(gas_nx,gas_ny,gas_nz)
  real*8,allocatable :: gas_siggrey(:,:,:)!(gas_nx,gas_ny,gas_nz)
  real*8,allocatable :: gas_fcoef(:,:,:)  !(gas_nx,gas_ny,gas_nz)
  real*8,allocatable :: gas_emitprob(:,:,:,:)           !(gas_ng,gas_nx,gas_ny,gas_nz)
!-- leakage opacities
  real*8,allocatable :: gas_opacleak(:,:,:,:) !(2,gas_nx,gas_ny,gas_nz)
  
  real*8,allocatable :: gas_eraddens(:,:,:) !(gas_nx,gas_ny,gas_nz)
!-- old Planck opacity for BDF-2 method
  real*8,allocatable :: gas_siggreyold(:,:,:) !(gas_nx,gas_ny,gas_nz)
!
!-- outbound grouped luminosity
  real*8,allocatable :: gas_luminos(:) !(gas_ng)
!-- sampled devation of group luminosity
  real*8,allocatable :: gas_lumdev(:) !(gas_ng)
!-- number of escaped particles per group
  integer,allocatable :: gas_lumnum(:) !(gas_ng)
!
!---
  integer,allocatable :: gas_nvol(:,:,:) !(gas_nx,gas_ny,gas_nz) number of thermal source particles generated per cell
  integer,allocatable :: gas_nvolex(:,:,:) !(gas_nx,gas_ny,gas_nz) number of external source particles generated per cell
  integer,allocatable :: gas_nvolinit(:,:,:) !(gas_nx,gas_ny,gas_nz) number of initial (t=tfirst) particles per cell
!  
  real*8,allocatable :: gas_emit(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of fictitious thermal energy emitted per cell in a time step
  real*8,allocatable :: gas_emitex(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of external energy emitted per cell per group in a time step
  real*8,allocatable :: gas_evolinit(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of initial energy per cell per group
!
  real*8,allocatable :: gas_temp(:,:,:) !(gas_nx,gas_ny,gas_nz)
!---
!
  integer,allocatable :: gas_methodswap(:,:,:) !(gas_nx,gas_ny,gas_nz) number of IMC-DDMC method changes per cell per time step
!
  type gas_secondary
    real*8 :: eraddens
    real*8 :: ur, rho, bcoef, nisource
       !real*8 :: temp       !gcell temperature
       real*8 :: volr       !gcell volume [rout=1 units]
       real*8 :: vol        !gcell volume [cm^3]
       real*8 :: volcrp     !effective volume (of linked rgrid cells) [cm^3]
       real*8 :: mass       !gcell mass
       real*8 :: mass0fr(-2:gas_nelem) = 0d0  !initial mass fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       real*8 :: natom      !gcell # atoms
       real*8 :: natom1fr(-2:gas_nelem) = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       real*8 :: natom0fr(-2:2) = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
       real*8 :: nelec=1d0  !gcell # electrons per atom
!-- energy reservoir
       real*8 :: engdep     !energy deposited by gamma rays
!-- material energy (temperature) source (may be manufactured), rev>244
       real*8 :: matsrc = 0d0
  end type gas_secondary
  type(gas_secondary),pointer :: gas_vals2(:,:,:) !(gas_nx,gas_ny,gas_nz)

  ! Picket-fence probabilities
  real*8 :: gas_ppick(2)

  integer,allocatable :: gas_numcensus(:,:,:) !(gas_nx,gas_ny,gas_nz)

  private read_temp_preset

  save

  contains


  subroutine gasgrid_init(nt,ng)
!-------------------------------------------------------
    use inputparmod
    implicit none
    integer,intent(in) :: nt,ng
!
    logical :: lexist
!
    gas_nx = in_nx
    gas_ny = in_ny
    gas_nz = in_nz
    gas_ng = ng
    !
    gas_isvelocity = in_isvelocity
    gas_novolsrc = in_novolsrc
    !power law heat capacity input:
    gas_cvcoef = in_cvcoef
    gas_cvtpwr = in_cvtpwr
    gas_cvrpwr = in_cvrpwr
    !power law scattering opacity input:
    gas_sigcoefs = in_sigcoefs
    gas_sigtpwrs = in_sigtpwrs
    gas_sigrpwrs = in_sigrpwrs
    !power law absorption opacity input:
    gas_sigcoef = in_sigcoef
    gas_sigtpwr = in_sigtpwr
    gas_sigrpwr = in_sigrpwr
    !group type:
    gas_opacanaltype = in_opacanaltype
    !picket fence input:
    gas_suol = in_suol
    gas_ppick(1) = in_suolpick1
    gas_ppick(2) = 1d0-in_suolpick1
    !group line disparities (strengths):
    gas_ldisp1 = in_ldisp1
    gas_ldisp2 = in_ldisp2
    !external analytic source input:
    gas_srctype = in_srctype
    gas_theav = in_theav
    gas_nheav = in_nheav
    gas_srcmax = in_srcmax

!-- primary
    allocate(gas_xarr(gas_nx+1)) !zone edge x position
    allocate(gas_yarr(gas_ny+1)) !zone edge y position
    allocate(gas_zarr(gas_nz+1)) !zone edge z position
    allocate(gas_dxarr(gas_nx))  !radial zone length
    allocate(gas_dyarr(gas_ny))  !radial zone length
    allocate(gas_dzarr(gas_nz))  !radial zone length
    allocate(gas_edep(gas_nx,gas_ny,gas_nz))   !energy absorbed by material
    allocate(gas_siggrey(gas_nx,gas_ny,gas_nz)) !Planck opacity (gray)
!-- rtw: using old gas_siggrey in bdf2 Fleck factor calculation
    allocate(gas_siggreyold(gas_nx,gas_ny,gas_nz)) !Planck opacity (gray)
!
!- Ryan W.: using power law to calculate gas_sig (similar to Planck opacity)
    allocate(gas_sig(gas_nx,gas_ny,gas_nz))    !grey scattering opacity
!----------------------------------------------------------------
    allocate(gas_fcoef(gas_nx,gas_ny,gas_nz))  !Fleck factor
    allocate(gas_emitprob(gas_ng,gas_nx,gas_ny,gas_nz))  !Probability of emission in a given zone and group
    allocate(gas_opacleak(gas_ng,gas_nx,gas_ny,gas_nz))
    allocate(gas_eraddens(gas_nx,gas_ny,gas_nz))  !radiation energy density in tsp_dt per group
    allocate(gas_luminos(gas_ng))  !outbound grouped luminosity array
    allocate(gas_lumdev(gas_ng))
    allocate(gas_lumnum(gas_ng))
!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
    allocate(gas_cap(gas_ng,gas_nx,gas_ny,gas_nz)) !Line+Cont extinction coeff
!--Ryan W: values below were formerly secondary (rev 183)
    allocate(gas_nvol(gas_nx,gas_ny,gas_nz))
    allocate(gas_nvolex(gas_nx,gas_ny,gas_nz))
    allocate(gas_nvolinit(gas_nx,gas_ny,gas_nz))
    allocate(gas_emit(gas_nx,gas_ny,gas_nz))
    allocate(gas_emitex(gas_nx,gas_ny,gas_nz))
    allocate(gas_evolinit(gas_nx,gas_ny,gas_nz))
!-- secondary
    allocate(gas_vals2(gas_nx,gas_ny,gas_nz))
    allocate(gas_capgam(gas_nx,gas_ny,gas_nz))
!------------------------------------
    allocate(gas_temphist(gas_nx,gas_ny,gas_nz,nt))

    allocate(gas_temp(gas_nx,gas_ny,gas_nz))  !cell average temperature
    allocate(gas_methodswap(gas_nx,gas_ny,gas_nz))
    allocate(gas_numcensus(gas_nx,gas_ny,gas_nz))  !# census prt_particles per cell
!
!-- read preset temperature profiles
    inquire(file='input.temp',exist=lexist)
    if(lexist) call read_temp_preset
  end subroutine gasgrid_init


  subroutine read_temp_preset
!-----------------------------!{{{
    use timestepmod
    integer :: istat
!
    open(4,file='input.temp',status='old',iostat=istat)
    if(istat/=0) stop 'rd_temp_preset: no file: input.temp'
!-- alloc and read
    allocate(gas_temppreset(gas_nx,gas_ny,gas_nz,tsp_nt))
    read(4,*,iostat=istat) gas_temppreset
    if(istat/=0) stop 'rd_temp_preset: file too short: input.temp'
!-- check EOF
    read(4,*,iostat=istat)
    if(istat==0) stop 'rd_temp_preset: file too long: input.temp'
    close(4)
    write(6,*) 'rd_temp_preset: custom temp profiles read successfully'
!}}}
  end subroutine read_temp_preset


end module gasgridmod
