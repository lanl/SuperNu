module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: gas_nr = 0
  integer :: gas_ng = 0

!--lumping index
  integer :: gas_epslump
!
  real*8 :: gas_velout = 0d0 !outer boundary velocity
  real*8 :: gas_v0 = 0d0 !inner boundary velocity (if in_isshell)

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
  real*8,allocatable :: gas_cap(:,:) !(gas_ng,gas_nr) Line+Cont extinction coeff
  real*8,allocatable :: gas_sig(:) !(gas_nr) scattering coefficient
!-(rev. 121): edge scattering coefficient (a secondary quantity)
  real*8,allocatable :: gas_sigbl(:) !(gas_nr), left
  real*8,allocatable :: gas_sigbr(:) !(gas_nr), right
!---------------------------------------------------------------
  real*8,allocatable :: gas_capgam(:) !(gas_nr) Gamma ray gray opacity
!
!-- temperature structure history
  real*8,allocatable :: gas_temphist(:,:) !(gas_nr,tim_nt)

  real*8 :: gas_lr = 0
  real*8 :: gas_l0 = 0  !innermost static radius
  logical :: gas_isvelocity = .false.
  logical :: gas_isshell = .false.  !domain is shell, innermost radius not zero
  logical :: gas_novolsrc = .false. !no external volume source (e.g. radioactivity)
  logical :: gas_depestimate = .true. !if true uses deposition estimator to update temperature
  real*8 :: gas_templ0=0 !surface temperature at innermost radius
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

  real*8, dimension(:), allocatable :: gas_rarr   !(gas_nr+1), left cell edge values
  real*8, dimension(:), allocatable :: gas_drarr  !(gas_nr)
  real*8, dimension(:), allocatable :: gas_curvcent !(gas_nr), multiplied by tauddmc for mfp threshold
  real*8, dimension(:), allocatable :: gas_edep, gas_siggrey, gas_fcoef !(gas_nr)
  real*8, dimension(:), allocatable :: gas_tempb  !(gas_nr+1), interpolated temperatures (keV)
  real*8, dimension(:), allocatable :: gas_rhob   !(gas_nr+1), interpolated densities
  real*8, allocatable :: gas_emitprob(:,:)           !(gas_ng,gas_nr)
  real*8, allocatable :: gas_ppl(:,:), gas_ppr(:,:)  !(gas_ng,gas_nr)
!-- temporary array used to compute leakage opacities
  real*8, allocatable :: gas_caprosl(:,:), gas_caprosr(:,:)  !(gas_ng,gas_nr)
!-- leakage opacities
  real*8, allocatable :: gas_opacleakl(:,:), gas_opacleakr(:,:) !(gas_ng,gas_nr)
  
  real*8, allocatable :: gas_eraddens(:,:) !(gas_ng,gas_nr)
!-- old Planck opacity for BDF-2 method
  real*8, allocatable :: gas_siggreyold(:) !(gas_nr)
!
!-- outbound grouped luminosity
  real*8, allocatable :: gas_luminos(:) !(gas_ng)
!
!---
  integer, allocatable :: gas_nvol(:) !(gas_nr) number of thermal source particles generated per cell
  integer, allocatable :: gas_nvolex(:) !(gas_nr) number of external source particles generated per cell
  integer, allocatable :: gas_nvolinit(:) !(gas_nr) number of initial (t=tfirst) particles per cell
!  
  real*8, allocatable :: gas_emit(:) !(gas_nr) amount of fictitious thermal energy emitted per cell in a time step
  real*8, allocatable :: gas_emitex(:,:) !(gas_ng,gas_nr) amount of external energy emitted per cell per group in a time step
  real*8, allocatable :: gas_evolinit(:,:) !(gas_ng,gas_nr) amount of initial energy per cell per group
!
  real*8, allocatable :: gas_temp(:)
!-- rtw: attempting method upgrade to order 2 backward difference formula in time (rev. 365)
!-- rtw: must store a temperature value from previous time step (only on rank 0)
  real*8, allocatable :: gas_tempold(:)
!---
!
!
  type gas_secondary
    sequence
    
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
!-- opacity invalidity flag
       logical :: opdirty=.true. !opacity needs recalculation
!-- energy reservoir
       real*8 :: engdep     !energy deposited by gamma rays
!-- material energy (temperature) source (may be manufactured), rev>244
       real*8 :: matsrc = 0d0
  end type gas_secondary
  type(gas_secondary),pointer :: gas_vals2(:)

  ! Picket-fence probabilities
  real*8 :: gas_ppick(2)

  integer, dimension(:), allocatable :: gas_numcensus !(gas_nr)

  save

  contains


  subroutine gasgrid_init(nt,ng)
!-------------------------------------------------------
    use inputparmod
    implicit none
!
    integer,intent(in) :: nt,ng
!
    gas_nr = in_nr
    gas_ng = ng
    !
    gas_isvelocity = in_isvelocity
    gas_isshell = in_isshell
    gas_novolsrc = in_novolsrc
    gas_depestimate = in_depestimate
    !inner edge radius (if in_isshell):
    !inner edge temp:
    gas_templ0 = in_templ0
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
    !lumping index
    gas_epslump = in_epslump

!-- primary
    allocate(gas_numcensus(gas_nr))  !# census prt_particles per cell
    allocate(gas_rarr(gas_nr+1)) !zone edge radii
    allocate(gas_drarr(gas_nr))  !radial zone length
    allocate(gas_curvcent(gas_nr))  ! gas_curvcent*tauddmc=mean fee path threshold for IMC-DDMC heuristic
    allocate(gas_edep(gas_nr))  !energy absorbed by material
    allocate(gas_siggrey(gas_nr)) !Planck opacity (gray)
!-- rtw: using old gas_siggrey in bdf2 Fleck factor calculation
    allocate(gas_siggreyold(gas_nr)) !Planck opacity (gray)
!
!- Ryan W.: using power law to calculate gas_sig (similar to Planck opacity)
    allocate(gas_sig(gas_nr))    !grey scattering opacity
!----------------------------------------------------------------
    allocate(gas_fcoef(gas_nr))  !Fleck factor
    allocate(gas_emitprob(gas_ng,gas_nr))  !Probability of emission in a given zone and group
    allocate(gas_opacleakl(gas_ng,gas_nr))
    allocate(gas_opacleakr(gas_ng,gas_nr))
    allocate(gas_ppl(gas_ng,gas_nr))  !can potentially be removed
    allocate(gas_ppr(gas_ng,gas_nr))  !can potentially be removed
    allocate(gas_eraddens(gas_ng,gas_nr))  !radiation energy density in tsp_dt per group
    allocate(gas_luminos(gas_ng))  !outbound grouped luminosity array
!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
    allocate(gas_cap(gas_ng,gas_nr)) !Line+Cont extinction coeff
!--Ryan W: values below were formerly secondary (rev 183)
    allocate(gas_nvol(gas_nr))
    allocate(gas_nvolex(gas_nr))
    allocate(gas_nvolinit(gas_nr))
    allocate(gas_emit(gas_nr))
    allocate(gas_emitex(gas_ng,gas_nr))
    allocate(gas_evolinit(gas_ng,gas_nr))
!-- secondary
    allocate(gas_vals2(gas_nr))
    allocate(gas_capgam(gas_nr))
!--Ryan W: added edge scattering opacities for leakage coefficients (rev. 121)
    allocate(gas_sigbl(gas_nr))
    allocate(gas_sigbr(gas_nr))
!------------------------------------
    allocate(gas_temphist(gas_nr,nt))

    allocate(gas_temp(gas_nr))  !cell average temperature
    allocate(gas_tempb(gas_nr+1))  !cell boundary temperature
    allocate(gas_tempold(gas_nr))  !temp used in BDF-2 update
    allocate(gas_rhob(gas_nr+1))   !cell boundary density
    
    allocate(gas_caprosl(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
    allocate(gas_caprosr(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||

  end subroutine gasgrid_init

end module gasgridmod
