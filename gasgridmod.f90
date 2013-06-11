module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: gas_nr = 0
  integer :: gas_ng = 0

  real*8 :: gas_velout = 0d0 !outer boundary velocity
  real*8 :: gas_v0 = 0d0 !inner boundary velocity (if in_isshell)

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
  real*8,allocatable :: gas_dwl(:) !(gas_ng) wavelength grid bin width
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

  real*8 :: gas_emat = 0d0

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

  !Ryan W.: External sources (currently used for analytic_source):
  real*8, dimension(:,:), allocatable :: gas_exsource !(gas_ng,gas_nr)

  type gas_secondary
    sequence
    integer :: nvol, nvolex
    real*8 :: emit, eraddens
    real*8 :: ur, rho, bcoef, nisource
       real*8 :: temp       !gcell temperature
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
  end type gas_secondary
  type(gas_secondary),pointer :: gas_vals2(:)

  real*8 :: gas_eleft, gas_eright, gas_erad, gas_eint, gas_etot, gas_esurf

  ! Picket-fence probabilities
  real*8 :: gas_ppick(2)

  integer, dimension(:), allocatable :: gas_numcensus !(gas_nr)

  save

  contains


  subroutine gasgrid_init(nt)
!-------------------------------------------------------
    use inputparmod
    implicit none
!
    integer,intent(in) :: nt
!
    gas_nr = in_nr
    gas_ng = in_ng
    gas_lr = in_lr
    gas_v0 = in_v0
    !
    gas_isvelocity = in_isvelocity
    gas_isshell = in_isshell
    gas_novolsrc = in_novolsrc
    !inner edge radius (if in_isshell):
    gas_l0 = in_l0
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


!-- primary
    allocate(gas_numcensus(gas_nr))  !# census prt_particles per cell
    allocate(gas_rarr(gas_nr+1)) !zone edge radii
    allocate(gas_drarr(gas_nr))  !radial zone length
    allocate(gas_curvcent(gas_nr))  ! gas_curvcent*tauddmc=mean fee path threshold for IMC-DDMC heuristic
    allocate(gas_edep(gas_nr))  !energy absorbed by material
    allocate(gas_siggrey(gas_nr)) !Planck opacity (gray)
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
!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
    allocate(gas_cap(gas_ng,gas_nr)) !Line+Cont extinction coeff

!-- secondary
    allocate(gas_vals2(gas_nr))
    !allocate(gas_sig(gas_nr))
    allocate(gas_capgam(gas_nr))
!--Ryan W: added edge scattering opacities for leakage coefficients (rev. 121)
    allocate(gas_sigbl(gas_nr))
    allocate(gas_sigbr(gas_nr))
!------------------------------------
    allocate(gas_temphist(gas_nr,nt))
    allocate(gas_dwl(gas_ng)) !wavelength grid bin width

    allocate(gas_tempb(gas_nr+1))  !cell boundary temperature
    allocate(gas_rhob(gas_nr+1))   !cell boundary density
    
    allocate(gas_caprosl(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
    allocate(gas_caprosr(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||

    allocate(gas_exsource(gas_ng,gas_nr))

  end subroutine gasgrid_init

end module gasgridmod
