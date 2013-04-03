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


  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
  real*8,allocatable :: gas_dwl(:) !(gas_ng) wavelength grid bin width
  real*8,allocatable :: gas_cap(:,:) !(gas_nr,gas_ng) Line+Cont extinction coeff
  real*8,allocatable :: gas_sig(:) !(gas_nr) scattering coefficient
  real*8,allocatable :: gas_capgam(:) !(gas_nr) Gamma ray gray opacity
!
!-- temperature structure history
  real*8,allocatable :: gas_temphist(:,:) !(gas_nr,tim_nt)

  real*8 :: gas_lr = 0
  real*8 :: gas_l0 = 0  !innermost static radius
  logical :: gas_isvelocity
  logical :: gas_isshell  !domain is shell, innermost radius not zero
  logical :: gas_novolsrc !no external volume source (e.g. radioactivity)
  logical :: gas_isanalgrp  !switch to use analytic_opacity
  logical :: gas_isanalsrc  !switch to use analytic_source
  integer :: gas_velno, gas_velyes
  real*8 :: gas_templ0=0 !surface temperature at innermost radius
  real*8 :: gas_sigcoef=0  !analytic opacity power law coefficient
  real*8 :: gas_sigtpwr=0  !analytic opacity power law temperature exponent
  real*8 :: gas_sigrpwr=0  !analytic opacity power law density exponent
  real*8 :: gas_cvcoef=0  !analytic heat capacity power law coefficient
  real*8 :: gas_cvtpwr=0  !analytic heat capacity power law temperature exponent
  real*8 :: gas_cvrpwr=0  !analytic heat capacity power law density exponent

  character(4) :: gas_grptype = 'grey' !analytic opacity dependence on group.
  !is used with power law to create group opacities each timestep (see 
  !inputparmod for possible values).
  character(4) :: gas_suol = 'tsta' !if gas_grptype='pick', sets picket
  !magnitudes with values from cases in literature (Su&Olson 1999).
  
  real*8 :: gas_ldisp !if gas_grptype='line',
  !ratio of strong line group opacity strength to weak group

  character(4) :: gas_srctype = 'heav' !analytic external source dependence
  !on space (and group if gas_srctype='manu')
  integer :: gas_nheav = 0 !outer cell bound of external heaviside ('heav') source
  real*8 :: gas_theav = 0d0 !duration of heaviside source
  real*8 :: gas_srcmax = 0d0 !peak strength (ergs/cm^3/s) of external source

  real*8 :: gas_emat

  real*8, dimension(:), allocatable :: gas_rarr  !(gas_nr+1)
  real*8, dimension(:), allocatable :: gas_drarr !(gas_nr)
  real*8, dimension(:), allocatable :: gas_edep, gas_sigmap, gas_fcoef !(gas_nr)
  real*8, dimension(:), allocatable :: gas_tempb  !(gas_nr+1)
  real*8, dimension(:,:), allocatable :: gas_sigmapg, gas_sigmargleft, gas_sigmargright, gas_emitprobg  !(gas_ng,gas_nr)
  real*8, dimension(:,:), allocatable :: gas_sigmal, gas_ppl, gas_sigmar, gas_ppr  !(gas_ng,gas_nr)
  !Ryan W.: External sources (currently used for analytic_source):
  real*8, dimension(:,:), allocatable :: gas_exsource !(gas_ng,gas_nr)

  type gas_secondary
    sequence
    integer :: nvol, nvolex
    real*8 :: dr3_34pi, tempkev, emit
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
    gas_velout = in_velout
    !
    gas_isvelocity = in_isvelocity
    gas_isshell = in_isshell
    gas_novolsrc = in_novolsrc
    gas_isanalgrp = in_isanalgrp
    !inner edge radius (if in_isshell):
    gas_l0 = in_l0
    !inner edge temp:
    gas_templ0 = in_templ0
    !power law grey opacity input:
    gas_sigcoef = in_sigcoef
    gas_sigtpwr = in_sigtpwr
    gas_sigrpwr = in_sigrpwr
    !power law heat capacity input:
    gas_cvcoef = in_cvcoef
    gas_cvtpwr = in_cvtpwr
    gas_cvrpwr = in_cvrpwr
    !group type:
    gas_grptype = in_grptype
    !picket fence input:
    gas_suol = in_suol
    gas_ppick(1) = in_suolpick1
    gas_ppick(2) = 1d0-in_suolpick1
    !group line disparity:
    gas_ldisp = in_ldisp
    !external analytic source input:
    gas_srctype = in_srctype
    gas_isanalsrc = in_isanalsrc
    gas_theav = in_theav
    gas_nheav = in_nheav
    gas_srcmax = in_srcmax

    ! Setting velocity option
    if (in_isvelocity.eqv..true.) then
       gas_velyes = 1
       gas_velno = 0
    else
       gas_velyes = 0
       gas_velno = 1
    endif

!-- primary
    allocate(gas_numcensus(gas_nr))  !# census prt_particles per cell
    allocate(gas_rarr(gas_nr+1)) !zone edge radii
    allocate(gas_drarr(gas_nr))  !radial zone length
    allocate(gas_edep(gas_nr))  !energy absorbed by material
    allocate(gas_sigmap(gas_nr)) !Planck opacity (gray)
    allocate(gas_fcoef(gas_nr))  !Fleck factor
    allocate(gas_sigmapg(gas_ng,gas_nr))  !group Planck opacities
    allocate(gas_emitprobg(gas_ng,gas_nr))  !Probability of emission in a given zone and group
    allocate(gas_sigmal(gas_ng,gas_nr))
    allocate(gas_sigmar(gas_ng,gas_nr))
    allocate(gas_ppl(gas_ng,gas_nr))  !can potentially be removed
    allocate(gas_ppr(gas_ng,gas_nr))  !can potentially be removed
    !allocate(gas_ppick(gas_ng))  !gas_ng=2 for to temp picket fence verification

!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
    allocate(gas_cap(gas_nr,gas_ng)) !Line+Cont extinction coeff

!-- secondary
    allocate(gas_vals2(gas_nr))
    allocate(gas_sig(gas_nr))
    allocate(gas_capgam(gas_nr))

    allocate(gas_temphist(gas_nr,nt))
    allocate(gas_dwl(gas_ng)) !wavelength grid bin width

    allocate(gas_tempb(gas_nr+1))  !cell boundary temperature

    allocate(gas_sigmargleft(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
    allocate(gas_sigmargright(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||

    allocate(gas_exsource(gas_ng,gas_nr))

  end subroutine gasgrid_init

end module gasgridmod
