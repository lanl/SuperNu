MODULE gasgridmod

  IMPLICIT NONE
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: gas_nr = 0
  integer :: gas_ng = 0

  REAL*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid
  REAL*8,allocatable :: gas_dwl(:) !(gas_ng) wavelength grid bin width
  REAL*8,allocatable :: gas_cap(:,:) !(gas_nr,gas_ng) Line+Cont extinction coeff
!
!-- temperature structure history
  REAL*8,allocatable :: gas_temphist(:,:) !(gas_nr,tim_nt)

  REAL*8 :: gas_lr = 0
  LOGICAL :: gas_isvelocity
  INTEGER :: gas_velno, gas_velyes

  REAL*8 :: gas_nidecay = 1.73*(1.6022e-6)  !erg/s/g  !this value is used in the first iteration

  REAL*8 :: gas_emat

  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_rarr  !(gas_nr+1)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_drarr !(gas_nr)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_edep, gas_sigmap, gas_fcoef !(gas_nr)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_tempb  !(gas_nr+1)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: gas_sigmapg, gas_sigmargleft, gas_sigmargright, gas_emitprobg  !(gas_nr,gas_ng)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: gas_sigmal, gas_ppl, gas_sigmar, gas_ppr  !(gas_nr,gas_ng)

  TYPE gas_secondary
    SEQUENCE
    INTEGER :: nvol
    REAL*8 :: dr3_34pi, tempkev, emit
    REAL*8 :: ur, rho, bcoef, nisource
       REAL*8 :: temp       !gcell temperature
       REAL*8 :: volr       !gcell volume [rout=1 units]
       REAL*8 :: vol        !gcell volume [cm^3]
       REAL*8 :: volcrp     !effective volume (of linked rgrid cells) [cm^3]
       REAL*8 :: mass       !gcell mass
       REAL*8 :: mass0fr(-2:gas_nelem) = 0d0  !initial mass fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       REAL*8 :: natom      !gcell # atoms
       REAL*8 :: natom1fr(-2:gas_nelem) = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       REAL*8 :: natom0fr(-2:2) = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
       REAL*8 :: nelec=1d0  !gcell # electrons per atom
!-- opacity invalidity flag
       LOGICAL :: opdirty=.true. !opacity needs recalculation
!-- energy reservoir
       REAL*8 :: engdep     !energy deposited by gamma rays
  END TYPE gas_secondary
  TYPE(gas_secondary),POINTER :: gas_vals2(:)

  REAL*8 :: gas_eleft, gas_eright, gas_erad, gas_einit, gas_einp, gas_eint, gas_etot, gas_esurf

  ! Picket-fence probabilities
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_ppick

  INTEGER, DIMENSION(:), ALLOCATABLE :: gas_numcensus !(gas_nr)

  save

  CONTAINS


  SUBROUTINE gasgrid_init(nr,ng,nt,lr,isvelocity)
!-------------------------------
    INTEGER,intent(in) :: nr,ng,nt
    REAL*8,intent(in) :: lr
    LOGICAL,intent(in) :: isvelocity
    gas_nr = nr
    gas_ng = ng
    gas_lr = lr
    gas_isvelocity = isvelocity

!-- primary
    ALLOCATE(gas_numcensus(gas_nr))  !# census prt_particles per cell
    ALLOCATE(gas_rarr(gas_nr+1)) !zone edge radii
    ALLOCATE(gas_drarr(gas_nr))  !radial zone length
    ALLOCATE(gas_edep(gas_nr))  !energy absorbed by material
    ALLOCATE(gas_sigmap(gas_nr)) !Planck opacity (gray)
    ALLOCATE(gas_fcoef(gas_nr))  !Fleck factor
    ALLOCATE(gas_sigmapg(gas_ng,gas_nr))  !group Planck opacities
    ALLOCATE(gas_emitprobg(gas_ng,gas_nr))  !Probability of emission in a given zone and group
    ALLOCATE(gas_sigmal(gas_ng,gas_nr))
    ALLOCATE(gas_sigmar(gas_ng,gas_nr))
    ALLOCATE(gas_ppl(gas_ng,gas_nr))  !can potentially be removed
    ALLOCATE(gas_ppr(gas_ng,gas_nr))  !can potentially be removed
    ALLOCATE(gas_ppick(gas_ng))  !gas_ng=2 for to temp picket fence verification

    allocate(gas_wl(gas_ng)) !wavelength grid
    allocate(gas_cap(gas_nr,gas_ng)) !Line+Cont extinction coeff

!-- secondary
    allocate(gas_vals2(gas_nr))

    allocate(gas_temphist(gas_nr,nt))
    allocate(gas_dwl(gas_ng)) !wavelength grid bin width

    ALLOCATE(gas_tempb(gas_nr+1))  !cell boundary temperature

    ALLOCATE(gas_sigmargleft(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
    ALLOCATE(gas_sigmargright(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||

  END SUBROUTINE gasgrid_init

END MODULE gasgridmod
