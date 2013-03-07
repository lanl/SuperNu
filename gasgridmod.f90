MODULE gasgridmod

  IMPLICIT NONE

  INTEGER :: gas_nr = 0
  INTEGER :: gas_ng = 0
  REAL*8 :: gas_lr = 0
  LOGICAL :: gas_isvelocity
  INTEGER :: gas_velno, gas_velyes

  INTEGER :: gas_nelem = 0

  REAL*8 :: gas_nidecay = 1.73*(1.6022e-6)  !erg/s/g  !this value is used in the first iteration

  REAL*8 :: gas_emat

  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_rarr  !(gas_nr+1)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_drarr, gas_dr3arr  !(gas_nr)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_edep, gas_temp, gas_sigmap, gas_fcoef, gas_ur !(gas_nr)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_tempb  !(gas_nr+1)
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_rhoarr, gas_bcoef, gas_emit, gas_nisource  !(gas_nr)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: gas_sigmapg, gas_sigmargleft, gas_sigmargright, gas_emitprobg  !(gas_nr,gas_ng)
  REAL*8, DIMENSION(:,:), ALLOCATABLE :: gas_sigmal, gas_ppl, gas_sigmar, gas_ppr  !(gas_nr,gas_ng)

  ! Picket-fence probabilities
  REAL*8, DIMENSION(:), ALLOCATABLE :: gas_ppick

  INTEGER, DIMENSION(:), ALLOCATABLE :: gas_numcensus, gas_nvol !(gas_nr)

  save

  CONTAINS


  SUBROUTINE gasgrid_init(nr,ng,lr,isvelocity)
!-------------------------------
    INTEGER,intent(in) :: nr, ng
    REAL*8,intent(in) :: lr
    LOGICAL,intent(in) :: isvelocity
    gas_nr = nr
    gas_ng = ng
    gas_lr = lr
    gas_isvelocity = isvelocity

    ALLOCATE(gas_numcensus(gas_nr))  !# census prt_particles per cell
    ALLOCATE(gas_nvol(gas_nr))   !Number source prt_particles in each cell (per tsp_time step)
    ALLOCATE(gas_rarr(gas_nr+1)) !zone edge radii
    ALLOCATE(gas_drarr(gas_nr))  !radial zone length
    ALLOCATE(gas_dr3arr(gas_nr))  !zone volume*3/(4*pc_pi)
    ALLOCATE(gas_edep(gas_nr))  !energy absorbed by material
    ALLOCATE(gas_temp(gas_nr))  !cell-centered temperature
    ALLOCATE(gas_tempb(gas_nr+1))  !cell boundary temperature
    ALLOCATE(gas_ur(gas_nr))  !equilibrium radiation energy density
    ALLOCATE(gas_rhoarr(gas_nr)) !density
    ALLOCATE(gas_bcoef(gas_nr))  !heat capacity (erg/keV/cm^3)
    ALLOCATE(gas_sigmap(gas_nr)) !Planck opacity
    ALLOCATE(gas_fcoef(gas_nr))  !Fleck factor
    ALLOCATE(gas_emit(gas_nr))   !Emission energy divided amongst new source prt_particles
    ALLOCATE(gas_nisource(gas_nr))  ! Nickel gamma source

    !ALLOCATE(gas_gammag(gas_ng)) !group integrated emission source
    ALLOCATE(gas_sigmapg(gas_ng,gas_nr))  !group Planck opacities
    ALLOCATE(gas_sigmargleft(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
    ALLOCATE(gas_sigmargright(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||
    ALLOCATE(gas_emitprobg(gas_ng,gas_nr))  !Probability of emission in a given zone and group
    ALLOCATE(gas_sigmal(gas_ng,gas_nr))
    ALLOCATE(gas_sigmar(gas_ng,gas_nr))
    ALLOCATE(gas_ppl(gas_ng,gas_nr))
    ALLOCATE(gas_ppr(gas_ng,gas_nr))
    ALLOCATE(gas_ppick(gas_ng))  !gas_ng=2 for to temp picket fence verification
  END SUBROUTINE gasgrid_init

END MODULE gasgridmod
