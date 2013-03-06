SUBROUTINE globalallocations

  USE particlemod
  USE gasgridmod
  USE physconstmod

  IMPLICIT NONE
  ! Opacity in cm^-1

  ALLOCATE(prt_numcensus(gas_nr))  !# census prt_particles per cell
  ALLOCATE(gas_rarr(gas_nr+1)) !zone edge radii
  ALLOCATE(gas_drarr(gas_nr))  !radial zone length
  ALLOCATE(gas_dr3arr(gas_nr))  !zone volume*3/(4*pi)
  ALLOCATE(gas_edep(gas_nr))  !energy absorbed by material
  ALLOCATE(gas_temp(gas_nr))  !cell-centered temperature
  ALLOCATE(gas_tempb(gas_nr+1))  !cell boundary temperature
  ALLOCATE(gas_ur(gas_nr))  !equilibrium radiation energy density
  ALLOCATE(gas_rhoarr(gas_nr)) !density
  ALLOCATE(gas_bcoef(gas_nr))  !heat capacity (erg/keV/cm^3)
  ALLOCATE(gas_sigmap(gas_nr)) !Planck opacity
  ALLOCATE(gas_fcoef(gas_nr))  !Fleck factor
  ALLOCATE(gas_emit(gas_nr))   !Emission energy divided amongst new source prt_particles
  ALLOCATE(prt_nvol(gas_nr))   !Number source prt_particles in each cell (per tsp_time step)
  ALLOCATE(gas_nisource(gas_nr))  ! Nickel gamma source

  !ALLOCATE(gammag(gas_ng)) !group integrated emission source
  ALLOCATE(gas_sigmapg(gas_ng,gas_nr))  !group Planck opacities
  ALLOCATE(gas_sigmargleft(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
  ALLOCATE(gas_sigmargright(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||
  ALLOCATE(gas_emitprobg(gas_ng,gas_nr))  !Probability of emission in a given zone and group
  ALLOCATE(gas_sigmal(gas_ng,gas_nr))
  ALLOCATE(gas_sigmar(gas_ng,gas_nr))
  ALLOCATE(gas_ppl(gas_ng,gas_nr))
  ALLOCATE(gas_ppr(gas_ng,gas_nr))

  ALLOCATE(prt_particles(prt_npartmax))

  ALLOCATE(Ppick(gas_ng))  !gas_ng=2 for to temp picket fence verification

END SUBROUTINE globalallocations
