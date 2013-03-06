SUBROUTINE globalallocations

  USE particlemod
  USE gasgridmod
  USE physconstmod

  IMPLICIT NONE
  ! Opacity in cm^-1

  ALLOCATE(numcensus(gas_nr))  !# census particles per cell
  ALLOCATE(rarr(gas_nr+1)) !zone edge radii
  ALLOCATE(drarr(gas_nr))  !radial zone length
  ALLOCATE(dr3arr(gas_nr))  !zone volume*3/(4*pi)
  ALLOCATE(Edep(gas_nr))  !energy absorbed by material
  ALLOCATE(Temp(gas_nr))  !cell-centered temperature
  ALLOCATE(Tempb(gas_nr+1))  !cell boundary temperature
  ALLOCATE(Ur(gas_nr))  !equilibrium radiation energy density
  ALLOCATE(rhoarr(gas_nr)) !density
  ALLOCATE(bcoef(gas_nr))  !heat capacity (erg/keV/cm^3)
  ALLOCATE(sigmap(gas_nr)) !Planck opacity
  ALLOCATE(fcoef(gas_nr))  !Fleck factor
  ALLOCATE(Emit(gas_nr))   !Emission energy divided amongst new source particles
  ALLOCATE(Nvol(gas_nr))   !Number source particles in each cell (per time step)
  ALLOCATE(nisource(gas_nr))  ! Nickel gamma source

  !ALLOCATE(gammag(gas_ng)) !group integrated emission source
  ALLOCATE(sigmapg(gas_ng,gas_nr))  !group Planck opacities
  ALLOCATE(sigmargleft(gas_ng,gas_nr))  !left cell edge group Rosseland opacities
  ALLOCATE(sigmargright(gas_ng,gas_nr)) !right ||   ||    ||     ||        ||
  ALLOCATE(EmitProbg(gas_ng,gas_nr))  !Probability of emission in a given zone and group
  ALLOCATE(sigmaL(gas_ng,gas_nr))
  ALLOCATE(sigmaR(gas_ng,gas_nr))
  ALLOCATE(PPL(gas_ng,gas_nr))
  ALLOCATE(PPR(gas_ng,gas_nr))

  ALLOCATE(particles(prt_Npartmax))

  ALLOCATE(gas_ppick(gas_ng))  !gas_ng=2 for to temp picket fence verification

END SUBROUTINE globalallocations
