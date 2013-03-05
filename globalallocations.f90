SUBROUTINE globalallocations

  USE data_mod
  IMPLICIT NONE
  ! Opacity in cm^-1

  ALLOCATE(numcensus(nr))  !# census particles per cell
  ALLOCATE(rarr(nr+1)) !zone edge radii
  ALLOCATE(drarr(nr))  !radial zone length
  ALLOCATE(dr3arr(nr))  !zone volume*3/(4*pi)
  ALLOCATE(Edep(nr))  !energy absorbed by material
  ALLOCATE(Temp(nr))  !cell-centered temperature
  ALLOCATE(Tempb(nr+1))  !cell boundary temperature
  ALLOCATE(Ur(nr))  !equilibrium radiation energy density
  ALLOCATE(rhoarr(nr)) !density
  ALLOCATE(bcoef(nr))  !heat capacity (erg/keV/cm^3)
  ALLOCATE(sigmap(nr)) !Planck opacity
  ALLOCATE(fcoef(nr))  !Fleck factor
  ALLOCATE(Emit(nr))   !Emission energy divided amongst new source particles
  ALLOCATE(Nvol(nr))   !Number source particles in each cell (per time step)
  ALLOCATE(nisource(nr))  ! Nickel gamma source

  !ALLOCATE(gammag(ng)) !group integrated emission source
  ALLOCATE(sigmapg(ng,nr))  !group Planck opacities
  ALLOCATE(sigmargleft(ng,nr))  !left cell edge group Rosseland opacities
  ALLOCATE(sigmargright(ng,nr)) !right ||   ||    ||     ||        ||
  ALLOCATE(EmitProbg(ng,nr))  !Probability of emission in a given zone and group
  ALLOCATE(sigmaL(ng,nr))
  ALLOCATE(sigmaR(ng,nr))
  ALLOCATE(PPL(ng,nr))
  ALLOCATE(PPR(ng,nr))

  ALLOCATE(particles(Npartmax))

  ALLOCATE(Ppick(ng))  !ng=2 for to temp picket fence verification

END SUBROUTINE globalallocations
