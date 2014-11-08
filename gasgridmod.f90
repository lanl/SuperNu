module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: dd_ncell=0

  integer :: gas_ng=0

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid

  logical :: gas_novolsrc = .false. !no external volume source (e.g. radioactivity)
!-(rev. 121)
  real*8 :: gas_sigcoefs=0  !analytic scattering opacity power law coefficient
  real*8 :: gas_sigtpwrs=0  !analytic scattering opacity power law temperature exponent
  real*8 :: gas_sigrpwrs=0  !analytic scattering opacity power law density exponent
!-
  real*8 :: gas_sigcoef=0  !analytic absorption opacity power law coefficient
  real*8 :: gas_sigtpwr=0  !analytic absorption opacity power law temperature exponent
  real*8 :: gas_sigrpwr=0  !analytic absorption opacity power law density exponent

  ! Picket-fence probabilities
  real*8 :: gas_ppick(2)

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


!-- temperature structure history
  real*8,allocatable :: dd_temppreset(:,:) !(ncell,tim_nt)

!
!
!-- DOMAIN DECOMPOSITION
!=======================
  real*8,allocatable :: dd_temp(:)       !(ncell)
  real*8,allocatable :: dd_eraddens(:)
  real*8,allocatable :: dd_ur(:)
  real*8,allocatable :: dd_rho(:)
  real*8,allocatable :: dd_bcoef(:)
  real*8,allocatable :: dd_nisource(:)
  real*8,allocatable :: dd_vol(:)        !cell volume [cm^3]
  real*8,allocatable :: dd_mass(:)       !cell mass [g]
  real*8,allocatable :: dd_natom(:)      !cell number of atoms
  real*8,allocatable :: dd_nelec(:)      !cell number of electrons per atom
  real*8,allocatable :: dd_natom1fr(:,:) !(-2:gas_nelem,ncell)  !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
  real*8,allocatable :: dd_natom0fr(:,:) !(-2:2,ncell)  !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
!-- mate,allocatablerial energy (temperature) source (may be manufactured), rev>244
  real*8,allocatable :: dd_matsrc(:)
!
  real*8,allocatable :: dd_edep(:)

!== DD copies
!-- Probability of emission in a given zone and group
  real*8,allocatable :: dd_emitprob(:,:) !(gas_ng,ncell)
!-- Line+Cont extinction coeff
  real*8,allocatable :: dd_cap(:,:) !(gas_ng,ncell)
!-- leakage opacities
! real*8,allocatable :: dd_opacleak(:,:) !(6,ncell)
!-- scattering coefficient
  real*8,allocatable :: dd_sig(:) !(ncell)
!-- Gamma ray gray opacity
  real*8,allocatable :: dd_capgam(:) !(ncell)
!-- Planck opacity (gray)
  real*8,allocatable :: dd_siggrey(:)!(ncell)
!-- Fleck factor
  real*8,allocatable :: dd_fcoef(:)  !(ncell)

  real*8 :: dd_emat = 0d0 !material energy

!  
  real*8,allocatable :: dd_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
  real*8,allocatable :: dd_emitex(:) !(ncell) amount of external energy emitted per cell per group in a time step
  real*8,allocatable :: dd_evolinit(:) !(ncell) amount of initial energy per cell per group

  save

  contains


  subroutine gasgrid_init(ltalk,ncell)
!-------------------------------------
    use inputparmod
    implicit none
    logical,intent(in) :: ltalk
    integer,intent(in) :: ncell

    integer :: n
    logical :: lexist

    dd_ncell = ncell

    gas_novolsrc = in_novolsrc
!-- power law scattering opacity input:
    gas_sigcoefs = in_sigcoefs
    gas_sigtpwrs = in_sigtpwrs
    gas_sigrpwrs = in_sigrpwrs
!-- power law absorption opacity input:
    gas_sigcoef = in_sigcoef
    gas_sigtpwr = in_sigtpwr
    gas_sigrpwr = in_sigrpwr
!-- group type:
    gas_opacanaltype = in_opacanaltype
!-- picket fence input:
    gas_suol = in_suol
    gas_ppick(1) = in_suolpick1
    gas_ppick(2) = 1d0-in_suolpick1
!-- group line disparities (strengths):
    gas_ldisp1 = in_ldisp1
    gas_ldisp2 = in_ldisp2
!-- external analytic source input:
    gas_srctype = in_srctype
    gas_theav = in_theav
    gas_nheav = in_nheav
    gas_srcmax = in_srcmax


!-- secondary
    allocate(dd_temp(dd_ncell)) !(dd_ncell)
    allocate(dd_ur(dd_ncell))
    allocate(dd_rho(dd_ncell))
    allocate(dd_bcoef(dd_ncell))
    allocate(dd_nisource(dd_ncell))
    allocate(dd_vol(dd_ncell))        !gcell volume [cm^3]
    allocate(dd_mass(dd_ncell))       !gcell mass
    allocate(dd_natom(dd_ncell))      !gcell # atoms
    allocate(dd_natom1fr(-2:gas_nelem,dd_ncell))
    allocate(dd_natom0fr(-2:2,dd_ncell))
    allocate(dd_nelec(dd_ncell))  !gcell # electrons per atom
    allocate(dd_matsrc(dd_ncell))
    dd_natom1fr = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
    dd_natom0fr = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
    dd_nelec = 1d0  !gcell # electrons per atom 
    dd_matsrc = 0d0  !-- material energy (temperature) source (may be manufactured)
    allocate(dd_emitprob(gas_ng,ncell))
    allocate(dd_cap(gas_ng,ncell))
!   allocate(dd_opacleak(6,ncell))
    allocate(dd_sig(ncell))
    allocate(dd_capgam(ncell))
    allocate(dd_siggrey(ncell))
    allocate(dd_fcoef(ncell))
!
    allocate(dd_eraddens(dd_ncell))
    allocate(dd_edep(ncell))
!
    allocate(dd_emit(ncell))
    allocate(dd_emitex(ncell))
    allocate(dd_evolinit(ncell))
!
!-- output
    if(ltalk) then
     n = dd_ncell*(11 + 5 + gas_nelem+3)/1024 !kB
     write(6,*) 'ALLOC gas    :',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif !ltalk
  end subroutine gasgrid_init


end module gasgridmod
