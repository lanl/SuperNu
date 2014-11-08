module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

!-- wavelength grid (gridmod has a copy as well)
  integer :: gas_ng=0
  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid

!-- domain decomposed grid variables used to calculate the state of the material (gas)
  integer :: gas_ncell=0
  real*8,allocatable :: gas_temp(:)       !(ncell)
  real*8,allocatable :: gas_eraddens(:)
  real*8,allocatable :: gas_ur(:)
  real*8,allocatable :: gas_rho(:)
  real*8,allocatable :: gas_bcoef(:)
  real*8,allocatable :: gas_nisource(:)
  real*8,allocatable :: gas_vol(:)        !cell volume [cm^3]
  real*8,allocatable :: gas_mass(:)       !cell mass [g]
  real*8,allocatable :: gas_natom(:)      !cell number of atoms
  real*8,allocatable :: gas_nelec(:)      !cell number of electrons per atom
  real*8,allocatable :: gas_natom1fr(:,:) !(-2:gas_nelem,ncell)  !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
  real*8,allocatable :: gas_natom0fr(:,:) !(-2:2,ncell)  !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
!-- mate,allocatablerial energy (temperature) source (may be manufactured), rev>244
  real*8,allocatable :: gas_matsrc(:)
!
  real*8,allocatable :: gas_edep(:)

!== DD copies
!-- Probability of emission in a given zone and group
  real*8,allocatable :: gas_emitprob(:,:) !(gas_ng,ncell)
!-- Line+Cont extinction coeff
  real*8,allocatable :: gas_cap(:,:) !(gas_ng,ncell)
!-- leakage opacities
! real*8,allocatable :: dd_opacleak(:,:) !(6,ncell)
!-- scattering coefficient
  real*8,allocatable :: gas_sig(:) !(ncell)
!-- Gamma ray gray opacity
  real*8,allocatable :: gas_capgam(:) !(ncell)
!-- Planck opacity (gray)
  real*8,allocatable :: gas_siggrey(:)!(ncell)
!-- Fleck factor
  real*8,allocatable :: gas_fcoef(:)  !(ncell)
!  
  real*8,allocatable :: gas_emit(:) !(ncell) amount of fictitious thermal energy emitted per cell in a time step
  real*8,allocatable :: gas_emitex(:) !(ncell) amount of external energy emitted per cell per group in a time step
  real*8,allocatable :: gas_evolinit(:) !(ncell) amount of initial energy per cell per group


!-- temperature structure history (only allocated when used)
  real*8,allocatable :: gas_temppreset(:,:) !(ncell,tim_nt)

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

    gas_ncell = ncell


!-- secondary
    allocate(gas_temp(gas_ncell)) !(gas_ncell)
    allocate(gas_ur(gas_ncell))
    allocate(gas_rho(gas_ncell))
    allocate(gas_bcoef(gas_ncell))
    allocate(gas_nisource(gas_ncell))
    allocate(gas_vol(gas_ncell))        !gcell volume [cm^3]
    allocate(gas_mass(gas_ncell))       !gcell mass
    allocate(gas_natom(gas_ncell))      !gcell # atoms
    allocate(gas_natom1fr(-2:gas_nelem,gas_ncell))
    allocate(gas_natom0fr(-2:2,gas_ncell))
    allocate(gas_nelec(gas_ncell))  !gcell # electrons per atom
    allocate(gas_matsrc(gas_ncell))
    gas_natom1fr = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
    gas_natom0fr = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
    gas_nelec = 1d0  !gcell # electrons per atom 
    gas_matsrc = 0d0  !-- material energy (temperature) source (may be manufactured)
    allocate(gas_emitprob(gas_ng,gas_ncell))
    allocate(gas_cap(gas_ng,gas_ncell))
!   allocate(dd_opacleak(6,gas_ncell))
    allocate(gas_sig(gas_ncell))
    allocate(gas_capgam(gas_ncell))
    allocate(gas_siggrey(gas_ncell))
    allocate(gas_fcoef(gas_ncell))
!
    allocate(gas_eraddens(gas_ncell))
    allocate(gas_edep(gas_ncell))
!
    allocate(gas_emit(gas_ncell))
    allocate(gas_emitex(gas_ncell))
    allocate(gas_evolinit(gas_ncell))
!
!-- output
    if(ltalk) then
     n = gas_ncell*(11 + 5 + gas_nelem+3)/1024 !kB
     write(6,*) 'ALLOC gas    :',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif !ltalk
  end subroutine gasgrid_init


end module gasgridmod
