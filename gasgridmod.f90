module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: gas_nx = 0
  integer :: gas_ny = 0
  integer :: gas_nz = 0
  integer :: gas_ng = 0

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid

  logical :: gas_isvelocity = .false.
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
!
!-- energy conservation check quantities
  real*8 :: gas_eext = 0d0 !time-integrated input energy from external source
  real*8 :: gas_emat = 0d0 !material energy
  real*8 :: gas_erad = 0d0 !census radiation energy
  real*8 :: gas_eleft = 0d0 !left (inward) leaked energy from domain
  real*8 :: gas_eright = 0d0 !right (outward) leaked energy from domain
  real*8 :: gas_etot = 0d0 !total source energy added per time step
  real*8 :: gas_evelo = 0d0 !total energy change to rad field from fluid
  real*8 :: gas_eerror = 0d0 !error in integral problem energy
!-- average quantities used for energy check with MPI
  real*8 :: gas_eextav = 0d0
  real*8 :: gas_eveloav = 0d0
!--
  real*8 :: gas_esurf

  real*8,allocatable :: gas_xarr(:)   !(gas_nx+1), left cell edge values
  real*8,allocatable :: gas_yarr(:)   !(gas_ny+1), left cell edge values
  real*8,allocatable :: gas_zarr(:)   !(gas_nz+1), left cell edge values


!-- Probability of emission in a given zone and group
  real*8,allocatable :: gas_emitprob(:,:,:,:) !(gas_ng,gas_nx,gas_ny,gas_nz)
!-- Line+Cont extinction coeff
  real*8,allocatable :: gas_cap(:,:,:,:) !(gas_ng,gas_nx,gas_ny,gas_nz)
!-- leakage opacities
  real*8,allocatable :: gas_opacleak(:,:,:,:) !(6,gas_nx,gas_ny,gas_nz)


!-- scattering coefficient
  real*8,allocatable :: gas_sig(:,:,:) !(gas_nx,gas_ny,gas_nz)
!-- Gamma ray gray opacity
  real*8,allocatable :: gas_capgam(:,:,:) !(gas_nx,gas_ny,gas_nz)
!-- Planck opacity (gray)
  real*8,allocatable :: gas_siggrey(:,:,:)!(gas_nx,gas_ny,gas_nz)
!-- Fleck factor
  real*8,allocatable :: gas_fcoef(:,:,:)  !(gas_nx,gas_ny,gas_nz)


!-- energy absorbed by material
  real*8,allocatable :: gas_edep(:,:,:)   !(gas_nx,gas_ny,gas_nz)
!-- radiation energy density in tsp_dt
  real*8,allocatable :: gas_eraddens(:,:,:) !(gas_nx,gas_ny,gas_nz)


!-- number of IMC-DDMC method changes per cell per time step
  integer,allocatable :: gas_methodswap(:,:,:) !(gas_nx,gas_ny,gas_nz)
!-- number of census prt_particles per cell
  integer,allocatable :: gas_numcensus(:,:,:) !(gas_nx,gas_ny,gas_nz)

!
!-- packet number and energy distribution
!========================================
  integer,allocatable :: gas_nvol(:,:,:) !(gas_nx,gas_ny,gas_nz) number of thermal source particles generated per cell
  integer,allocatable :: gas_nvolex(:,:,:) !(gas_nx,gas_ny,gas_nz) number of external source particles generated per cell
  integer,allocatable :: gas_nvolinit(:,:,:) !(gas_nx,gas_ny,gas_nz) number of initial (t=tfirst) particles per cell
!  
  real*8,allocatable :: gas_emit(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of fictitious thermal energy emitted per cell in a time step
  real*8,allocatable :: gas_emitex(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of external energy emitted per cell per group in a time step
  real*8,allocatable :: gas_evolinit(:,:,:) !(gas_nx,gas_ny,gas_nz) amount of initial energy per cell per group
!
  real*8,allocatable :: gas_temp(:,:,:) !(gas_nx,gas_ny,gas_nz)

  private read_temp_preset
!
!
!-- DOMAIN DECOMPOSITION
!=======================
  real*8,allocatable :: dd_temp(:,:,:) !(gas_nx,gas_ny,gas_nz)
  real*8,allocatable :: dd_eraddens(:,:,:)
  real*8,allocatable :: dd_ur(:,:,:)
  real*8,allocatable :: dd_rho(:,:,:)
  real*8,allocatable :: dd_bcoef(:,:,:)
  real*8,allocatable :: dd_nisource(:,:,:)
  real*8,allocatable :: dd_vol(:,:,:)        !gcell volume [cm^3]
  real*8,allocatable :: dd_mass(:,:,:)       !gcell mass
  real*8,allocatable :: dd_natom(:,:,:)      !gcell # atoms
  real*8,allocatable :: dd_natom1fr(:,:,:,:) !(-2:gas_nelem,nx,ny,nz)  !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
  real*8,allocatable :: dd_natom0fr(:,:,:,:) !(-2:2,nx,ny,nz)  !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
  real*8,allocatable :: dd_nelec(:,:,:)      !gcell # electrons per atom
!-- mate,allocatablerial energy (temperature) source (may be manufactured), rev>244
  real*8,allocatable :: dd_matsrc(:,:,:)


!-- temperature structure history
  real*8,allocatable :: gas_temppreset(:,:,:,:) !(gas_nx,gas_ny,gas_nz,tim_nt)


  save

  contains


  subroutine gasgrid_init
!------------------------
    use inputparmod
    implicit none

    integer :: n,nx,ny,nz
    logical :: lexist

    gas_nx = in_ndim(1)
    gas_ny = in_ndim(2)
    gas_nz = in_ndim(3)

    gas_isvelocity = in_isvelocity
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

!-- primary
    allocate(gas_xarr(gas_nx+1)) !zone edge x position
    allocate(gas_yarr(gas_ny+1)) !zone edge y position
    allocate(gas_zarr(gas_nz+1)) !zone edge z position

    nx = gas_nx !shortcut
    ny = gas_ny !shortcut
    nz = gas_nz !shortcut
    allocate(gas_edep(nx,ny,nz))
    allocate(gas_siggrey(nx,ny,nz))
    allocate(gas_capgam(nx,ny,nz))
!
!- Ryan W.: using power law to calculate gas_sig (similar to Planck opacity)
    allocate(gas_sig(nx,ny,nz))    !grey scattering opacity
!----------------------------------------------------------------
    allocate(gas_fcoef(nx,ny,nz))
     allocate(gas_emitprob(gas_ng,nx,ny,nz))
    allocate(gas_opacleak(6,nx,ny,nz))
    allocate(gas_eraddens(nx,ny,nz))

!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
     allocate(gas_cap(gas_ng,nx,ny,nz)) !Line+Cont extinction coeff

!--Ryan W: values below were formerly secondary (rev 183)
    allocate(gas_temp(nx,ny,nz))  !cell average temperature
    allocate(gas_nvol(nx,ny,nz))
    allocate(gas_nvolex(nx,ny,nz))
    allocate(gas_nvolinit(nx,ny,nz))
    allocate(gas_emit(nx,ny,nz))
    allocate(gas_emitex(nx,ny,nz))
    allocate(gas_evolinit(nx,ny,nz))

    allocate(gas_methodswap(nx,ny,nz))
    allocate(gas_numcensus(nx,ny,nz))

!-- secondary
    allocate(dd_temp(nx,ny,nz)) !(nx,ny,nz)
    allocate(dd_eraddens(nx,ny,nz))
    allocate(dd_ur(nx,ny,nz))
    allocate(dd_rho(nx,ny,nz))
    allocate(dd_bcoef(nx,ny,nz))
    allocate(dd_nisource(nx,ny,nz))
    allocate(dd_vol(nx,ny,nz))        !gcell volume [cm^3]
    allocate(dd_mass(nx,ny,nz))       !gcell mass
    allocate(dd_natom(nx,ny,nz))      !gcell # atoms
    allocate(dd_natom1fr(-2:gas_nelem,nx,ny,nz))
    allocate(dd_natom0fr(-2:2,nx,ny,nz))
    allocate(dd_nelec(nx,ny,nz))  !gcell # electrons per atom
    allocate(dd_matsrc(nx,ny,nz))
    dd_natom1fr = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
    dd_natom0fr = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
    dd_nelec = 1d0  !gcell # electrons per atom 
    dd_matsrc = 0d0  !-- material energy (temperature) source (may be manufactured)
!
!-- output
    n = nx*ny*nz
    n = int(int(n,8)*8*(21 + 11 + 5 + gas_nelem+3)/1024) !kB
    write(6,*) 'ALLOC gasgrid:',n,"kB",n/1024,"MB",n/1024**2,"GB"
    n = nx*ny*nz
    n = int(((8+8)*int(n,8)*gas_ng)/1024) !kB
    write(6,*) 'ALLOC gas_cap:',n,"kB",n/1024,"MB",n/1024**2,"GB"
!
!-- read preset temperature profiles
    inquire(file='input.temp',exist=lexist)
    if(lexist) call read_temp_preset
  end subroutine gasgrid_init


  subroutine read_temp_preset
!-----------------------------!{{{
    use timestepmod
    integer :: istat
!
    open(4,file='input.temp',status='old',iostat=istat)
    if(istat/=0) stop 'rd_temp_preset: no file: input.temp'
!-- alloc and read
    allocate(gas_temppreset(gas_nx,gas_ny,gas_nz,tsp_nt))
    read(4,*,iostat=istat) gas_temppreset
    if(istat/=0) stop 'rd_temp_preset: file too short: input.temp'
!-- check EOF
    read(4,*,iostat=istat)
    if(istat==0) stop 'rd_temp_preset: file too long: input.temp'
    close(4)
    write(6,*) 'rd_temp_preset: custom temp profiles read successfully'
!}}}
  end subroutine read_temp_preset


end module gasgridmod
