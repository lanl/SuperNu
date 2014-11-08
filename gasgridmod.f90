module gasgridmod

  implicit none
!***********************************************************************
! gas grid structure
!***********************************************************************
  integer,parameter :: gas_nelem=30
  integer,parameter :: gas_ini56=-1, gas_ico56=-2 !positions in mass0fr and natom1fr arrays

  integer :: dd_ncell=0

  integer :: grd_nx=0
  integer :: grd_ny=0
  integer :: grd_nz=0
  integer :: gas_ng=0

  real*8,allocatable :: gas_wl(:) !(gas_ng) wavelength grid

  logical :: grd_isvelocity = .false.
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
  real*8 :: gas_evelo = 0d0 !total energy change to rad field from fluid
  real*8 :: gas_eerror = 0d0 !error in integral problem energy
!-- average quantities used for energy check with MPI
  real*8 :: gas_eextav = 0d0
  real*8 :: gas_eveloav = 0d0
!--
  real*8 :: gas_esurf

  real*8,allocatable :: grd_xarr(:)   !(grd_nx+1), left cell edge values
  real*8,allocatable :: grd_yarr(:)   !(grd_ny+1), left cell edge values
  real*8,allocatable :: grd_zarr(:)   !(grd_nz+1), left cell edge values


!-- Probability of emission in a given zone and group
  real*8,allocatable :: grd_emitprob(:,:,:,:) !(gas_ng,grd_nx,grd_ny,grd_nz)
!-- Line+Cont extinction coeff
  real*8,allocatable :: grd_cap(:,:,:,:) !(gas_ng,grd_nx,grd_ny,grd_nz)
!-- leakage opacities
  real*8,allocatable :: grd_opacleak(:,:,:,:) !(6,grd_nx,grd_ny,grd_nz)


!-- scattering coefficient
  real*8,allocatable :: grd_sig(:,:,:) !(grd_nx,grd_ny,grd_nz)
!-- Gamma ray gray opacity
  real*8,allocatable :: grd_capgam(:,:,:) !(grd_nx,grd_ny,grd_nz)
!-- Planck opacity (gray)
  real*8,allocatable :: grd_siggrey(:,:,:)!(grd_nx,grd_ny,grd_nz)
!-- Fleck factor
  real*8,allocatable :: grd_fcoef(:,:,:)  !(grd_nx,grd_ny,grd_nz)


!-- energy absorbed by material
  real*8,allocatable :: grd_edep(:,:,:)   !(grd_nx,grd_ny,grd_nz)
!-- radiation energy density in tsp_dt
  real*8,allocatable :: grd_eraddens(:,:,:) !(grd_nx,grd_ny,grd_nz)


!-- number of IMC-DDMC method changes per cell per time step
  integer,allocatable :: grd_methodswap(:,:,:) !(grd_nx,grd_ny,grd_nz)
!-- number of census prt_particles per cell
  integer,allocatable :: grd_numcensus(:,:,:) !(grd_nx,grd_ny,grd_nz)

!
!-- packet number and energy distribution
!========================================
  integer,allocatable :: grd_nvol(:,:,:) !(grd_nx,grd_ny,grd_nz) number of thermal source particles generated per cell
  integer,allocatable :: grd_nvolex(:,:,:) !(grd_nx,grd_ny,grd_nz) number of external source particles generated per cell
  integer,allocatable :: grd_nvolinit(:,:,:) !(grd_nx,grd_ny,grd_nz) number of initial (t=tfirst) particles per cell
!  
  real*8,allocatable :: grd_emit(:,:,:) !(grd_nx,grd_ny,grd_nz) amount of fictitious thermal energy emitted per cell in a time step
  real*8,allocatable :: grd_emitex(:,:,:) !(grd_nx,grd_ny,grd_nz) amount of external energy emitted per cell per group in a time step
  real*8,allocatable :: grd_evolinit(:,:,:) !(grd_nx,grd_ny,grd_nz) amount of initial energy per cell per group
!
  real*8,allocatable :: grd_temp(:,:,:) !(grd_nx,grd_ny,grd_nz)
  real*8,allocatable :: grd_vol(:,:,:) !(grd_nx,grd_ny,grd_nz)

  private read_temp_preset


!-- temperature structure history
  real*8,allocatable :: gas_temppreset(:,:,:,:) !(grd_nx,grd_ny,grd_nz,tim_nt)
  real*8,allocatable :: dd_temppreset(:,:) !(grd_nx,grd_ny,grd_nz,tim_nt)

!
!
!-- DOMAIN DECOMPOSITION
!=======================
  real*8,allocatable :: dd_temp(:) !(grd_nx,grd_ny,grd_nz)
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

    integer :: n,nx,ny,nz
    logical :: lexist

    dd_ncell = ncell

    grd_nx = in_ndim(1)
    grd_ny = in_ndim(2)
    grd_nz = in_ndim(3)

    grd_isvelocity = in_isvelocity
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

!!-- primary
!    allocate(grd_xarr(grd_nx+1)) !zone edge x position
!    allocate(grd_yarr(grd_ny+1)) !zone edge y position
!    allocate(grd_zarr(grd_nz+1)) !zone edge z position

    nx = grd_nx !shortcut
    ny = grd_ny !shortcut
    nz = grd_nz !shortcut
    allocate(grd_edep(nx,ny,nz))
    allocate(grd_siggrey(nx,ny,nz))
    allocate(grd_capgam(nx,ny,nz))
!
!- Ryan W.: using power law to calculate grd_sig (similar to Planck opacity)
    allocate(grd_sig(nx,ny,nz))    !grey scattering opacity
!----------------------------------------------------------------
    allocate(grd_fcoef(nx,ny,nz))
     allocate(grd_emitprob(gas_ng,nx,ny,nz))
    allocate(grd_opacleak(6,nx,ny,nz))
    allocate(grd_eraddens(nx,ny,nz))

!-Ryan W: gas_wl being allocated in gasgrid_setup now--
    !allocate(gas_wl(gas_ng)) !wavelength grid
!------------------------------------------------------
     allocate(grd_cap(gas_ng,nx,ny,nz)) !Line+Cont extinction coeff

!--Ryan W: values below were formerly secondary (rev 183)
    allocate(grd_temp(nx,ny,nz))  !cell average temperature
    allocate(grd_vol(nx,ny,nz))  !cell average temperature

    allocate(grd_nvol(nx,ny,nz))
    allocate(grd_nvolex(nx,ny,nz))
    allocate(grd_nvolinit(nx,ny,nz))
    allocate(grd_emit(nx,ny,nz))
    allocate(grd_emitex(nx,ny,nz))
    allocate(grd_evolinit(nx,ny,nz))

    allocate(grd_methodswap(nx,ny,nz))
    allocate(grd_numcensus(nx,ny,nz))

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
     write(6,*) 'ALLOC dd     :',n,"kB",n/1024,"MB",n/1024**2,"GB"
     n = nx*ny*nz
     n = int((int(n,8)*8*21)/1024) !kB
     write(6,*) 'ALLOC gasgrid:',n,"kB",n/1024,"MB",n/1024**2,"GB"
     n = nx*ny*nz
     n = int(((8+8)*int(n,8)*gas_ng)/1024) !kB
     write(6,*) 'ALLOC grd_cap:',n,"kB",n/1024,"MB",n/1024**2,"GB"
    endif !ltalk
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
    allocate(gas_temppreset(grd_nx,grd_ny,grd_nz,tsp_nt))
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
