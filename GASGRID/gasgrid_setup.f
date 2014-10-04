      subroutine gasgrid_setup
c     ------------------------
      use inputstrmod
      use physconstmod
      use inputparmod
      use timestepmod
      use gasgridmod
      use manufacmod
      use miscmod, only:warn
      use profiledatamod
      implicit none
************************************************************************
* Initialize the gas grid, the part that is constant with time and
* temperature. The part that changes is done in gas_grid_update.
************************************************************************
      integer :: l,ll,i
      real*8 :: help
c
c--
      write(6,*)
      if(gas_isvelocity) then
       write(6,*) 'setup velocity grid:'
       write(6,*) '==========================='
      else
       write(6,*) 'setup spatial grid:'
       write(6,*) '==========================='
      endif
c
c----
c-- agnostic grid setup (rev. 200) ----------------------------------
      gas_xarr = str_velleft  !TODO: str_velleft -> 3D
      gas_dxarr(:) = gas_xarr(2:)-gas_xarr(:gas_nx)
c--------------------------------------------------------------------
c
c-- agnostic mass setup (rev. 200) ----------------------------------
      if(gas_ny>1) stop 'gg_setup: str_mass 1D'
      gas_vals2(:,1,1)%mass = str_mass
c--------------------------------------------------------------------
c----
c
c-- scale gas cell volumes to unit sphere depending on expanding or static
      if(gas_isvelocity) then
       help = gas_velout
      else
       if(gas_ny>1) stop "gg_setup: help: no 2D"
       help = gas_lx
      endif
c
c-- volume of unit-radius sphere shells
!TODO: multiple geometries
      if(gas_ny>1) stop 'gasgrid_setup: 1D volume only'
      forall(i=1:gas_nx) gas_vals2(i,1,1)%volr = 
     &  pc_pi4/3d0*(gas_xarr(i+1)**3 - gas_xarr(i)**3)/help**3  !volume in outer radius units
c
c
c--
      write(6,*)
      write(6,*) 'setup gas grid:'
      write(6,*) '==========================='
c
c-- temperature
      if(in_srctype=='manu') then!{{{
       call init_manutemp
      elseif(in_consttemp==0d0) then
       call read_restart_file
      elseif(in_consttemp<0d0) then
       if(gas_ny>1 .or. gas_nz>1) stop 'gg_setup: trad_prof no 2D/3D'
       gas_temp(:,1,1) = trad_profile(tsp_t)
       write(6,*) 'Trad applied to initial gas temp'
      else
       gas_temp = in_consttemp
      endif!}}}
c
c
c-- temp and ur
      gas_vals2%ur = pc_acoef*gas_temp**4 !initial guess, may be overwritten by read_temp_str
c
c-- adopt partial masses from input file
      if(.not.in_noreadstruct) then
       if(.not.allocated(str_massfr)) stop 'no input.str read'
       if(gas_ny>1) stop 'gg_setup: str_massfr: no 2D'
       do l=1,str_nabund
        ll = str_iabund(l)
        if(ll>gas_nelem) ll = 0 !divert to container
        gas_vals2(:,1,1)%mass0fr(ll) = str_massfr(l,:)
       enddo
      elseif(.not.in_novolsrc) then
       if(gas_ny>1) stop 'gg_setup: str_massfr: no 2D'
        gas_vals2%mass0fr(28) = 1d0 !stable+unstable Ni abundance
        gas_vals2(1:nint(4d0*gas_nx/5d0),1,1)%mass0fr(-1) = 1d0 !Ni56 core
      else
       stop 'gg_setup: no input.str and in_novolsrc=true!'
      endif
c
c-- convert mass fractions to # atoms
      call massfr2natomfr
c
c-- output
C$$$      write(6,*) 'mass fractions'
C$$$      write(6,'(1p,33i12)') (l,l=-2,30)
C$$$      write(6,'(1p,33e12.4)') (gas_vals2(l)%mass0fr,l=1,gas_nx)
C$$$      write(6,*) 'number fractions'
C$$$      write(6,'(1p,33i12)') (l,l=-2,30)
C$$$      write(6,'(1p,33e12.4)') (gas_vals2(l)%natom1fr,l=1,gas_nx)
c
      end subroutine gasgrid_setup
c
c
c
      subroutine massfr2natomfr
c     -------------------------
      use physconstmod
      use elemdatamod, only:elem_data
      use gasgridmod
      implicit none
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,j,k,l
      real*8 :: help
c
      do k=1,gas_nz
      do j=1,gas_ny
      do i=1,gas_nx
c!{{{
c-- sanity test
       if(all(gas_vals2(i,j,k)%mass0fr(1:)==0d0)) stop
     &    'massfr2natomfr: all mass fractions zero'
       if(any(gas_vals2(i,j,k)%mass0fr(1:)<0d0)) stop
     &    'massfr2natomfr: negative mass fractions'
c
c-- renormalize (the container fraction (unused elements) is taken out)
       gas_vals2(i,j,k)%mass0fr(:) = gas_vals2(i,j,k)%mass0fr(:)/
     &   sum(gas_vals2(i,j,k)%mass0fr(1:))
c
c-- partial mass
       gas_vals2(i,j,k)%natom1fr = gas_vals2(i,j,k)%mass0fr*
     &   gas_vals2(i,j,k)%mass
c-- only stable nickel and cobalt
       gas_vals2(i,j,k)%natom1fr(28) = gas_vals2(i,j,k)%natom1fr(28) -
     &   gas_vals2(i,j,k)%natom1fr(gas_ini56)
       gas_vals2(i,j,k)%natom1fr(27) = gas_vals2(i,j,k)%natom1fr(27) -
     &   gas_vals2(i,j,k)%natom1fr(gas_ico56)
c
c-- convert to natoms
       do l=1,gas_nelem
        gas_vals2(i,j,k)%natom1fr(l) = gas_vals2(i,j,k)%natom1fr(l)/
     &    (elem_data(l)%m*pc_amu)
       enddo !j
c-- special care for ni56 and co56
!      help = elem_data(26)%m*pc_amu
       help = elem_data(28)%m*pc_amu !phoenix compatible
       gas_vals2(i,j,k)%natom1fr(gas_ini56) =
     &   gas_vals2(i,j,k)%natom1fr(gas_ini56)/help
       help = elem_data(27)%m*pc_amu !phoenix compatible
       gas_vals2(i,j,k)%natom1fr(gas_ico56) =
     &   gas_vals2(i,j,k)%natom1fr(gas_ico56)/help
c-- store initial fe/co/ni
       gas_vals2(i,j,k)%natom0fr(-2:-1)=gas_vals2(i,j,k)%natom1fr(-2:-1)!unstable
       gas_vals2(i,j,k)%natom0fr(0:2) = gas_vals2(i,j,k)%natom1fr(26:28)!stable
c-- add unstable to stable again
       gas_vals2(i,j,k)%natom1fr(28) = gas_vals2(i,j,k)%natom1fr(28) +
     &   gas_vals2(i,j,k)%natom1fr(gas_ini56)
       gas_vals2(i,j,k)%natom1fr(27) = gas_vals2(i,j,k)%natom1fr(27) +
     &   gas_vals2(i,j,k)%natom1fr(gas_ico56)
c
c-- total natom
       gas_vals2(i,j,k)%natom = sum(gas_vals2(i,j,k)%natom1fr(1:))
c
c-- convert natoms to natom fractions
       gas_vals2(i,j,k)%natom1fr = gas_vals2(i,j,k)%natom1fr/
     &   gas_vals2(i,j,k)%natom
       gas_vals2(i,j,k)%natom0fr = gas_vals2(i,j,k)%natom0fr/
     &   gas_vals2(i,j,k)%natom
c!}}}
      enddo !i
      enddo !j
      enddo !k
c
      end subroutine massfr2natomfr
