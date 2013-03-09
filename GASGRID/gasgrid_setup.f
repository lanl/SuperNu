      subroutine gasgrid_setup
c     --------------------------------------
      use physconstmod
      use inputparmod
      use gasgridmod
      use miscmod, only:warn
      IMPLICIT NONE
************************************************************************
* Initialize the gas grid, the part that is constant with time and
* temperature. The part that changes is done in gas_grid_update.
************************************************************************
      integer :: i
      REAL*8 :: help
c
c--
      write(6,*)
      write(6,*) 'setup gas grid:'     
      write(6,*) '==========================='
c
c-- init temp structure
      gas_vals2%temp = 1d4  !initial guess, may be overwritten by read_temp_str
c
c-- convert mass fractions to # atoms
      call massfr2natomfr
!c
!c-- output
!      write(6,*) 'mass fractions'
!      write(6,'(1p,33i12)') (i,i=-2,30)
!      write(6,'(1p,33e12.4)') (gas_vals2(i)%mass0fr,i=1,gas_nr)
!      write(6,*) 'number fractions'
!      write(6,'(1p,33i12)') (i,i=-2,30)
!      write(6,'(1p,33e12.4)') (gas_vals2(i)%natom1fr,i=1,gas_nr)
c
c-- gas wavelength grid
      forall(i=1:gas_ng) gas_wl(i) =
     &  in_wlmin*(in_wlmax/dble(in_wlmin))**((i-1d0)/(gas_ng-1d0))
      gas_dwl = pc_ang*gas_wl*log(in_wlmax/dble(in_wlmin)) /
     &  (gas_ng-1d0) !wl grid bin width
c-- sanity test
      help = sum(gas_dwl)/pc_ang
      if(abs(help/(in_wlmax-in_wlmin) - 1d0) .gt. 1d-3) then
       call warn('gasgrid_setup','ggrid_dwl not accurate')
       write(6,*) help,in_wlmax-in_wlmin
      endif
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
      IMPLICIT NONE
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,j
      REAL*8 :: help
c
      do i=1,gas_nr
c
c-- renormalize (the container fraction (unused elements) is taken out)
       gas_vals2(i)%mass0fr(:) = gas_vals2(i)%mass0fr(:)/
     &   sum(gas_vals2(i)%mass0fr(1:))
c
c-- partial mass
       gas_vals2(i)%natom1fr = gas_vals2(i)%mass0fr*gas_vals2(i)%mass
c-- only stable nickel and cobalt
       gas_vals2(i)%natom1fr(28) = gas_vals2(i)%natom1fr(28) -
     &   gas_vals2(i)%natom1fr(gas_ini56)
       gas_vals2(i)%natom1fr(27) = gas_vals2(i)%natom1fr(27) -
     &   gas_vals2(i)%natom1fr(gas_ico56)
c
c-- convert to natoms
       do j=1,gas_nelem
        gas_vals2(i)%natom1fr(j) = gas_vals2(i)%natom1fr(j)/
     &    (elem_data(j)%m*pc_amu) 
       enddo !j
c-- special care for ni56 and co56
       help = elem_data(26)%m*pc_amu
!      help = elem_data(28)%m*pc_amu !phoenix compatible
       gas_vals2(i)%natom1fr(gas_ini56) =
     &   gas_vals2(i)%natom1fr(gas_ini56)/help
!      help = elem_data(27)%m*pc_amu !phoenix compatible
       gas_vals2(i)%natom1fr(gas_ico56) =
     &   gas_vals2(i)%natom1fr(gas_ico56)/help
c-- store initial fe/co/ni
       gas_vals2(i)%natom0fr(-2:-1) = gas_vals2(i)%natom1fr(-2:-1) !unstable
       gas_vals2(i)%natom0fr(0:2) = gas_vals2(i)%natom1fr(26:28) !stable
c-- add unstable to stable again
       gas_vals2(i)%natom1fr(28) = gas_vals2(i)%natom1fr(28) +
     &   gas_vals2(i)%natom1fr(gas_ini56)
       gas_vals2(i)%natom1fr(27) = gas_vals2(i)%natom1fr(27) +
     &   gas_vals2(i)%natom1fr(gas_ico56)
c
c-- total natom
       gas_vals2(i)%natom = sum(gas_vals2(i)%natom1fr(1:))
c
c-- convert natoms to natom fractions
       gas_vals2(i)%natom1fr = gas_vals2(i)%natom1fr/gas_vals2(i)%natom
       gas_vals2(i)%natom0fr = gas_vals2(i)%natom0fr/gas_vals2(i)%natom
c
      enddo !i
c
      end subroutine massfr2natomfr
