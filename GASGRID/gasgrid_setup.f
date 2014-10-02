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
      integer :: i,j,ir
      real*8 :: help
      real*8 :: help2
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
      gas_rarr = str_velleft
      gas_drarr(:) = gas_rarr(2:)-gas_rarr(:gas_nr)
c--------------------------------------------------------------------
c
c-- agnostic mass setup (rev. 200) ----------------------------------
      gas_vals2%mass = str_mass
c--------------------------------------------------------------------
c----
c
c-- scale gas cell volumes to unit sphere depending on expanding or static
      if(gas_isvelocity) then
       help = gas_velout
      else
       help = gas_lr
      endif
c
c-- volume of unit-radius sphere shells
      forall(ir=1:gas_nr) gas_vals2(ir)%volr = 
     &  pc_pi4/3d0*(gas_rarr(ir+1)**3 - gas_rarr(ir)**3)/help**3  !volume in outer radius units
c
c
c--
      write(6,*)
      write(6,*) 'setup gas grid:'
      write(6,*) '==========================='
c
c-- temperature
      if(in_srctype=='manu') then!{{{
       call init_manutemp(tsp_t)
      elseif(in_consttemp==0d0) then
       call read_restart_file
      elseif(in_consttemp<0d0) then
       gas_temp = trad_profile(tsp_t)
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
       do i=1,str_nabund
        j = str_iabund(i)
        if(j>gas_nelem) j = 0 !divert to container
        gas_vals2(:)%mass0fr(j) = str_massfr(i,:)
       enddo
      else
       stop 'gg_setup: no input.str and no solidni56!'
      endif
c
c-- convert mass fractions to # atoms
      call massfr2natomfr
c
c-- output
C$$$      write(6,*) 'mass fractions'
C$$$      write(6,'(1p,33i12)') (i,i=-2,30)
C$$$      write(6,'(1p,33e12.4)') (gas_vals2(i)%mass0fr,i=1,gas_nr)
C$$$      write(6,*) 'number fractions'
C$$$      write(6,'(1p,33i12)') (i,i=-2,30)
C$$$      write(6,'(1p,33e12.4)') (gas_vals2(i)%natom1fr,i=1,gas_nr)
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
      integer :: i,j
      real*8 :: help
c
      do i=1,gas_nr
c!{{{
c-- sanity test
       if(all(gas_vals2(i)%mass0fr(1:)==0d0)) stop 'massfr2natomfr: '//
     &   'all mass fractions zero'
       if(any(gas_vals2(i)%mass0fr(1:)<0d0)) stop 'massfr2natomfr: '//
     &   'negative mass fractions'
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
!      help = elem_data(26)%m*pc_amu
       help = elem_data(28)%m*pc_amu !phoenix compatible
       gas_vals2(i)%natom1fr(gas_ini56) =
     &   gas_vals2(i)%natom1fr(gas_ini56)/help
       help = elem_data(27)%m*pc_amu !phoenix compatible
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
c!}}}
      enddo !i
c
      end subroutine massfr2natomfr
