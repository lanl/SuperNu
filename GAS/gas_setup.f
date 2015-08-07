      subroutine gas_setup
c     --------------------------
      use inputstrmod
      use physconstmod
      use inputparmod
      use gasmod
      use manufacmod
      use miscmod, only:warn
      implicit none
************************************************************************
* Initialize the gas grid, the part that is constant with time and
* temperature. The part that changes is done in gas_grid_update.
************************************************************************
      integer :: l,i
      real*8 :: mass0fr(-2:gas_nelem,gas_ncell)
c
c-- agnostic mass setup
      gas_mass = str_massdd
c
c-- decompose idcell
      i = 0
      do l=gas_icell1,gas_icell1+gas_ncell-1
       i = i+1
       gas_idcell(i) = str_idcell(l)
      enddo !l
c
c-- temperature
      if(in_srctype=='manu') then
       call init_manutemp
      elseif(in_consttemp==0d0) then
       if(in_ntres>1) then
        call read_restart_file
       elseif(str_ltemp) then
         gas_temp = str_tempdd
       else
         stop 'gas_setup: input.temp_str not avail'
       endif
      else
       gas_temp = in_consttemp
      endif
c
c
c-- used in fleck_factor
      gas_eraddens = pc_acoef*in_tempradinit**4
c
c
c-- temp and ur
      gas_ur = pc_acoef*gas_temp**4 !initial guess, may be overwritten by read_temp_str
c
c-- adopt partial masses from input file
      mass0fr = 0d0
      if(.not.in_noreadstruct) then
       if(.not.allocated(str_massfrdd)) stop 'input.str data not avail'
       do l=1,str_nabund
        i = str_iabund(l)
        if(i>gas_nelem) i = 0 !divert to container
        mass0fr(i,:) = str_massfrdd(l,:)
       enddo
      elseif(.not.in_novolsrc) then
        mass0fr(28,:) = 1d0 !stable+unstable Ni abundance
        mass0fr(-1,:) = 1d0
!      else
!       stop 'gg_setup: no input.str and in_novolsrc=true!'
      endif
c
c-- convert mass fractions to # atoms
      if(.not.in_noreadstruct.or..not.in_novolsrc) 
     &     call massfr2natomfr(mass0fr)
c
      end subroutine gas_setup
c
c
c
      subroutine massfr2natomfr(mass0fr)
c     ----------------------------------!{{{
      use physconstmod
      use elemdatamod, only:elem_data
      use gasmod
      implicit none
      real*8,intent(inout) :: mass0fr(-2:gas_nelem,gas_ncell)
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,l
      real*8 :: help
c
      do i=1,gas_ncell
       if(gas_mass(i)<=0) cycle
c-- sanity test
       if(all(mass0fr(1:,i)==0d0)) stop
     &    'massfr2natomfr: all mass fractions zero'
       if(any(mass0fr(1:,i)<0d0)) stop
     &    'massfr2natomfr: negative mass fractions'
c
c-- renormalize (the container fraction (unused elements) is taken out)
       mass0fr(:,i) = mass0fr(:,i)/sum(mass0fr(1:,i))
c
c-- partial mass
       gas_natom1fr(:,i) = mass0fr(:,i)*gas_mass(i)
c-- only stable nickel and cobalt
       gas_natom1fr(28,i) = gas_natom1fr(28,i) -
     &   gas_natom1fr(gas_ini56,i)
       gas_natom1fr(27,i) = gas_natom1fr(27,i) -
     &   gas_natom1fr(gas_ico56,i)
c
c-- convert to natoms
       do l=1,gas_nelem
        gas_natom1fr(l,i) = gas_natom1fr(l,i)/
     &    (elem_data(l)%m*pc_amu)
       enddo !j
c-- special care for ni56 and co56
!      help = elem_data(26)%m*pc_amu
       help = elem_data(28)%m*pc_amu !phoenix compatible
       gas_natom1fr(gas_ini56,i) =
     &   gas_natom1fr(gas_ini56,i)/help
       help = elem_data(27)%m*pc_amu !phoenix compatible
       gas_natom1fr(gas_ico56,i) =
     &   gas_natom1fr(gas_ico56,i)/help
c-- store initial ni/co/fe
       gas_natom0fr(-2,i,1) = gas_natom1fr(gas_ini56,i)!unstable
       gas_natom0fr(-1,i,1) = gas_natom1fr(gas_ico56,i)!unstable
       gas_natom0fr(0:2,i,1) = gas_natom1fr(26:28,i)!stable
c-- store initial fe/mn/cr
       gas_natom0fr(-2,i,2) = gas_natom1fr(gas_ife52,i)!unstable
       gas_natom0fr(-1,i,2) = gas_natom1fr(gas_imn52,i)!unstable
       gas_natom0fr(0:2,i,2) = gas_natom1fr(24:26,i)!stable
c-- store initial cr/v/ti
       gas_natom0fr(-2,i,2) = gas_natom1fr(gas_icr48,i)!unstable
       gas_natom0fr(-1,i,2) = gas_natom1fr(gas_iv48,i)!unstable
       gas_natom0fr(0:2,i,2) = gas_natom1fr(22:24,i)!stable
c
c-- add unstable to stable again
       gas_natom1fr(28,i) = gas_natom1fr(28,i) +
     &   gas_natom1fr(gas_ini56,i)
       gas_natom1fr(27,i) = gas_natom1fr(27,i) +
     &   gas_natom1fr(gas_ico56,i)
       gas_natom1fr(26,i) = gas_natom1fr(26,i) +
     &   gas_natom1fr(gas_ife52,i)
       gas_natom1fr(25,i) = gas_natom1fr(25,i) +
     &   gas_natom1fr(gas_imn52,i)
       gas_natom1fr(24,i) = gas_natom1fr(24,i) +
     &   gas_natom1fr(gas_icr48,i)
       gas_natom1fr(23,i) = gas_natom1fr(23,i) +
     &   gas_natom1fr(gas_iv48,i)
c
c-- total natom
       gas_natom(i) = sum(gas_natom1fr(1:,i))
c
c-- convert natoms to natom fractions
       gas_natom1fr(:,i) = gas_natom1fr(:,i)/gas_natom(i)
       gas_natom0fr(:,i,1) = gas_natom0fr(:,i,1)/gas_natom(i)
c
      enddo !i
c!}}}
      end subroutine massfr2natomfr
