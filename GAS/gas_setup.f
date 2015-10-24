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
      real*8 :: mass0fr(-2*gas_nchain:gas_nelem,gas_ncell)
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
      elseif(in_consttemp/=0d0) then
       gas_temp = in_consttemp
      elseif(str_ltemp) then
       gas_temp = str_tempdd
      else
       stop 'gas_setup: initial gas temp not specified'
      endif
c
c
c-- used in fleck_factor
      if(in_tempradinit>=0d0) then
        gas_eraddens = pc_acoef*in_tempradinit**4
      else
        gas_eraddens = pc_acoef*gas_temp**4
      endif
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
     &  call massfr2natomfr(mass0fr)
c
c-- electron fraction at t=0
      if(str_lye) then
c-- from input.str
       gas_ye0 = str_yedd
      else
c-- as calculated in massfr2natomfr
       gas_ye0 = gas_ye
      endif
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
      real*8,intent(inout) :: mass0fr(-2*gas_nchain:gas_nelem,gas_ncell)
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,l
      real*8 :: help
      logical :: lwarn
c
      lwarn = .true.
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
       help = sum(mass0fr(1:,i))
       mass0fr(:,i) = mass0fr(:,i)/help
c-- warn if renormalization factor is abnormal
       if(abs(help-1d0)>.01 .and. lwarn) then
        write(0,*) 'WARNING: large abundance normalization factor!',help
        write(6,*) 'WARNING: large abundance normalization factor!',help
        lwarn = .false. !warn once
       endif
c
c-- partial mass
       gas_natom1fr(:,i) = mass0fr(:,i)*gas_mass(i)
c
c-- take out radioactive part
       gas_natom1fr(28,i) = gas_natom1fr(28,i) -
     &   gas_natom1fr(gas_ini56,i)
       gas_natom1fr(27,i) = gas_natom1fr(27,i) -
     &   gas_natom1fr(gas_ico56,i)
       gas_natom1fr(26,i) = gas_natom1fr(26,i) -
     &   gas_natom1fr(gas_ife52,i)
       gas_natom1fr(25,i) = gas_natom1fr(25,i) -
     &   gas_natom1fr(gas_imn52,i)
       gas_natom1fr(24,i) = gas_natom1fr(24,i) -
     &   gas_natom1fr(gas_icr48,i)
       gas_natom1fr(23,i) = gas_natom1fr(23,i) -
     &   gas_natom1fr(gas_iv48,i)
c
c-- verify that 
       if(any(gas_natom1fr(23:28,i)<0d0)) then
        write(0,*) gas_natom1fr(:,i)
        stop 'massfr2natomfr: Check input structure: ' //
     &    'abund(unstable) > abund(stable+unstable)'
       endif
c
c-- convert to natoms
       do l=1,gas_nelem
        gas_natom1fr(l,i) = gas_natom1fr(l,i)/(elem_data(l)%m*pc_amu)
       enddo !j
c-- ni56/co56
       help = elem_data(26)%m*pc_amu
       gas_natom1fr(gas_ini56,i) = gas_natom1fr(gas_ini56,i)/help
       gas_natom1fr(gas_ico56,i) = gas_natom1fr(gas_ico56,i)/help
c-- fe52/mn52
       help = elem_data(24)%m*pc_amu
       gas_natom1fr(gas_ife52,i) = gas_natom1fr(gas_ife52,i)/help
       gas_natom1fr(gas_imn52,i) = gas_natom1fr(gas_imn52,i)/help
c-- fe52/mn52
       help = elem_data(22)%m*pc_amu
       gas_natom1fr(gas_icr48,i) = gas_natom1fr(gas_icr48,i)/help
       gas_natom1fr(gas_iv48,i) = gas_natom1fr(gas_iv48,i)/help
c
c-- calculate Ye
       gas_ye(i) = gas_natom1fr(0,i)*.5d0 !container with unknown elements
c-- stable abundances
       do l=1,gas_nelem
        gas_ye(i) = gas_ye(i) + gas_natom1fr(l,i)*l/elem_data(l)%m
       enddo
c-- unstable abundances
       if(gas_nchain/=3) stop 'massfr2natomfr: gas_nchain updated'
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ini56,i)*(28/56d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ico56,i)*(27/56d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_ife52,i)*(26/52d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_imn52,i)*(25/52d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_icr48,i)*(24/48d0)
       gas_ye(i) = gas_ye(i) + gas_natom1fr(gas_iv48,i)*(23/48d0)
c-- normalize
       gas_ye(i) = gas_ye(i)/sum(gas_natom1fr(:,i))
c
c-- store initial numbers
c-- ni/co/fe
       gas_natom0fr(-2,i,1) = gas_natom1fr(gas_ini56,i)!unstable
       gas_natom0fr(-1,i,1) = gas_natom1fr(gas_ico56,i)!unstable
       gas_natom0fr(0:2,i,1) = gas_natom1fr(26:28,i)!stable
c-- fe/mn/cr
       gas_natom0fr(-2,i,2) = gas_natom1fr(gas_ife52,i)!unstable
       gas_natom0fr(-1,i,2) = gas_natom1fr(gas_imn52,i)!unstable
       gas_natom0fr(0:2,i,2) = gas_natom1fr(24:26,i)!stable
c-- cr/v/ti
       gas_natom0fr(-2,i,3) = gas_natom1fr(gas_icr48,i)!unstable
       gas_natom0fr(-1,i,3) = gas_natom1fr(gas_iv48,i)!unstable
       gas_natom0fr(0:2,i,3) = gas_natom1fr(22:24,i)!stable
c
c-- add radioactive part to stable again
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
       gas_natom0fr(:,i,:) = gas_natom0fr(:,i,:)/gas_natom(i)
c
      enddo !i
c!}}}
      end subroutine massfr2natomfr
