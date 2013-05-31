      subroutine gasgrid_setup
c     ------------------------
      use inputstrmod
      use physconstmod
      use inputparmod
      use timestepmod
      use gasgridmod
      use miscmod, only:warn
      implicit none
************************************************************************
* Initialize the gas grid, the part that is constant with time and
* temperature. The part that changes is done in gas_grid_update.
************************************************************************
      integer :: i,ir
      real*8 :: help, rrcenter, uudd, masv
      real*8 :: help2
      real*8, allocatable :: wlstore(:)
      integer :: ng,ngm,nrm,irr
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
c-- WARNING: rarr first is [0,1], later [0,vout]. Grid with unit radius: scaled up below by gas_lr for static, or gas_velout for velocity grid
c-- inner shell radius
      if(.not.gas_isshell) then
       gas_rarr(1) = 0.0d0    !Initial inner most radius
      else
       if(gas_isvelocity) then
        gas_rarr(1) = gas_v0/(gas_velout-gas_v0)
       else
        gas_rarr(1) = gas_l0/gas_lr
       endif
      endif
c
c-- outer shells
      do ir=1,gas_nr
       gas_drarr(ir) = 1d0/real(gas_nr)
       gas_rarr(ir+1) = gas_rarr(ir)+gas_drarr(ir)
       gas_vals2(ir)%dr3_34pi = gas_rarr(ir+1)**3-gas_rarr(ir)**3
      enddo
c
c-- volume of unit-radius sphere shells
      forall(ir=1:gas_nr) gas_vals2(ir)%volr = 
     &  pc_pi4/3d0*(gas_rarr(ir+1)**3 - gas_rarr(ir)**3)  !volume in outer radius units
c
c-- set grid size to velout or lr depending on expanding or static
      if(gas_isvelocity) then
       help = gas_velout - gas_v0
      else
       help = gas_lr
      endif
      gas_rarr = gas_rarr*help
      gas_drarr = gas_drarr*help
      gas_vals2%dr3_34pi = gas_vals2%dr3_34pi*help**3
c
c
c
c
c-- IMC-DDMC heuristic coefficient for spherical geometry (rev 146)
c==================================================================
      do ir = 1, gas_nr
       help2=gas_rarr(ir+1)**2+gas_rarr(ir)*gas_rarr(ir+1)
       help2=help2+gas_rarr(ir)**2
       help2 = sqrt(1d0/help2)
       gas_curvcent(ir) = sqrt((gas_rarr(ir+1)**2+gas_rarr(ir)**2))
       gas_curvcent(ir) = help2*gas_curvcent(ir)
      enddo
c
c
c--
      write(6,*)
      write(6,*) 'setup gas grid:'
      write(6,*) '==========================='
c
c-- temperature
      if(in_istempflat) then
       do ir=1,gas_nr
        gas_vals2(ir)%tempkev = in_consttempkev
       enddo
      else
       call read_restart_file
       do ir=1,gas_nr
        if(gas_vals2(ir)%tempkev<1d-6) then
         gas_vals2(ir)%tempkev=1d-6
        endif
       enddo
      endif
c-- Ryan W.: temporary override of initial temperature
c-- to Gaussian for manufacture tests
      if(in_srctype=='manu') then
       uudd = 2.5d8
       do ir=1,gas_nr
        rrcenter=(gas_rarr(ir+1)+gas_rarr(ir))/2d0
        gas_vals2(ir)%tempkev = in_templ0*exp(-0.5*(rrcenter/uudd)**2)
        if(gas_vals2(ir)%tempkev<1d-3) then
         gas_vals2(ir)%tempkev=1d-3
        endif
       enddo
      endif
c
c
c-- mass
      if(in_totmass==0d0) then
       if(.not.allocated(str_mass)) error stop 'gg_setup: str_mass allc'
       gas_vals2%mass = str_mass
      elseif(in_srctype=='manu') then
       do ir=1,gas_nr
        masv = in_totmass*gas_vals2(ir)%dr3_34pi
        masv = masv/(gas_rarr(gas_nr+1)**3-gas_rarr(1)**3)
        gas_vals2(ir)%mass = masv
       enddo
      else
       error stop 'gg_setup: in_totmass was not implemented yet'
      endif
c
c-- temp and ur
      gas_vals2%temp = gas_vals2%tempkev * 1e3*pc_ev/pc_kb !initial guess, may be overwritten by read_temp_str
      gas_vals2%ur = pc_acoef*gas_vals2%tempkev**4
c
c-- set flat composition if selected
      if(in_solidni56) then
       gas_vals2%mass0fr(28) = 1d0 !stable+unstable Ni abundance
       gas_vals2(1:nint(4d0*gas_nr/5d0))%mass0fr(-1) = 1d0 !Ni56 core
      else
       if(.not.allocated(str_mass)) error stop 'no input.str read'
c-- no alternative implemented yet
       stop 'gg_setup: in_solidni56==.false. no mass fractions defined'
      endif
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
c
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
