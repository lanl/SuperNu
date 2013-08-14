      module manufacmod
c     -----------------
      implicit none
************************************************************************
* any constants or factors from particular manufactured solutions
* to be used in source
************************************************************************
c
c
c-- max and min values for constant linear profile function
      real*8,parameter :: man_aa11 = 1.371d14*2.997924562d10 !erg/cm^2/s
      real*8,parameter :: man_aa22 = 1.371d12*2.997924562d10 !erg/cm^2/s
c
c-- a uniform temperature value (or possibly nominal temperature value)
!      real*8,parameter :: man_temp0 = 1.1602621d7 !K
      real*8,parameter :: man_temp0 = 1.160237998048407d7 !K
c
c
      contains
c
c
c
      subroutine check_manufacpars
c     ----------------------------
      use miscmod, only:warn
      use inputparmod, only: in_nobbopac,in_nobfopac,in_noffopac,
     &     in_nothmson,in_dentype
      implicit none
************************************************************************
* check relevant input parameters for consistency with manufactured
* solution constraints
************************************************************************
c
      if(in_opacanaltype=='none')
     &     stop 'check_manufacpars: invalid in_opacanaltype'
      if(.not.in_nobbopac)
     &     stop 'check_manufacpars: invalid in_nobbopac'
      if(.not.in_nobfopac)
     &     stop 'check_manufacpars: invalid in_nobfopac'
      if(.not.in_noffopac)
     &     stop 'check_manufacpars: invalid in_noffopac'
      if(.not.in_nothmson)
     &     stop 'check_manufacpars: invalid in_nothmson'
      if(.not.in_dentype=='mass')
     &     call warn('check_manufacpars','in_dentype/=unif')
c
      end subroutine check_manufacpars
c
c
      subroutine generate_manuradsrc
c     ------------------------------
      use physconstmod
      use gasgridmod
      implicit none
************************************************************************
* calculate finite volume manufactured radiation source in ergs/cm^3/s
* with manufactured parameters
************************************************************************
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(gas_isvelocity) then
c         
c-- implement/modify velocity dependent manufactured radiation source
         select case(gas_opacanaltype)
         case('none')
         case('grey')
         case('mono')
         case('pick')
         case('line')
         case default
            stop 'gas_opacanaltype unknown'
         end select         
c
c
      else
c
c-- implement/modify static manufactured radiation source
         stop 'generate_manuradsrc: no static sources'
c
c
      endif
c
      end subroutine generate_manuradsrc
c
c
      subroutine generate_manutempsrc
c     -------------------------------
      use physconstmod
      use gasgridmod
      implicit none
************************************************************************
* calculate finite volume manufactured temperature source in
* ergs/cm^3/s with manufactured parameters
************************************************************************
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(gas_isvelocity) then
c
c-- implement/modify velocity dependent manufactured temperature source
         select case(gas_opacanaltype)
         case('none')
         case('grey')
         case('mono')
         case('pick')
         case('line')
         case default
            stop 'gas_opacanaltype unknown'
         end select
c
c
      else
c
c-- implement/modify static manufactured temperature source
         stop 'generate_manutempsrc: no static sources'
c
c
      endif
c
      end subroutine generate_manutempsrc
c
      end module manufacmod
