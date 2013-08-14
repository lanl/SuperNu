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
      use inputparmod
      use miscmod, only:warn
      implicit none
************************************************************************
* check relevant input parameters for consistency with manufactured
* solution constraints
************************************************************************
      if(in_srttype=='manu') then
         if(in_opacanaltype=='none')
     &        stop 'check_manufacpars: invalid in_opacanaltype'
         if(.not.in_nobbopac) then
     &        stop 'check_manufacpars: invalid in_nobbopac'
         if(.not.in_nobfopac) then
     &        stop 'check_manufacpars: invalid in_nobfopac'
         if(.not.in_noffopac) then
     &        stop 'check_manufacpars: invalid in_noffopac'
         if(.not.in_nothmson) then
     &        stop 'check_manufacpars: invalid in_nothmson'
         if(.not.in_dentype=='mass') then
     &        call warn('check_manufacpars','in_dentype/=unif')
      endif
      end subroutine check_manufacpars
c
c
      subroutine generate_manuradsrc
c     ------------------------------
      implicit none
************************************************************************
* calculate finite volume manufactured radiation source in ergs/cm^3/s
* with manufactured parameters
************************************************************************
      end subroutine generate_manuradsrc
c
c
      subroutine generate_manutempsrc
c     ------------------------------
      implicit none
************************************************************************
* calculate finite volume manufactured temperature source in
* ergs/cm^3/s with manufactured parameters
************************************************************************      
      end subroutine generate_manutempsrc
c
      end module manufacmod
