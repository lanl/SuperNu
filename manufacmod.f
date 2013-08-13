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
      implicit none
************************************************************************
* check relevant input parameters for consistency with manufactured
* solution constraints
************************************************************************
      end subroutine check_manufacpars
c
c
      end module manufacmod
