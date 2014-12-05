      pure function specint3(xarr,n,mode) result(ss)
c     -----------------------------------------
      implicit none
      real*8,intent(in) :: xarr(n)
      integer,intent(in) :: n,mode
      real*8 :: ss(n-1)
************************************************************************
* Integrate normalized Planck spectrum using Newton-Cotes formulae of
* different degrees.
************************************************************************
      real*8,parameter :: one6th=1d0/6d0
      real*8,parameter :: one90th=1d0/90d0
      real*8 :: farr(n)
      real*8 :: f2arr(n-1)
c
c-- edge values are always used
      farr = f3(xarr)
c
      select case(mode)
c-- linear
      case(1)
       ss = .5d0*abs(xarr(2:)-xarr(:n-1)) * (farr(2:)+farr(:n-1))
c
c-- quadratic
      case(2)
c-- midpoints
       f2arr = f3(.5d0*(xarr(2:) + xarr(:n-1)))
c-- simpson's rule
       ss = one6th*abs(xarr(2:)-xarr(:n-1)) *
     &  (farr(2:) + farr(:n-1) + 4d0*f2arr)
c
c-- cubic
      case(4)
c-- quarter points
       f2arr = 12d0*f3(.5d0*(xarr(2:) + xarr(:n-1))) + 32d0*(
     &   f3(.25d0*xarr(2:) + .75d0*xarr(:n-1)) +
     &   f3(.75d0*xarr(2:) + .25d0*xarr(:n-1)))
c-- Boole's rule
       ss = one90th*abs(xarr(2:)-xarr(:n-1)) *
     &   (7d0*(farr(2:) + farr(:n-1)) + f2arr)
c
c-- invalid
      case default
       ss = 0d0
      endselect
c
      contains
c
      elemental real*8 function f3(x)
      real*8,intent(in) :: x
      f3 = x**3/(exp(x) - 1d0)
      end function
c
      end function specint3
