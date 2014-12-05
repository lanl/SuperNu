      module miscmod
c     --------------
      implicit none
************************************************************************
* Avoid explicit interfaces for these subroutines.
************************************************************************
      interface
c
      function memusg() result(mbsize)
      integer :: mbsize(2)
      end function memusg
c
      subroutine warn(source,mesg,sunt)
      character(*),intent(in) :: source
      character(*),intent(in) :: mesg
      character(*),intent(in),optional :: sunt
      end subroutine warn
c
      function lcase(input_string) result(output_string)
      character(*),intent(in) :: input_string
      character(len(input_string)) :: output_string
      end function lcase
c
      elemental function specint(t1,t2,n,m) result(ss)
      real*8 :: ss
      integer,intent(in) :: n
      real*8,intent(in) :: t1,t2
      integer,intent(in),optional :: m
      end function specint
c
      pure function specint3(xarr,n,mode) result(ss)
c     -----------------------------------------
      real*8 :: ss(n-1)
      real*8,intent(in) :: xarr(n)
      integer,intent(in) :: n,mode
      end function specint3

      end interface
c
      end module miscmod
