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
      elemental function specint(t1,t2,n) result(ss)
      real*8 :: ss
      integer, intent(in) :: n
      real*8, intent(in) :: t1, t2
      end function specint

      end interface
c
      end module miscmod
