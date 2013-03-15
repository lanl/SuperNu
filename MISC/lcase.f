      function lcase(input_string) result(output_string)
c     --------------------------------------------------
      implicit none
      character(*),intent(in) :: input_string
      character(len(input_string)) :: output_string
************************************************************************
* return lowercase string
************************************************************************
c-- local variables
      integer :: i,n
      character(*),parameter :: lower_case =
     &  'abcdefghijklmnopqrstuvwxyz'
      character(*),parameter :: upper_case =
     &  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
c
      output_string = input_string
      do i=1,len(output_string)
       n = index(upper_case,output_string(i:i))
       if(n /= 0) output_string(i:i) = lower_case(n:n)
      enddo
      end function lcase
