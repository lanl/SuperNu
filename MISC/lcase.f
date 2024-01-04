* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
c vim: fdm=marker
