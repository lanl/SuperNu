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
      pure function binsrch(x,arr,ng,widerange)
      integer :: binsrch
      integer,intent(in) :: ng
      real*8,intent(in) :: x
      real*8,intent(in) :: arr(ng)
      logical,intent(in) :: widerange
      end function binsrch
c
      end interface
c
      end module miscmod
c vim: fdm=marker
