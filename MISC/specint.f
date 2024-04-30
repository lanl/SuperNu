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
      elemental function specint(x1,x2,m,nn) result(ss)
c     -------------------------------------------------
      implicit none
      real*8 :: ss
      integer,intent(in) :: m
      real*8,intent(in) :: x1,x2
      integer,intent(in),optional :: nn
!#########################################
! For m=3, this function integrates normalized
! Planck spectrum from 0 to t
! Generally, m>=2
!#########################################
      integer :: n,i
      real*8 :: h,x
c
      ss = 0d0
      if (x2 == 0.0) return
c
      if(present(nn)) then
        n = nn
      else
        n = 100
      endif
c
      h = (x2-x1)/n
c      
      ss = x1**m/(exp(x1) - 1d0) + x2**m/(exp(x2) - 1d0)
c
      do i=1,n,2
       x = x1 + i*h
       ss = ss + 4d0*x**m/(exp(x) - 1d0)
      enddo
      do i=2,n-1,2
       x = x1 + i*h
       ss = ss + 2d0*x**m/(exp(x) - 1d0)
      enddo
c
      ss = ss*h/3d0
c
      !x = (x2+x1)/2d0
      !dd = (x**m)/(exp(x)-1d0)
      !ss = (x2-x1)*dd
      end function specint
c vim: fdm=marker
