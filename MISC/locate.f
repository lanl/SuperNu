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
      pure integer function locate(xarray,n,x)
c     -------------------------------
      implicit none
      integer,intent(in) :: n
      real*8,intent(in) :: xarray(n),x
****************************************************************
* search in an ordered table.
* from numerical recipes
****************************************************************
      integer jl,jm,ju
c
      jl = 0
      ju = n + 1
      do
       jm = (ju + jl)/2
       if((xarray(n)>=xarray(1)) .eqv. (x>=xarray(jm))) then
        jl=jm
       else
        ju=jm
       endif
       if(ju-jl .le. 1) exit
      enddo
c--
      locate = jl
      end function locate
c vim: fdm=marker
