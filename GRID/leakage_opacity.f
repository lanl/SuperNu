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
      subroutine leakage_opacity
c     ----------------------------------------
      use gridmod
      use timingmod
      implicit none
************************************************************************
* wrapper
************************************************************************
      real*8 :: t0,t1
c
      t0 = t_time()
c
      select case(grd_igeom)
      case(1)
         call leakage_opacity1
      case(2)
         call leakage_opacity2
      case(3)
         call leakage_opacity3
      case(11)
         call leakage_opacity11
      case default
         stop 'leakage_opacity: invalid igeom'
      endselect
c
      t1 = t_time()
      call timereg(t_opacleak,t1-t0) 
c
      end subroutine leakage_opacity
c vim: fdm=marker
