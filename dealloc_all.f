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
      subroutine dealloc_all
c     ----------------------
      use mpimod
      use ionsmod
      use bbxsmod
      use gridmod
      use groupmod
      use gasmod
      use particlemod
      use fluxmod
      use sourcemod
      use randommod
      use tbxsmod
      use tbsrcmod
      use timestepmod
      implicit none
************************************************************************
* deallocate all that was used till the end of the program. Any
* allocatable arrays remaining allocated after this had to be dealt
* with earlier.  This helps to catch memory leaks! (drr)
************************************************************************
c-- ionsmod
      call ions_dealloc
      call gas_dealloc
      call grid_dealloc
      call flux_dealloc
      call tbxs_dealloc
      call tbsrc_dealloc
      deallocate(prt_particles,prt_isvacant)
      call mpimod_dealloc
      deallocate(grp_wl,grp_wlinv)
      deallocate(src_nvacantall)
      deallocate(rnd_states)
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f
      deallocate(tsp_tarr)

      end subroutine dealloc_all
c vim: fdm=marker
