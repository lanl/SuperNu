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
      subroutine opacity_planckmean
c     -----------------------------
      use inputparmod
      use groupmod
      use physconstmod
      use gasmod
      use miscmod
      implicit none
************************************************************************
* Convert opacity into grey opacities
************************************************************************
      integer :: i
!     integer :: ig
!     real*8 :: x1,x2
!     real*8 :: help,hlparr(grp_ng+1)
      real*8 :: specarr(grp_ng)
c
c-- Planck opacity
      do i=1,gas_ncell
       call specintv(1/gas_temp(i),grp_ng,specarr)
       gas_capgrey(i) = sum(gas_cap(:,i)*specarr)
      enddo !i
c-- compute Planck weighted integral of emission opacity
      if (in_doemiss) then
        do i=1,gas_ncell
           call specintv(1/gas_temp(i),grp_ng,specarr)
           gas_capgrey(i) = sum(gas_cap(:,i)*specarr)
        enddo !i
      endif
c
c
cc-- old
c      do i=1,gas_ncell
c       help = pc_h*pc_c/(pc_kb*gas_temp(i))
c       hlparr = help*grp_wlinv
c       gas_capgrey(i) = 15d0/pc_pi**4 *
cc-- use the same specint resolution as in emission_probability!
c     &   sum(gas_cap(:,i)*specint(hlparr(2:),hlparr(:grp_ng),3,10))
c      enddo !i
c
c
cc-- older
c     gas_capgrey = 0d0
c     do i=1,gas_ncell
c      do ig=1,grp_ng
c       x1 = pc_h*pc_c*grp_wlinv(ig + 1)/(pc_kb*gas_temp(i))
c       x2 = pc_h*pc_c*grp_wlinv(ig)/(pc_kb*gas_temp(i))
c       gas_capgrey(i) = gas_capgrey(i) + 15d0/pc_pi**4*
c    &    gas_cap(ig,i)*specint(x1,x2,3)
c      enddo
c     enddo !i
c
      end subroutine opacity_planckmean
c vim: fdm=marker
