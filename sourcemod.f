*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module sourcemod
c     ----------------
      implicit none
      integer :: src_ns,src_ninit
      integer :: src_nnew
      integer :: src_nnonth
      integer :: src_nsurf
      integer :: src_ninitnew=0
      integer*8 :: src_nvacantmin=huge(src_ns)
      integer :: src_nflux=0
      integer :: src_nsurftot
c
      integer,allocatable :: src_ivacant(:) !array of vacant particle array locations
      integer*8,allocatable :: src_nvacantall(:) !(nmpi) number of vacancies on each of the ranks
c
      save
c
      contains

      subroutine sourcemod_init(nmpi)
c     ----------------------------------------
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* init particle module
************************************************************************
      allocate(src_nvacantall(nmpi))
      src_nvacantall = 0
      end subroutine sourcemod_init
c
      end module sourcemod
c vim: fdm=marker
