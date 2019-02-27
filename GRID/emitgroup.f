*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      pure function emitgroup(r,ic) result(ig)
c     --------------------------------------
      use miscmod
      use groupmod
      use gridmod, dummy=>emitgroup
      use physconstmod
      implicit none
      integer :: ig
      real*8,intent(in) :: r
      integer,intent(in) :: ic
************************************************************************
* Determine the group in which to emit a particle.
************************************************************************
      real*8 :: r1
      integer :: l,iep,nepg,igp1
      real*8 :: specval(grd_nepg)
      real*8 :: emitprob
c
c-- search unnormalized cumulative emission probability values
      r1 = r*grd_capgrey(ic)
      iep = binsrch(r1,grd_emitprob(:,ic),grd_nep,.true.)
      ig = iep*grd_nepg + 1
      igp1 = min(ig + grd_nepg - 1, grp_ng)
      nepg = igp1 - ig + 1
      call specintv(grd_tempinv(ic),nepg,specval,offset=ig)
c
c-- start value
      if(iep==0) then
       emitprob = 0d0
      else
       emitprob = grd_emitprob(iep,ic)
      endif
c
c-- step up until target r1 is reached
      l = 0
      do ig=ig,igp1-1
       l = l + 1
       emitprob = emitprob + specval(l)*grd_cap(ig,ic)
       if(emitprob>r1) exit
      enddo
!     if(ig>grp_ng) stop 'transport1: ig not valid'
c
      end function emitgroup
c vim: fdm=marker
