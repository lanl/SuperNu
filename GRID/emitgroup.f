      function emitgroup(i,j,k,r) result(ig)
c     --------------------------------------
      use miscmod
      use groupmod
      use gridmod
      use physconstmod
      implicit none
      integer :: ig
      integer,intent(in) :: i,j,k
      real*8,intent(in) :: r
************************************************************************
* Determine the group in which to emit a particle.
************************************************************************
      real*8 :: r1
      integer :: l,iep,nepg,igp1
      integer,external :: binsrch
      real*8 :: specval(grd_nepg)
      real*8 :: emitprob
c
c-- search unnormalized cumulative emission probability values
      r1 = r*grd_capgrey(i,j,k)
      iep = binsrch(r1,grd_emitprob2(:,i,j,k),grd_nep)
      ig = iep*grd_nepg + 1
      igp1 = min(ig + grd_nepg - 1, grp_ng)
      nepg = igp1 - ig + 1
      specval(:nepg) = specintvp(1d0/grd_temp(i,j,k),ig,igp1)
c
c-- start value
      if(iep==0) then
       emitprob = 0d0
      else
       emitprob = grd_emitprob2(iep,i,j,k)
      endif
c
c-- step up until target r1 is reached
      l = 0
      do ig=ig,igp1-1
       l = l + 1
       emitprob = emitprob + specval(l)*grd_cap(ig,i,j,k)
       if(emitprob>r1) exit
      enddo
      if(ig>grp_ng) stop 'transport1: ig not valid'
c
      end function emitgroup
