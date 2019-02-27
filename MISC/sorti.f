*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine sorti(n,arrin,indx)
c     ------------------------------
      implicit none
      integer,intent(in):: n
      real*8,intent(in):: arrin(n)
      integer,intent(inout):: indx(n)
************************************************************************
* index-based sorting routine taken from the routine 'indexx' by
* press etal: 1986, 'numerical recipes', cambridge university press
************************************************************************
      integer :: i,j,l,ir,indxt
      real*8 :: q
      if(n.lt.1) return
      if(n.eq.1) then
       indx(1) = 1
       return
      endif
c
      l = n/2+1
      ir = n
      do
       if(l.gt.1)then
        l = l-1
        indxt = indx(l)
        q = arrin(indxt)
       else
        indxt = indx(ir)
        q = arrin(indxt)
        indx(ir) = indx(1)
        ir = ir-1
        if(ir.eq.1)then
          indx(1) = indxt
          return
        endif
       endif
c
       i = l
       j = l+l
c
       do
        if(j.gt.ir) exit
        if(j.lt.ir)then
         if(arrin(indx(j)).lt.arrin(indx(j+1))) j = j+1
        endif
        if(q.lt.arrin(indx(j)))then
         indx(i) = indx(j)
         i = j
         j = j+j
        else
         j = ir+1
        endif
       enddo
       indx(i) = indxt
      enddo
      end subroutine sorti
c vim: fdm=marker
