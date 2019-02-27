!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure function binsrch(x,arr,n,widerange)
  implicit none
  integer :: binsrch

  integer,intent(in) :: n
  real*8,intent(in) :: x
  real*8,intent(in) :: arr(n) !array
  logical,intent(in) :: widerange
!---------------------------------------------------
! Binary search of a real*8 array (arr) of size n.
! Finds index of interval containing x.
! Returns an integer between 1 and n-1, inclusive, if not widerange.
! Returns an integer between 0 and n, inclusive, if widerange.
!---------------------------------------------------
  integer :: imin, imax, imid

!-- quick return for array of size 2
  if(n==2) then
     if(x>=arr(1).and.x<=arr(n)) then
        binsrch = 1
     elseif(x<arr(1)) then
        binsrch = 0
     elseif(x>arr(n)) then
        binsrch = n
     else
!--  invalid inputs
        binsrch = -1 !invalid
        return
     endif
!-- limit result to within array bounds
     if(.not.widerange) then
        binsrch = 1
     endif
     return
  endif
  
  imin = 1
  imax = n
  imid = (n+1)/2

  do while(imax - imin > 1)
     if(x>=arr(imin).and.x<arr(imid)) then
        imax = imid
        imid = (imax+imin)/2
     elseif(x>=arr(imid).and.x<=arr(imax)) then
        imin = imid
        imid = (imax+imin)/2
     else
        if(x<arr(1)) then
           imin = 0
           imid = 0
           exit
        elseif(x>arr(n)) then
           imin = n
           imid = n
           exit
        else
!--  invalid inputs
           binsrch = -1 !invalid
           return
        endif
     endif
  enddo

  if(imid/=imin) imid = -1 !invalid

!-- limit result to within array bounds
  if(.not.widerange) then
     imid = min(imid,n-1)
     imid = max(imid,1)
  endif

  binsrch = imid

end function binsrch
! vim: fdm=marker
