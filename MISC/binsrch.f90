function binsrch(x,arr,n,widerange)
  implicit none
  integer :: binsrch

  integer,intent(in) :: n
  real*8,intent(in) :: x
  real*8,intent(in) :: arr(n) !array
  logical,intent(in) :: widerange
!---------------------------------------------------
! Binary search of a real*8 array (arr) of size n.
! Finds index of interval containing x.
! Returns an integer between 1 and n-1, inclusive.
!---------------------------------------------------
  integer :: imin, imax, imid

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
