function binsrch(x,arr,ng,widerange)
  implicit none
  integer :: binsrch

  integer,intent(in) :: ng
  real*8,intent(in) :: x
  real*8,intent(in) :: arr(ng) !array
  logical,intent(in) :: widerange
!---------------------------------------------------
! Binary search of a real*8 array (arr) of size ng.
! Finds index of interval containing x.
! Returns an integer between 1 and ng-1, inclusive.
!---------------------------------------------------
  integer :: imin, imax, imid

  imin = 1
  imax = ng
  imid = (ng+1)/2

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
        elseif(x>arr(ng)) then
           imin = ng
           imid = ng
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
     imid = min(imid,ng)
     imid = max(imid,1)
  endif

  binsrch = imid

end function binsrch
