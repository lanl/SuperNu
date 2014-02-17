function binsrch(lamp,wl,ng,ngin)
!---------------------------------------------------
! Binary search of a real*8 array (wl) of size ng.
! Finds index of interval containing lamp.
! Returns an integer between 1 and ng-1, inclusive.
!---------------------------------------------------
  implicit none
  integer :: binsrch

  integer, intent(in) :: ng, ngin
  real*8, intent(in) :: lamp
  real*8, intent(in) :: wl(ng) !array
  !
  integer :: imin, imax, imid

  if(ngin==0) then
!-- initialize binary indexes and key
     imin = 1
     imax = ng
     imid = (ng+1)/2

     do while(imax - imin > 1)
        if(lamp>=wl(imin).and.lamp<wl(imid)) then
           imax = imid
           imid = (imax+imin)/2
        elseif(lamp>=wl(imid).and.lamp<=wl(imax)) then
           imin = imid
           imid = (imax+imin)/2
        else
           if(lamp<wl(1)) then
              imin = 0
              imid = 0
              exit
           elseif(lamp>wl(ng)) then
              imin = ng
              imid = ng
              exit
           else
              write(*,*) ng, wl(imin), wl(imid), lamp
              stop 'binsrch: invalid inputs'
           endif
        endif
     enddo

     if(imid/=imin) stop 'binsrch: no index'
     binsrch = imid
!
  else
!-- logarithmic wavelength groups
     binsrch=floor(1d0+(ng-1)*log(lamp/wl(1))/log(wl(ng)/wl(1)))
  endif

end function binsrch
