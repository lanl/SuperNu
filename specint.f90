function specint(t1,t2,n)

!#########################################
! For n=3, this function integrates normalized
! Planck spectrum from 0 to t
! Generally, n>=2
!#########################################
  real*8 :: specint
  integer, intent(in) :: n
  real*8, intent(in) :: t1, t2

  real*8 :: dx, ss, dd, x
  integer :: mm, im

  mm = 1000
  dx = (t2-t1)/real(2*mm)
  
  ss = 0.0
  if (t2 .ne. 0.0) then
     do im = 1, mm
        x = (2*im-1)*dx+t1
        dd = (x**n)/(exp(x)-1.0)
        ss = ss+4.0*dd*dx/3.0
     enddo
     do im = 1, mm-1
        x = 2*im*dx+t1
        dd = (x**n)/(exp(x)-1.0)
        ss = ss+2.0*dd*dx/3.0
     enddo
     dd = (t2**n)/(exp(t2)-1.0)
     ss = ss+dx*dd/3.0
  endif
  !x = (t2+t1)/2d0
  !dd = (x**n)/(exp(x)-1d0)
  !ss = (t2-t1)*dd

  specint = ss


end function specint
