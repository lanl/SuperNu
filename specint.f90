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
  if (t .ne. 0.0) then
     do im = 1, mm
        x = (2*im-1)*dx
        dd = (x**n)/(exp(x)-1.0)
        ss = ss+4.0*dd*dx/3.0
     enddo
     do im = 1, mm-1
        x = 2*im*dx
        dd = (x**n)/(exp(x)-1.0)
        ss = ss+2.0*dd*dx/3.0
     enddo
     dd = (t**n)/(exp(t)-1.0)
     ss = ss+dx*dd/3.0
  endif

  specint = ss


end function specint
