      function locate(xarray,n,x)
c     -------------------------------
      integer,intent(in) :: n
      REAL*8,intent(in) :: xarray(n),x
****************************************************************
* search in an ordered table.
* from numerical recipes
****************************************************************
      integer jl,jm,ju
c
      jl = 0
      ju = n + 1
      do
       jm = (ju + jl)/2
       if((xarray(n)>=xarray(1)) .eqv. (x>=xarray(jm))) then
        jl=jm
       else
        ju=jm
       endif
       if(ju-jl .le. 1) exit
      enddo
c--
      locate = jl
      return
      end function locate
