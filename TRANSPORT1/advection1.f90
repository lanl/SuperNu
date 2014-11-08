subroutine advection1(pretrans,ig,zsrc,rsrc)
  use timestepmod
  use gasgridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  integer,intent(in) :: ig
  integer,intent(inout) :: zsrc
  real*8,intent(inout) :: rsrc
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid.  It is geometry dependent
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer,external :: binsrch
  integer :: zholder,zfdiff
  real*8 :: help
  integer :: ir
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
!
!-- different values are used before and after transport
  if(pretrans) then
    rsrc = rsrc*tsp_t/(tsp_t+alph2*tsp_dt)
  else
    rsrc = rsrc*(tsp_t + alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif
!
  if (rsrc < grd_xarr(zsrc)) then
!
    zholder = binsrch(rsrc,grd_xarr,grd_nx+1,0)
!
    if(.not.in_puretran.and.partstopper) then
       zfdiff = -1
       if(grd_isvelocity) then
          help = tsp_t
       else
          help = 1d0
       endif
       do ir = zsrc-1,zholder,-1
          if((grd_sig(ir,1,1)+grd_cap(ig,ir,1,1))*dx(ir) &
               *help>=prt_tauddmc) then
             zfdiff = ir
             exit
          endif
       enddo
       if(zfdiff.ne.-1) then
!--
          zsrc = zfdiff+1
          rsrc = grd_xarr(zsrc)
!--
       else
          zsrc = zholder
       endif
     else
       zsrc = zholder
     endif
!
  endif
end subroutine advection1
