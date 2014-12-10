subroutine advection11(pretrans,ig,ix,r)
  use timestepmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  integer,intent(in) :: ig
  integer,intent(inout) :: ix
  real*8,intent(inout) :: r
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
!-- setting tentative new position
  if(pretrans) then
    r = r*tsp_t/(tsp_t+alph2*tsp_dt)
  else
    r = r*(tsp_t + alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

!
!-- nothing to do
  if (r>=grd_xarr(ix)) return

!
!-- finding tentative new index
  zholder = binsrch(r,grd_xarr,grd_nx+1,0)

!
!-- quick exit if DDMC is active
  if(in_puretran .or. .not.partstopper) then
    ix = zholder
    return
  endif

  zfdiff = -1
  if(grd_isvelocity) then
     help = tsp_t
  else
     help = 1d0
  endif
  do ir = ix-1,zholder,-1
     if((grd_sig(ir,1,1)+grd_cap(ig,ir,1,1))*dx(ir) &
          *help>=prt_tauddmc) then
        zfdiff = ir
        exit
     endif
  enddo

!--
  if(zfdiff.ne.-1) then
     ix = zfdiff+1
     r = grd_xarr(ix)
  else
     ix = zholder
  endif

end subroutine advection11
