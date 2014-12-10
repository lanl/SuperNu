subroutine advection11(pretrans,ptcl,ig)

  use miscmod
  use timestepmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  type(packet),target,intent(inout) :: ptcl
  integer,intent(in) :: ig
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid.  It is geometry dependent
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer :: zholder,zfdiff
  real*8 :: help
  integer :: i
!-- pointers
  integer,pointer :: ix
  real*8,pointer :: x 
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

  ix => ptcl%ix
  x => ptcl%x

!
!-- setting tentative new position
  if(pretrans) then
    x = x*tsp_t/(tsp_t+alph2*tsp_dt)
  else
    x = x*(tsp_t + alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

!
!-- nothing to do
  if (x>=grd_xarr(ix)) return

!
!-- finding tentative new index
  zholder = binsrch(x,grd_xarr,grd_nx+1)

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
  do i=ix-1,zholder,-1
     if((grd_sig(i,1,1)+grd_cap(ig,i,1,1))*dx(i) &
          *help>=prt_tauddmc) then
        zfdiff = i
        exit
     endif
  enddo

!--
  if(zfdiff.ne.-1) then
     ix = zfdiff+1
     x = grd_xarr(ix)
  else
     ix = zholder
  endif

end subroutine advection11
