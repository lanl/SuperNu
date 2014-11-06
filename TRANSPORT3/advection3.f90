subroutine advection3(pretrans,ig,ix,iy,iz,x,y,z)
  use timestepmod
  use gasgridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  integer,intent(in) :: ig
  integer,intent(inout) :: ix,iy,iz
  real*8,intent(inout) :: x,y,z
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid in 3D planar geometry.
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer,external :: binsrch
  integer :: ixholder,iyholder,izholder
  real*8 :: help, xold,yold,zold,xtil,ytil,ztil
  integer :: i,j,k
  integer :: imove,nmove
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)

!-- storing initial position
  xold = x
  yold = y
  zold = z
!-- setting tentative new position
  if(pretrans) then
     x = x*tsp_t/(tsp_t+alph2*tsp_dt)
     y = y*tsp_t/(tsp_t+alph2*tsp_dt)
     z = z*tsp_t/(tsp_t+alph2*tsp_dt)
  else
     x = x*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
     y = y*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
     z = z*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

  if(abs(x)<min(abs(gas_xarr(ix)),abs(gas_xarr(ix+1))) .or. &
       abs(y)<min(abs(gas_yarr(iy)),abs(gas_yarr(iy+1))) .or. &
       abs(z)<min(abs(gas_zarr(iz)),abs(gas_zarr(iz+1)))) then
!-- finding tentative new index
     ixholder = binsrch(x,gas_xarr,gas_nx+1,0)
     iyholder = binsrch(y,gas_yarr,gas_ny+1,0)
     izholder = binsrch(z,gas_zarr,gas_nz+1,0)
!-- checking if DDMC is active
     if(.not.in_puretran.and.partstopper) then
!-- initializing tracking cells
        i = ix
        j = iy
        k = iz
!-- number of cell moves
        nmove = abs(ix-ixholder)+abs(iy-iyholder) + &
             abs(iz-izholder)
        do imove=1,nmove

        enddo

     else
!-- DDMC inactive, setting index to tentative value
        ix = ixholder
        iy = iyholder
        iz = izholder
     endif
  endif

end subroutine advection3
