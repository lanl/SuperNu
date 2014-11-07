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
  real*8 :: rold, xold,yold,zold,drx,dry,drz
  real*8 :: help
  integer :: i,j,k
  integer :: imove,nmove
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xnext,ynext,znext
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)
  xmag(l) = min(abs(gas_xarr(l)),abs(gas_xarr(l+1)))
  ymag(l) = min(abs(gas_yarr(l)),abs(gas_yarr(l+1)))
  zmag(l) = min(abs(gas_zarr(l)),abs(gas_zarr(l+1)))

!-- storing initial position
  xold = x
  yold = y
  zold = z
  rold = sqrt(x**2+y**2+z**2)
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

  if(abs(x)<xmag(ix).or.abs(y)<ymag(iy) .or. &
       abs(z)<zmag(iz)) then
!
!-- sanity check
     if(xold==0d0.and.yold==0d0.and.zold==0d0) &
          stop 'advection3: invalid position update'
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

           if(xold==0d0) then
              drx = rold
           else
              drx = rold*(1d0-xmag(i)/abs(xold))
           endif
           if(yold==0d0) then
              dry = rold
           else
              dry = rold*(1d0-ymag(j)/abs(yold))
           endif
           if(zold==0d0) then
              drz = rold
           else
              drz = rold*(1d0-zmag(k)/abs(zold))
           endif

           help = min(drx,dry,drz)
           if(drx == help) then
              if(xmag(i)==abs(gas_xarr(i+1))) then
                 i = i+1
              else
                 i = i-1
              endif
           endif
           if(dry == help) then
              if(ymag(j)==abs(gas_yarr(j+1))) then
                 j = j+1
              else
                 j = j-1
              endif
           endif
           if(drz == help) then
              if(zmag(k)==abs(gas_zarr(k+1))) then
                 k = k+1
              else
                 k = k-1
              endif
           endif

        enddo

     else
!-- DDMC inactive, setting index to tentative value
        ix = ixholder
        iy = iyholder
        iz = izholder
     endif
  endif

end subroutine advection3
