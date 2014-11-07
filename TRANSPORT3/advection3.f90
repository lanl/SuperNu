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
  real*8 :: rold, xold,yold,zold,rx,ry,rz
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
        help = 0d0
!-- number of cell moves
        nmove = abs(ix-ixholder)+abs(iy-iyholder) + &
             abs(iz-izholder)
        do imove=1,nmove

!-- radial displacements toward grid planes
           if(xold==0d0) then
              rx = 0d0
           else
              rx = rold*xmag(i)/abs(xold)
           endif
           if(yold==0d0) then
              ry = 0d0
           else
              ry = rold*ymag(j)/abs(yold)
           endif
           if(zold==0d0) then
              rz = 0d0
           else
              rz = rold*zmag(k)/abs(zold)
           endif

!-- using min radial displacement to move index
           help = max(rx,ry,rz)

!-- x-plane
           if(rx == help) then
              if(xmag(i)==abs(gas_xarr(i+1))) then
!-- x<0
                 if((gas_sig(i+1,j,k)+gas_cap(ig,i+1,j,k)) * &
                      min(dy(j),dx(i+1),dz(k))*tsp_t >= &
                      prt_tauddmc) then
                    x = gas_xarr(i+1)
                    y = (yold/xold)*gas_xarr(i+1)
                    z = (zold/xold)*gas_xarr(i+1)
                    exit
                 else
                    i = i+1
                 endif
              else
!-- x>0
                 if((gas_sig(i-1,j,k)+gas_cap(ig,i-1,j,k)) * &
                      min(dy(j),dx(i-1),dz(k))*tsp_t >= &
                      prt_tauddmc) then
                    x = gas_xarr(i)
                    y = (yold/xold)*gas_xarr(i)
                    z = (zold/xold)*gas_xarr(i)
                    exit
                 else
                    i = i-1
                 endif
              endif

!-- y-plane
           elseif(ry == help) then
              if(ymag(j)==abs(gas_yarr(j+1))) then
!-- y<0
                 if((gas_sig(i,j+1,k)+gas_cap(ig,i,j+1,k)) * &
                      min(dy(j+1),dx(i),dz(k))*tsp_t >= &
                      prt_tauddmc) then
                    x = (xold/yold)*gas_yarr(j+1)
                    y = gas_yarr(j+1)
                    z = (zold/yold)*gas_yarr(j+1)
                    exit
                 else
                    j = j+1
                 endif
              else
!-- y>0
                 if((gas_sig(i,j-1,k)+gas_cap(ig,i,j-1,k)) * &
                      min(dy(j-1),dx(i),dz(k))*tsp_t >= &
                      prt_tauddmc) then
                    x = (xold/yold)*gas_yarr(j)
                    y = gas_yarr(j)
                    z = (zold/yold)*gas_yarr(j)
                    exit
                 else
                    j = j-1
                 endif
              endif

!-- z-plane
           elseif(rz == help) then
              if(zmag(k)==abs(gas_zarr(k+1))) then
!-- z<0
                 if((gas_sig(i,j,k+1)+gas_cap(ig,i,j,k+1)) * &
                      min(dy(j),dx(i),dz(k+1))*tsp_t >= &
                      prt_tauddmc) then
                    x = (xold/zold)*gas_zarr(k+1)
                    y = (yold/zold)*gas_zarr(k+1)
                    z = gas_zarr(k+1)
                    exit
                 else
                    k = k+1
                 endif
              else
!-- z>0
                 if((gas_sig(i,j,k-1)+gas_cap(ig,i,j,k-1)) * &
                      min(dy(j),dx(i),dz(k-1))*tsp_t >= &
                      prt_tauddmc) then
                    x = (xold/zold)*gas_zarr(k)
                    y = (yold/zold)*gas_zarr(k)
                    z = gas_zarr(k)
                    exit
                 else
                    k = k-1
                 endif
              endif
           endif
        enddo
        ix = i
        iy = j
        iz = k

     else
!-- DDMC inactive, setting index to tentative value
        ix = ixholder
        iy = iyholder
        iz = izholder
     endif
  endif

end subroutine advection3
