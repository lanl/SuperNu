subroutine advection3(pretrans,ig,ix,iy,iz,x,y,z)
  use timestepmod
  use gridmod
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
  real*8 :: dx,dy,dz,xmag,ymag,zmag
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xmag(l) = min(abs(grd_xarr(l)),abs(grd_xarr(l+1)))
  ymag(l) = min(abs(grd_yarr(l)),abs(grd_yarr(l+1)))
  zmag(l) = min(abs(grd_zarr(l)),abs(grd_zarr(l+1)))

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

!-- nothing to do
  if(abs(x)>=xmag(ix) .and. abs(y)>=ymag(iy) .and. abs(z)>=zmag(iz)) return

!
!-- sanity check
  if(xold==0d0.and.yold==0d0.and.zold==0d0) &
       stop 'advection3: invalid position update'

!-- finding tentative new index
  ixholder = binsrch(x,grd_xarr,grd_nx+1,0)
  iyholder = binsrch(y,grd_yarr,grd_ny+1,0)
  izholder = binsrch(z,grd_zarr,grd_nz+1,0)
!-- correcting new index
  if(x==0d0) ixholder = ix !-- on yz plane
  if(x<0d0.and.grd_xarr(ix)==x) ixholder=ixholder-1 !-- moved to negative x-plane (rare)
  if(y==0d0) iyholder = iy !-- on xz plane
  if(y<0d0.and.grd_yarr(iy)==y) iyholder=iyholder-1 !-- moved to negative y-plane (rare)
  if(z==0d0) izholder = iz !-- on xy plane
  if(z<0d0.and.grd_zarr(iz)==z) izholder=izholder-1 !-- moved to negative z-plane (rare)

!
!-- quick exit if DDMC is active
  if(in_puretran .or. .not.partstopper) then
!-- DDMC inactive, setting index to tentative value
     ix = ixholder
     iy = iyholder
     iz = izholder
     return
  endif

!
!-- DDMC active
!-- initializing tracking cells
  i = ix
  j = iy
  k = iz
  help = 0d0

!-- number of cell moves
  nmove = abs(ix-ixholder) + abs(iy-iyholder) + abs(iz-izholder)

  do imove=1,nmove
!-- speed at grid planes
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
     if(help == rx) then
        if(xmag(i)==abs(grd_xarr(i+1))) then!{{{
!-- x<0
           if((grd_sig(i+1,j,k)+grd_cap(ig,i+1,j,k)) * &
                min(dy(j),dx(i+1),dz(k))*tsp_t >= &
                prt_tauddmc) then
              x = grd_xarr(i+1)
              y = (yold/xold)*grd_xarr(i+1)
              z = (zold/xold)*grd_xarr(i+1)
              exit
           else
              i = i+1
           endif
        else
!-- x>0
           if((grd_sig(i-1,j,k)+grd_cap(ig,i-1,j,k)) * &
                min(dy(j),dx(i-1),dz(k))*tsp_t >= &
                prt_tauddmc) then
              x = grd_xarr(i)
              y = (yold/xold)*grd_xarr(i)
              z = (zold/xold)*grd_xarr(i)
              exit
           else
              i = i-1
           endif
        endif
!}}}
!-- y-plane
     elseif(help == ry) then
        if(ymag(j)==abs(grd_yarr(j+1))) then!{{{
!-- y<0
           if((grd_sig(i,j+1,k)+grd_cap(ig,i,j+1,k)) * &
                min(dy(j+1),dx(i),dz(k))*tsp_t >= &
                prt_tauddmc) then
              x = (xold/yold)*grd_yarr(j+1)
              y = grd_yarr(j+1)
              z = (zold/yold)*grd_yarr(j+1)
              exit
           else
              j = j+1
           endif
        else
!-- y>0
           if((grd_sig(i,j-1,k)+grd_cap(ig,i,j-1,k)) * &
                min(dy(j-1),dx(i),dz(k))*tsp_t >= &
                prt_tauddmc) then
              x = (xold/yold)*grd_yarr(j)
              y = grd_yarr(j)
              z = (zold/yold)*grd_yarr(j)
              exit
           else
              j = j-1
           endif
        endif
!}}}
!-- z-plane
     elseif(help == rz) then
        if(zmag(k)==abs(grd_zarr(k+1))) then
!-- z<0
           if((grd_sig(i,j,k+1)+grd_cap(ig,i,j,k+1)) * &
                min(dy(j),dx(i),dz(k+1))*tsp_t >= &
                prt_tauddmc) then
              x = (xold/zold)*grd_zarr(k+1)
              y = (yold/zold)*grd_zarr(k+1)
              z = grd_zarr(k+1)
              exit
           else
              k = k+1
           endif
        else
!-- z>0
           if((grd_sig(i,j,k-1)+grd_cap(ig,i,j,k-1)) * &
                min(dy(j),dx(i),dz(k-1))*tsp_t >= &
                prt_tauddmc) then
              x = (xold/zold)*grd_zarr(k)
              y = (yold/zold)*grd_zarr(k)
              z = grd_zarr(k)
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

end subroutine advection3
