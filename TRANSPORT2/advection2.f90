!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine advection2(pretrans,ptcl,ptcl2)

  use miscmod
  use timestepmod
  use gridmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid in cylindrical geometry.
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer :: iyholder,ixholder
  real*8 :: rold,xold,yold,rx,ry
  real*8 :: help,dz
  integer :: i,j
  integer :: imove,nmove
!-- pointers
  integer,pointer :: ix, iy, iz, ic, ig
  real*8,pointer :: x 
  real*8,pointer :: y 
!-- statement functions
  integer :: l
  real*8 :: dx,dy,ymag,xm
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  ymag(l) = min(abs(grd_yarr(l)),abs(grd_yarr(l+1)))
  xm(l) = .5d0*(grd_xarr(l+1) + grd_xarr(l))

  ix => ptcl2%ix
  iy => ptcl2%iy
  iz => ptcl2%iz
  ic => ptcl2%ic
  ig => ptcl2%ig
  x => ptcl%x
  y => ptcl%y

!-- storing initial position
  xold = x
  yold = y
  rold = sqrt(x**2+y**2)
  dz=grd_zarr(iz+1)-grd_zarr(iz)
!-- setting tentative new position
  if(pretrans) then
     x = x*tsp_t/(tsp_t+alph2*tsp_dt)
     y = y*tsp_t/(tsp_t+alph2*tsp_dt)
  else
     x = x*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
     y = y*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

!-- nothing to do
  if(x>=grd_xarr(ix) .and. abs(y)>=ymag(iy)) return

!!
!!-- sanity check
!  if(xold==0d0.and.yold==0d0) &
!       stop 'advection2: invalid position update'

!-- finding tentative new index
  ixholder = binsrch(x,grd_xarr,grd_nx+1,.false.)
  iyholder = binsrch(y,grd_yarr,grd_ny+1,.false.)
!--correcting new index
  if(x==0d0) ixholder = ix !-- on y axis
  if(y==0d0) iyholder = iy !-- on x axis
  if(y<0d0.and.grd_yarr(iyholder)==y) iyholder = iyholder-1 !-- moved to negative y-line (rare)

!
!-- quick exit if DDMC is active
  if(in_puretran .or. .not.partstopper) then
!-- DDMC inactive, setting index to tentative value
     ix = ixholder
     iy = iyholder
     ic = grd_icell(ix,iy,iz)
     return
  endif

!
!-- DDMC active
!-- initializing tracking cells
  i = ix
  j = iy
  help = 0d0

!-- number of cell moves
  nmove = ix-ixholder+abs(iy-iyholder)

  do imove=1,nmove
!-- speed at grid edges
     if(xold==0d0) then
        rx = 0d0
     else
        rx = rold*grd_xarr(i)/xold
     endif
     if(yold==0d0) then
        ry = 0d0
     else
        ry = rold*ymag(j)/abs(yold)
     endif

!-- using min displacement to move index
     help = max(rx,ry)

!-- x-edge
     if(help == rx) then
        l = grd_icell(i-1,j,iz)
        if((grd_sig(l)+grd_cap(ig,l)) * &!{{{
             min(dy(j),dx(i-1),xm(i-1)*dz)*tsp_t >= trn_tauddmc) then
           x = grd_xarr(i)
           y = (yold/xold)*grd_xarr(i)
           exit
        else
           i = i-1
        endif
!}}}
!-- y-edge
     else
        if(ymag(j)==abs(grd_yarr(j+1))) then!{{{
!-- y<0
           l = grd_icell(i,j+1,iz)
           if((grd_sig(l)+grd_cap(ig,l)) * &
                min(dy(j+1),dx(i),xm(i)*dz)*tsp_t >= &
                trn_tauddmc) then
              x = (xold/yold)*grd_yarr(j+1)
              y = grd_yarr(j+1)
              exit
           else
              j = j+1
           endif
        else
!-- y>0
           l = grd_icell(i,j-1,iz)
           if((grd_sig(l)+grd_cap(ig,l)) * &
                min(dy(j-1),dx(i),xm(i)*dz)*tsp_t >= &
                trn_tauddmc) then
              x = (xold/yold)*grd_yarr(j)
              y = grd_yarr(j)
              exit
           else
              j = j-1
           endif
        endif!}}}
     endif
  enddo
  ix = i
  iy = j
  ic = grd_icell(ix,iy,iz)

end subroutine advection2
! vim: fdm=marker
