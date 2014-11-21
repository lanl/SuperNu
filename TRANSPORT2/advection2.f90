subroutine advection2(pretrans,ig,ix,iy,x,y)
  use timestepmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  integer,intent(in) :: ig
  integer,intent(inout) :: ix,iy
  real*8,intent(inout) :: x,y
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid in cylindrical geometry.
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer,external :: binsrch
  integer :: iyholder,ixholder
  real*8 :: xold, yold, xtil, ytil
  integer :: ir,iz
!-- statement functions
  integer :: l
  real*8 :: dx,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

!-- storing initial position
  xold = x
  yold = y
!-- setting tentative new position
  if(pretrans) then
     x = x*tsp_t/(tsp_t+alph2*tsp_dt)
     y = y*tsp_t/(tsp_t+alph2*tsp_dt)
  else
     x = x*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
     y = y*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

  if(x<grd_xarr(ix).or.abs(y) < &
       min(abs(grd_yarr(iy)),abs(grd_yarr(iy+1)))) then
!-- finding tentative new index
     ixholder = binsrch(x,grd_xarr,grd_nx+1,0)
     iyholder = binsrch(y,grd_yarr,grd_ny+1,0)
!-- checking if DDMC is active
     if(.not.in_puretran.and.partstopper) then
        ir = ix
        iz = iy

        if(iz>grd_ny/2) then
!-- z>0
           do while(iz>iyholder.or.ir>ixholder)
!-- moving towards tentative index
              if(yold>0d0) then
                 xtil = (xold/yold)*grd_yarr(iz)
              else
                 xtil = x
              endif
              if(xold>0d0) then
                 ytil = (yold/xold)*grd_xarr(ir)
              else
                 ytil = y
              endif
!
              if(xtil<grd_xarr(ir)) then

!-- meeting inner radial bound
                 if((grd_sig(ir-1,iz,1)+grd_cap(ig,ir-1,iz,1)) * &
                      min(dy(iz),dx(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    x = grd_xarr(ir)
                    y = ytil
                    exit
                 else
                    ir=ir-1
                 endif
              elseif(ytil<grd_yarr(iz)) then

!-- meeting lower z-axis bound
                 if((grd_sig(ir,iz-1,1)+grd_cap(ig,ir,iz-1,1)) * &
                      min(dy(iz-1),dx(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    y = grd_yarr(iz)
                    x = xtil
                    exit
                 else
                    iz = iz-1
                 endif
              endif
           enddo
!
        else

!
!-- z<0
           do while(iz<iyholder.or.ir>ixholder)
!-- moving towards tentative index
              if(yold<0d0) then
                 xtil = (xold/yold)*grd_yarr(iz+1)
              else
                 xtil = x
              endif
              if(xold>0d0) then
                 ytil = (yold/xold)*grd_xarr(ir)
              else
                 ytil = y
              endif
!
              if(xtil<grd_xarr(ir)) then

!-- meeting inner radial bound
                 if((grd_sig(ir-1,iz,1)+grd_cap(ig,ir-1,iz,1)) * &
                      min(dy(iz),dx(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    x = grd_xarr(ir)
                    y = ytil
                    exit
                 else
                    ir=ir-1
                 endif
              elseif(ytil>grd_yarr(iz+1)) then

!-- meeting upper z-axis bound (z<0)
                 if((grd_sig(ir,iz+1,1)+grd_cap(ig,ir,iz+1,1)) * &
                      min(dy(iz+1),dx(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    y = grd_yarr(iz+1)
                    x = xtil
                    exit
                 else
                    iz = iz+1
                 endif
              endif
           enddo
        endif
        ix = ir
        iy = iz

     else
!-- DDMC inactive, setting index to tentative value
        ix = ixholder
        iy = iyholder
     endif
  endif

end subroutine advection2
