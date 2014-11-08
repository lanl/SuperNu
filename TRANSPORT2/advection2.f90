subroutine advection2(pretrans,ig,zrsrc,zzsrc,rsrc,zsrc)
  use timestepmod
  use gasgridmod
  use particlemod
  use inputparmod
  implicit none
  logical,intent(in) :: pretrans
  integer,intent(in) :: ig
  integer,intent(inout) :: zrsrc,zzsrc
  real*8,intent(inout) :: rsrc,zsrc
!-----------------------------------------------------------------------
! This routine computes the advection of IMC particles through the
! velocity grid in cylindrical geometry.
!-----------------------------------------------------------------------
!-- advection split parameter
  real*8,parameter :: alph2 = .5d0
  logical,parameter :: partstopper = .true.
!
  integer,external :: binsrch
  integer :: zzholder,zrholder
  real*8 :: help, rold, zold, rtil, ztil
  integer :: ir,iz
!-- statement functions
  integer :: l
  real*8 :: dx,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

!-- storing initial position
  rold = rsrc
  zold = zsrc
!-- setting tentative new position
  if(pretrans) then
     rsrc = rsrc*tsp_t/(tsp_t+alph2*tsp_dt)
     zsrc = zsrc*tsp_t/(tsp_t+alph2*tsp_dt)
  else
     rsrc = rsrc*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
     zsrc = zsrc*(tsp_t+alph2*tsp_dt)/(tsp_t+tsp_dt)
  endif

  if(rsrc<grd_xarr(zrsrc).or.abs(zsrc) < &
       min(abs(grd_yarr(zzsrc)),abs(grd_yarr(zzsrc+1)))) then
!-- finding tentative new index
     zrholder = binsrch(rsrc,grd_xarr,grd_nx+1,0)
     zzholder = binsrch(zsrc,grd_yarr,grd_ny+1,0)
!-- checking if DDMC is active
     if(.not.in_puretran.and.partstopper) then
        ir = zrsrc
        iz = zzsrc

        if(iz>grd_ny/2) then
!-- z>0
           do while(iz>zzholder.or.ir>zrholder)
!-- moving towards tentative index
              if(zold>0d0) then
                 rtil = (rold/zold)*grd_yarr(iz)
              else
                 rtil = rsrc
              endif
              if(rold>0d0) then
                 ztil = (zold/rold)*grd_xarr(ir)
              else
                 ztil = zsrc
              endif
!
              if(rtil<grd_xarr(ir)) then

!-- meeting inner radial bound
                 if((grd_sig(ir-1,iz,1)+grd_cap(ig,ir-1,iz,1)) * &
                      min(dy(iz),dx(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    rsrc = grd_xarr(ir)
                    zsrc = ztil
                    exit
                 else
                    ir=ir-1
                 endif
              elseif(ztil<grd_yarr(iz)) then

!-- meeting lower z-axis bound
                 if((grd_sig(ir,iz-1,1)+grd_cap(ig,ir,iz-1,1)) * &
                      min(dy(iz-1),dx(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    zsrc = grd_yarr(iz)
                    rsrc = rtil
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
           do while(iz<zzholder.or.ir>zrholder)
!-- moving towards tentative index
              if(zold<0d0) then
                 rtil = (rold/zold)*grd_yarr(iz+1)
              else
                 rtil = rsrc
              endif
              if(rold>0d0) then
                 ztil = (zold/rold)*grd_xarr(ir)
              else
                 ztil = zsrc
              endif
!
              if(rtil<grd_xarr(ir)) then

!-- meeting inner radial bound
                 if((grd_sig(ir-1,iz,1)+grd_cap(ig,ir-1,iz,1)) * &
                      min(dy(iz),dx(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    rsrc = grd_xarr(ir)
                    zsrc = ztil
                    exit
                 else
                    ir=ir-1
                 endif
              elseif(ztil>grd_yarr(iz+1)) then

!-- meeting upper z-axis bound (z<0)
                 if((grd_sig(ir,iz+1,1)+grd_cap(ig,ir,iz+1,1)) * &
                      min(dy(iz+1),dx(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    zsrc = grd_yarr(iz+1)
                    rsrc = rtil
                    exit
                 else
                    iz = iz+1
                 endif
              endif
           enddo
        endif
        zrsrc = ir
        zzsrc = iz

     else
!-- DDMC inactive, setting index to tentative value
        zrsrc = zrholder
        zzsrc = zzholder
     endif
  endif

end subroutine advection2
