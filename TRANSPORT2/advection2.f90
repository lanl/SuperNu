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

  if(rsrc<gas_rarr(zrsrc).or.zsrc<gas_zarr(zzsrc)) then
!-- finding tentative new index
     zrholder = binsrch(rsrc,gas_rarr,gas_nr+1,0)
     zzholder = binsrch(zsrc,gas_zarr,gas_nz+1,0)
!-- checking if DDMC is active
     if(.not.in_puretran.and.partstopper) then
        ir = zrsrc
        iz = zzsrc

        if(iz>gas_nz/2) then
!-- z>0
           do while(iz>zzholder.or.ir>zrholder)
!-- moving towards tentative index
              if(zold>0d0) then
                 rtil = (rold/zold)*gas_zarr(iz)
              else
                 rtil = gas_rarr(ir+1)
              endif
              if(rold>0d0) then
                 ztil = (zold/rold)*gas_rarr(ir)
              else
                 ztil = gas_zarr(iz+1)
              endif
!
              if(rtil<gas_rarr(ir)) then

!-- meeting inner radial bound
                 if((gas_sig(ir-1,iz)+gas_cap(g,ir-1,iz)) * &
                      min(gas_dzarr(iz),gas_drarr(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    rsrc = gas_rarr(ir)
                    zsrc = ztil
                    exit
                 else
                    ir=ir-1
                 endif
              else

!-- meeting lower z-axis bound
                 if((gas_sig(ir,iz-1)+gas_cap(g,ir,iz-1)) * &
                      min(gas_dzarr(iz-1),gas_drarr(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    zsrc = gas_zarr(iz)
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
                 rtil = (rold/zold)*gas_zarr(iz+1)
              else
                 rtil = gas_rarr(ir+1)
              endif
              if(rold>0d0) then
                 ztil = (zold/rold)*gas_rarr(ir)
              else
                 ztil = gas_zarr(iz)
              endif
!
              if(rtil<gas_rarr(ir)) then

!-- meeting inner radial bound
                 if((gas_sig(ir-1,iz)+gas_cap(g,ir-1,iz)) * &
                      min(gas_dzarr(iz),gas_drarr(ir-1)) * &
                      tsp_t >= prt_tauddmc) then
                    rsrc = gas_rarr(ir)
                    zsrc = ztil
                    exit
                 else
                    ir=ir-1
                 endif
              elseif(ztil>gas_zarr(iz+1)) then

!-- meeting upper z-axis bound (z<0)
                 if((gas_sig(ir,iz+1)+gas_cap(g,ir,iz+1)) * &
                      min(gas_dzarr(iz+1),gas_drarr(ir)) * &
                      tsp_t >= prt_tauddmc) then
                    zsrc = gas_zarr(iz+1)
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
