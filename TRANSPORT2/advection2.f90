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
  integer :: zzholder,zrholder,zfdiff
  real*8 :: help, rold, zold
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
        zfdiff = -1
!-- moving particle towards tentative position
        do while((ir/=zrholder.or.iz/=zzholder).and. &
             zfdiff==-1)
!-- checking nearest boundary (r or z)
           if(max(gas_rarr(ir)/rold,abs(gas_zarr(iz+1)/zold))== &
                gas_rarr(ir)/rold) then
!-- checking if cell at inner radius is DDMC
              if((gas_sig(ir-1,iz)+gas_cap(g,ir-1,iz))* &
                   min(gas_drarr(ir-1),gas_dzarr(iz))*help>=&
                   prt_tauddmc) then
!-- placing particle on outer surface of DDMC cell
                 zfdiff = ir-1
                 zrsrc = zfdiff+1
                 zsrc = zold*(1d0-gas_rarr(zrsrc)/rold)
                 rsrc = gas_rarr(zrsrc)
              else
!-- updating ir index and position to inner cell
                 ir = ir-1
                 zold = zold*(1d0-gas_rarr(ir+1)/rold)
                 rold = gas_rarr(ir+1)
              endif
           else
!-- checking if z<0
              if(iz<gas_nz/2) then
!-- checking if upper cell is DDMC
                 if((gas_sig(ir,iz+1)+gas_cap(g,ir,iz+1))* &
                      min(gas_drarr(ir),gas_dzarr(iz+1))*help>=&
                      prt_tauddmc) then
!-- placing particle on lower surface of DDMC cell
                    zfdiff = iz+1
                    zzsrc = zfdiff-1
                    rsrc = rold*(1d0-abs(gas_zarr(zzsrc)/zold))
                    zsrc = gas_zarr(zzsrc+1)
                 else
!-- updating iz index and position to upper cell
                    iz = iz+1
                    rold = rold*(1d0-abs(gas_zarr(iz)/zold))
                    zold = gas_zarr(iz)
                 endif
!-- else z>0
              else
!-- checking if lower cell is DDMC
                 if((gas_sig(ir,iz-1)+gas_cap(g,ir,iz-1))* &
                      min(gas_drarr(ir),gas_dzarr(iz-1))*help>=&
                      prt_tauddmc) then
!-- placing particle on upperer surface of DDMC cell
                    zfdiff = iz-1
                    zzsrc = zfdiff+1
                    rsrc = rold*(1d0-abs(gas_zarr(zzsrc)/zold))
                    zsrc = gas_zarr(zzsrc)
                 else
!-- updating iz index and position to lower cell
                    iz = iz-1
                    rold = rold*(1d0-abs(gas_zarr(iz+1)/zold))
                    zold = gas_zarr(iz+1)
                 endif
              endif
           endif
        enddo
!-- if no DDMC encountered, setting index to tentative value
        if(zfdiff==-1) then
           zrsrc = zrholder
           zzsrc = zzholder
        endif
     else
!-- DDMC inactive, setting index to tentative value
        zrsrc = zrholder
        zzsrc = zzholder
     endif
  endif

end subroutine advection2
