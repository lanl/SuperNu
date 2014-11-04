subroutine temperature_update

  use gasgridmod
  use timestepmod
  use physconstmod
  use manufacmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine updates the material state.  At the moment
  !it updates material temperature and the approximate amount
  !of gamma ray energy introduced in the time step.
!##################################################
  integer :: i,j,k
  real*8 :: dtemp, dtemp2

  !calculating radiation energy density
  gas_eraddens=gas_eraddens/dd_vol
  dd_eraddens = gas_eraddens

!-- calculating temperature
  if(allocated(gas_temppreset)) then
!-- apply read-in temperature profile
   dd_temp = gas_temppreset(:,:,:,tsp_it)
  else
!-- calculate temp correction
     do k=1,gas_nz
     do j=1,gas_ny
     do i=1,gas_nx
        if(dd_bcoef(i,j,k)>0d0) then
           dtemp = gas_edep(i,j,k)/dd_vol(i,j,k) !new
           dtemp = (dtemp - tsp_dt*gas_fcoef(i,j,k)*gas_siggrey(i,j,k)* &
                pc_c*dd_ur(i,j,k))/dd_bcoef(i,j,k)
           dtemp2 = (gas_fcoef(i,j,k)/dd_bcoef(i,j,k))*tsp_dt* &
                dd_matsrc(i,j,k)
           dd_temp(i,j,k) = dd_temp(i,j,k)+dtemp+dtemp2
        elseif(dd_rho(i,j,k)>0d0.and.any(gas_cap(:,i,j,k)>0d0)) then
           dd_temp(i,j,k) = max(1000d0,(gas_eraddens(i,j,k)/pc_acoef)**(.25d0))
        else
!-- void?
           dd_temp(i,j,k) = 1000d0
        endif
     enddo !i
     enddo !j
     enddo !k
  endif

  dd_ur = pc_acoef*dd_temp**4
!
!-- summing comoving material energy
  gas_emat = sum(dd_bcoef*dd_temp*dd_vol)

end subroutine temperature_update
