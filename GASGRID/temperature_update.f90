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
  integer :: i
  real*8 :: dtemp, dtemp2

  !calculating radiation energy density
  dd_eraddens = dd_eraddens/dd_vol
  dd_eraddens = dd_eraddens

!!-- calculating temperature
!  if(allocated(gas_temppreset)) then
!!-- apply read-in temperature profile
!   dd_temp = gas_temppreset(:,:,:,tsp_it)
!  else
!!-- calculate temp correction
     do i=1,dd_ncell
        if(dd_bcoef(i)>0d0) then
           dtemp = dd_edep(i)/dd_vol(i) !new
           dtemp = (dtemp - tsp_dt*dd_fcoef(i)*dd_siggrey(i)* &
                pc_c*dd_ur(i))/dd_bcoef(i)
           dtemp2 = (dd_fcoef(i)/dd_bcoef(i))*tsp_dt* &
                dd_matsrc(i)
           dd_temp(i) = dd_temp(i)+dtemp+dtemp2
        elseif(dd_rho(i)>0d0.and.any(dd_cap(:,i)>0d0)) then
           dd_temp(i) = max(1000d0,(dd_eraddens(i)/pc_acoef)**(.25d0))
        else
!-- void?
           dd_temp(i) = 1000d0
        endif
     enddo !i
!  endif

  dd_ur = pc_acoef*dd_temp**4
!
!-- summing comoving material energy
  dd_emat = sum(dd_bcoef*dd_temp*dd_vol)

end subroutine temperature_update
