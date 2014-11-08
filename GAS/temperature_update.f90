subroutine temperature_update

  use gasmod
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
  gas_eraddens = gas_eraddens/gas_vol
  gas_eraddens = gas_eraddens

!-- calculating temperature
  if(allocated(gas_temppreset)) then
!-- apply read-in temperature profile
   gas_temp = gas_temppreset(:,tsp_it)
  else
!!-- calculate temp correction
     do i=1,gas_ncell
        if(gas_bcoef(i)>0d0) then
           dtemp = gas_edep(i)/gas_vol(i) !new
           dtemp = (dtemp - tsp_dt*gas_fcoef(i)*gas_siggrey(i)* &
                pc_c*gas_ur(i))/gas_bcoef(i)
           dtemp2 = (gas_fcoef(i)/gas_bcoef(i))*tsp_dt* &
                gas_matsrc(i)
           gas_temp(i) = gas_temp(i)+dtemp+dtemp2
        elseif(gas_rho(i)>0d0.and.any(gas_cap(:,i)>0d0)) then
           gas_temp(i) = max(1000d0,(gas_eraddens(i)/pc_acoef)**(.25d0))
        else
!-- void?
           gas_temp(i) = 1000d0
        endif
     enddo !i
  endif

end subroutine temperature_update
