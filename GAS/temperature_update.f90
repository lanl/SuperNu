!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine temperature_update

  use gasmod
  use gridmod
  use totalsmod
  use timestepmod
  use physconstmod
  use manufacmod
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

!-- calculating temperature
  if(allocated(grd_temppreset)) then
!-- apply read-in temperature profile
     gas_temp = grd_temppreset(grd_idd1:grd_idd1+grd_ndd-1,tsp_it)
  else
!!-- calculate temp correction
     do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        if(gas_bcoef(i)>0d0) then
           dtemp = gas_edep(i)/gas_vol(i) !new
           dtemp = (dtemp - tsp_dt*gas_fcoef(i)*gas_capgrey(i)* &
                pc_c*gas_ur(i))/gas_bcoef(i)
           dtemp2 = (gas_fcoef(i)/gas_bcoef(i))*tsp_dt* &
                gas_matsrc(i)
           gas_temp(i) = gas_temp(i)+dtemp+dtemp2
        elseif(gas_capgrey(i)>0d0) then
           gas_temp(i) = max(1000d0,(gas_eraddens(i)/pc_acoef)**(.25d0))
        else
!-- void?
           gas_temp(i) = 1000d0
        endif
     enddo !i
  endif

!-- total comoving material energy
  tot_emat = sum(gas_bcoef*gas_temp*gas_vol)

end subroutine temperature_update
! vim: fdm=marker
