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
  integer :: ir, ig
  real*8 :: dtemp, dtemp2
  real*8,parameter :: tauni = 8.8d0*86400.0d0
  real*8,parameter :: tauco = 111.3d0*86400.0d0
  real*8 :: temphelp

  !calculating radiation energy density
  do ir = 1, gas_nr
     gas_vals2(ir)%eraddens = 0d0
     do ig = 1, gas_ng
        gas_eraddens(ig,ir)=gas_eraddens(ig,ir)/gas_vals2(ir)%vol
        gas_vals2(ir)%eraddens = gas_vals2(ir)%eraddens+ &
             gas_eraddens(ig,ir)
     enddo
  enddo

  !calculating deposition estimate
  if(gas_depestimate) then
     gas_edep = 0d0
     do ir = 1, gas_nr
        do ig=1,gas_ng
           gas_edep(ir) = gas_edep(ir)+ &
             pc_c*tsp_dt*gas_fcoef(ir)*gas_cap(ig,ir)* &
             gas_eraddens(ig,ir)*gas_vals2(ir)%vol
        enddo
     enddo
  endif

!-- calculating temperature
  do ir = 1, gas_nr
     dtemp = gas_edep(ir)/gas_vals2(ir)%vol !new
     !write(6,*) gas_edep(ir), gas_vals2(ir)%vol
     dtemp = (dtemp - tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)*pc_c*gas_vals2(ir)%ur)/gas_vals2(ir)%bcoef
     
     dtemp2 = (gas_fcoef(ir)/gas_vals2(ir)%bcoef)*tsp_dt* &
          gas_vals2(ir)%matsrc

     gas_temp(ir)=gas_tempold(ir)+dtemp+dtemp2

     if(in_isbdf2.and.tsp_it>1) then
        temphelp = gas_tempold(ir)/3d0
        gas_tempold(ir)=gas_temp(ir)
        gas_temp(ir)=4d0*gas_tempold(ir)/3d0-temphelp
     elseif(.not.in_isbdf2) then
        gas_tempold(ir)=gas_temp(ir)
     endif

     gas_vals2(ir)%ur = pc_acoef*gas_temp(ir)**4

  enddo

!-- reset physical gas_siggreyold (used in fleck_factor, only approximate)
  if(gas_isvelocity.and.in_opacanaltyp=='none') then
     gas_siggreyold=gas_siggrey*(tsp_t/(tsp_t+tsp_dt))**3
  endif
!
!-- summing comoving material energy
  gas_emat = 0d0
  do ir = 1, gas_nr
     gas_emat=gas_emat+ &
          gas_vals2(ir)%bcoef*gas_temp(ir)* &
          gas_vals2(ir)%vol
  enddo

end subroutine temperature_update
