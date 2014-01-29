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
  !write(*,*) gas_eraddens(:,5), gas_temp(5), gas_vals2(5)%bcoef, gas_nvol(5)
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

  !calculating temperature
!  if(tsp_it/=1) then
  do ir = 1, gas_nr
     dtemp = gas_edep(ir)/gas_vals2(ir)%vol !new
     !write(6,*) gas_edep(ir), gas_vals2(ir)%vol
     dtemp = (dtemp - tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)*pc_c*gas_vals2(ir)%ur)/gas_vals2(ir)%bcoef
     
     dtemp2 = (gas_fcoef(ir)/gas_vals2(ir)%bcoef)*tsp_dt* &
          gas_vals2(ir)%matsrc

     gas_temp(ir)=gas_temp(ir)+dtemp+dtemp2
     !gas_temp(ir)=gas_temp(ir)*1d0
     if(.not.in_isbdf2.or.tsp_it==1) then
        gas_tempold(ir)=gas_temp(ir)
     else
        temphelp = gas_tempold(ir)/3d0
        gas_tempold(ir)=gas_temp(ir)
        gas_temp(ir)=4d0*gas_tempold(ir)/3d0-temphelp
     endif

     gas_siggreyold(ir)=gas_siggrey(ir)
     gas_vals2(ir)%ur = pc_acoef*gas_temp(ir)**4
     
  enddo
!  endif
!
!-- summing comoving material energy
  gas_emat = 0d0
  do ir = 1, gas_nr
     gas_emat=gas_emat+ &
          gas_vals2(ir)%bcoef*gas_tempold(ir)* &
          gas_vals2(ir)%vol
  enddo

end subroutine temperature_update
