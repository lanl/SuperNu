subroutine temperature_update

  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

!##################################################
  !This subroutine updates the material state.  At the moment
  !it updates material temperature and the approximate amount
  !of gamma ray energy introduced in the time step.
!##################################################
  integer :: ir, ig
  real*8 :: dtemp
  real*8,parameter :: tauni = 8.8d0*86400.0d0
  real*8,parameter :: tauco = 111.3d0*86400.0d0

  !calculating radiation energy density
  do ir = 1, gas_nr
     gas_vals2(ir)%eraddens = 0d0
     !gas_edep(ir) = 0d0
     do ig = 1, gas_ng
        gas_eraddens(ig,ir)=gas_eraddens(ig,ir)/gas_vals2(ir)%vol
        gas_vals2(ir)%eraddens = gas_vals2(ir)%eraddens+ &
             gas_eraddens(ig,ir)
        !gas_edep(ir) = gas_edep(ir)+pc_c*tsp_dt*gas_fcoef(ir)* &
        !     gas_cap(ig,ir)*gas_eraddens(ig,ir)*gas_vals2(ir)%vol
     enddo
  enddo
  !write(*,*) gas_edep(1), gas_vals2(1)%eraddens

  !calculating temperature
  do ir = 1, gas_nr

     dtemp = gas_edep(ir)/gas_vals2(ir)%vol !new
     !write(6,*) gas_edep(ir), gas_vals2(ir)%vol
     dtemp = (dtemp - tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)*pc_c*gas_vals2(ir)%ur)/gas_vals2(ir)%bcoef
     !if(tsp_it==17) then
        !write(6,*) dtemp, gas_edep(ir),gas_vals2(ir)%vol, gas_vals2(ir)%bcoef
     !endif
     if(gas_srctype=='manu') then
        !this may cause drift
        gas_vals2(ir)%tempkev=gas_vals2(ir)%tempkev*tsp_texp/(tsp_texp+tsp_dt)
     else
        gas_vals2(ir)%tempkev = gas_vals2(ir)%tempkev + dtemp
     endif

     gas_vals2(ir)%temp = gas_vals2(ir)%tempkev * 1e3*pc_ev/pc_kb  !initial guess, may be overwritten by read_temp_str
     !gas_vals2(ir)%ur=dtemp/(tsp_dt*pc_c*gas_siggrey(ir))
     !gas_vals2(ir)%tempkev = (gas_vals2(ir)%ur/pc_acoef)**(0.25d0)
     !gas_vals2(ir)%bcoef = 2.0*pc_acoef*gas_vals2(ir)%tempkev**3
     gas_vals2(ir)%ur = pc_acoef*gas_vals2(ir)%tempkev**4
     !gas_edep(ir) = gas_edep(ir)/gas_vals2(ir)%vol
  enddo

end subroutine temperature_update
