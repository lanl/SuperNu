subroutine temperature_update

  use gasgridmod
  use timestepmod
  use physconstmod
  use manufacmod
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
  real*8 :: ddrr3, ddrr4

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
        forall(ig=1:gas_ng) gas_edep(ir)=gas_edep(ir)+ &
             pc_c*tsp_dt*gas_fcoef(ir)*gas_cap(ig,ir)* &
             gas_eraddens(ig,ir)*gas_vals2(ir)%vol
     enddo
  endif

  !calculating temperature
  do ir = 1, gas_nr
     dtemp = gas_edep(ir)/gas_vals2(ir)%vol !new
     !write(6,*) gas_edep(ir), gas_vals2(ir)%vol
     dtemp = (dtemp - tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)*pc_c*gas_vals2(ir)%ur)/gas_vals2(ir)%bcoef
     
     !if(tsp_it==17) then
        !write(6,*) dtemp, gas_edep(ir),gas_vals2(ir)%vol, gas_vals2(ir)%bcoef
     !endif
     if(gas_srctype=='manu') then
        if(gas_isvelocity) then
           !this may cause drift
           gas_temp(ir)=gas_temp(ir)*tsp_texp/(tsp_texp+tsp_dt)
        else
           !gas_temp(ir)=0d0
           ddrr3 = gas_rarr(ir+1)**3-gas_rarr(ir)**3
           ddrr4 = gas_rarr(ir+1)**4-gas_rarr(ir)**4
           gas_temp(ir)=gas_temp(ir)+dtemp - &
                (tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)/gas_vals2(ir)%bcoef)*&
                (man_aa11-0.75d0*(man_aa11-man_aa22)*ddrr4/(gas_rarr(gas_nr+1)*ddrr3))
           write(*,*) gas_temp(ir)
           if(gas_temp(ir)<0d0) then
              gas_temp(ir)=0d0
           endif
        endif
     else
        gas_temp(ir) = gas_temp(ir) + dtemp
     endif

     gas_vals2(ir)%ur = pc_acoef*gas_temp(ir)**4
     
  enddo

end subroutine temperature_update
