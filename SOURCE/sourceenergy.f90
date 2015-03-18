subroutine sourceenergy(nmpi)

  use gasmod
  use sourcemod
  use totalsmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  use manufacmod
  implicit none
  integer,intent(in) :: nmpi

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################
! src_nsurf = number of surface prt_particles
! src_nnew = total number of new prt_particles~=src_ns
  
!-- prepare manufactured solution temperature source
  if(in_srctype=='manu') then
     call generate_manutempsrc(in_totmass,in_gas_capcoef,tsp_t,tsp_dt)
  endif

! Calculating fictitious emission energy per cell
!-- thermal source
  gas_emit =  tsp_dt*gas_vol*gas_fcoef*gas_capgrey*pc_c*gas_ur
  tot_sthermal = sum(gas_emit)
!
!-- manufactured solution material source
  gas_emit = gas_emit + tsp_dt*gas_vol*(1d0-gas_fcoef)*gas_matsrc
  tot_smanufac = sum(tsp_dt*gas_vol*(1d0-gas_fcoef)*gas_matsrc)
!
!-- accounting for total material source energy:
!-- gas_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-gas_fcoef is thermal radiation.
  if(tsp_it==1) then
     tot_eext = tot_eext + tsp_dt*sum(gas_vol*gas_matsrc)
  else
!-- rtw: tot_eext is only broadcast for tsp_it==1
     tot_eext = tot_eext + tsp_dt*sum(gas_vol*gas_matsrc)*nmpi
  endif
!
!-- non-thermal decay radiation source energy
  if(.not.in_novolsrc .and. in_srctype=='none') then
     gas_emitex = gas_decaygamma  !grey transport
     gas_emit = gas_emit + gas_decaybeta  !local deposition
!-- totals
     tot_sdecaygamma = sum(gas_decaygamma)
     tot_sdecaybeta = sum(gas_decaybeta)
  else
     gas_emitex = 0d0
     tot_sdecaygamma = 0d0
     tot_sdecaybeta = 0d0
  endif

end subroutine sourceenergy
