subroutine sourceenergy(nmpi)

  use gasmod
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
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: i
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns
  
!-- prepare manufactured solution temperature source
  if(in_srctype=='manu') then
     call generate_manutempsrc(in_totmass,in_gas_capcoef,tsp_t,tsp_dt)
  endif

  gas_emit = 0d0
  gas_emitex = 0d0
  ! Calculating fictitious emission energy per cell: loop
  do i=1,gas_ncell
!
!-- thermal source
     gas_emit(i) =  tsp_dt*gas_fcoef(i)*gas_capgrey(i)*pc_c* &
          gas_ur(i)*gas_vol(i)
!
!-- manufactured solution material source
     gas_emit(i) = gas_emit(i) + &
          tsp_dt*gas_vol(i)*(1d0-gas_fcoef(i))*&
          gas_matsrc(i)
!-- accounting for total material source energy:
!-- gas_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-gas_fcoef is thermal radiation.
     if(tsp_it==1) then
        tot_eext = tot_eext+tsp_dt*gas_vol(i)*&
             gas_matsrc(i)
     else
!-- rtw: tot_eext is only broadcast for tsp_it==1
        tot_eext = tot_eext+tsp_dt*gas_vol(i)*&
             gas_matsrc(i)*dble(nmpi)
     endif
!
!-- non-thermal decay radiation source energy
     if(.not.in_novolsrc .and. in_srctype=='none') then
        gas_emitex(i) = gas_nisource(i)
     endif
  enddo

end subroutine sourceenergy
