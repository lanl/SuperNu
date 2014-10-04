subroutine sourcenumbers

  use mpimod
  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: ir
! gas_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

  gas_nvol = 0
  gas_nvolex = 0
  gas_esurf = 0d0
  gas_emit = 0d0
  gas_emitex = 0d0
  
  gas_etot = 0d0
  gas_esurf = 0d0

  if(gas_isvelocity) then
     gas_etot = gas_esurf*tsp_t**2
     gas_esurf = gas_esurf*tsp_t**2
  else
     gas_etot = gas_esurf
  endif
!
!  gas_eext=gas_eext+gas_etot
!
  ! Calculating gas_emitex from analytic distribution
  call analytic_source

  ! Calculating fictitious emission energy per cell: loop
  do ir = 1, gas_nx
     
     gas_emit(ir,1,1) =  tsp_dt*gas_fcoef(ir,1,1)*gas_siggrey(ir,1,1)*pc_c* &
          gas_vals2(ir,1,1)%ur*gas_vals2(ir,1,1)%vol
!
     gas_emit(ir,1,1) = gas_emit(ir,1,1)+&
          tsp_dt*gas_vals2(ir,1,1)%vol*(1d0-gas_fcoef(ir,1,1))*&
          gas_vals2(ir,1,1)%matsrc
!-- accounting for total material source energy:
!-- gas_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-gas_fcoef is thermal radiation.
     if(tsp_it==1) then
        gas_eext = gas_eext+tsp_dt*gas_vals2(ir,1,1)%vol*&
             gas_vals2(ir,1,1)%matsrc
     else
!-- rtw: gas_eext is only broadcast for tsp_it==1
        gas_eext = gas_eext+tsp_dt*gas_vals2(ir,1,1)%vol*&
             gas_vals2(ir,1,1)%matsrc*dble(nmpi)
     endif
!
!-- accounting for thermal decay radiation source energy
     if(.not.gas_novolsrc .and. gas_srctype=='none') then
        gas_emit(ir,1,1) = gas_emit(ir,1,1) + gas_vals2(ir,1,1)%nisource
        if(tsp_it==1) then
           gas_eext = gas_eext + gas_vals2(ir,1,1)%nisource
        else
!-- rtw: gas_eext is only broadcast for tsp_it==1
           gas_eext = gas_eext + dble(nmpi)*gas_vals2(ir,1,1)%nisource
        endif
     endif
     gas_etot = gas_etot + gas_emit(ir,1,1)
  enddo
  
  ! Adding external energy to gas_etot
  gas_etot = gas_etot + sum(gas_emitex)
  
  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf = nint(gas_esurf*prt_ns/gas_etot)
  prt_nnew = prt_nsurf

  ! Calculating number of particles per cell (gas_vals2%nvol): loop
  prt_nexsrc=0
  do ir = 1, gas_nx
     if(gas_emit(ir,1,1)<=0d0) then
        gas_nvol(ir,1,1)=0
     else
        gas_nvol(ir,1,1)=nint(abs(gas_emit(ir,1,1))*prt_ns/gas_etot)+50
     endif
     prt_nnew = prt_nnew + gas_nvol(ir,1,1)
     !external source volume numbers
     gas_nvolex(ir,1,1)=nint(gas_emitex(ir,1,1)*prt_ns/gas_etot)
     prt_nexsrc = prt_nexsrc + gas_nvolex(ir,1,1)
     prt_nnew = prt_nnew + gas_nvolex(ir,1,1)
  enddo

end subroutine sourcenumbers
