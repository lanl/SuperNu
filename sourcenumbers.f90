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

  integer :: i,j,k
  integer :: ihelp
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
  do k = 1, gas_nz
  do j = 1, gas_ny
  do i = 1, gas_nx
     
     gas_emit(i,j,k) =  tsp_dt*gas_fcoef(i,j,k)*gas_siggrey(i,j,k)*pc_c* &
          dd_ur(i,j,k)*dd_vol(i,j,k)
!
     gas_emit(i,j,k) = gas_emit(i,j,k)+&
          tsp_dt*dd_vol(i,j,k)*(1d0-gas_fcoef(i,j,k))*&
          dd_matsrc(i,j,k)
!-- accounting for total material source energy:
!-- gas_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-gas_fcoef is thermal radiation.
     if(tsp_it==1) then
        gas_eext = gas_eext+tsp_dt*dd_vol(i,j,k)*&
             dd_matsrc(i,j,k)
     else
!-- rtw: gas_eext is only broadcast for tsp_it==1
        gas_eext = gas_eext+tsp_dt*dd_vol(i,j,k)*&
             dd_matsrc(i,j,k)*dble(nmpi)
     endif
!
!-- accounting for thermal decay radiation source energy
     if(.not.gas_novolsrc .and. gas_srctype=='none') then
        gas_emit(i,j,k) = gas_emit(i,j,k) + dd_nisource(i,j,k)
        if(tsp_it==1) then
           gas_eext = gas_eext + dd_nisource(i,j,k)
        else
!-- rtw: gas_eext is only broadcast for tsp_it==1
           gas_eext = gas_eext + dble(nmpi)*dd_nisource(i,j,k)
        endif
     endif
     gas_etot = gas_etot + gas_emit(i,j,k)
  enddo
  enddo
  enddo

  ! Adding external energy to gas_etot
  gas_etot = gas_etot + sum(gas_emitex)
  
  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf = nint(gas_esurf*prt_ns/gas_etot)
  prt_nnew = prt_nsurf

  ! Calculating number of particles per cell (dd_nvol): loop
  select case(in_igeom)
  case(1)
     ihelp = 50
  case(2)
     ihelp = 50/nmpi+1
  case(3)
     ihelp = 1
  endselect
  prt_nexsrc=0
  do k = 1, gas_nz
  do j = 1, gas_ny
  do i = 1, gas_nx
     if(gas_emit(i,j,k)<=0d0) then
        gas_nvol(i,j,k)=0
     else
        gas_nvol(i,j,k)=nint(gas_emit(i,j,k)*prt_ns/gas_etot) + &
             ihelp
     endif
     prt_nnew = prt_nnew + gas_nvol(i,j,k)
     !external source volume numbers
     if(gas_emitex(i,j,k)<=0d0) then
        gas_nvolex(i,j,k)=0
     else
        gas_nvolex(i,j,k)=nint(gas_emitex(i,j,k)*prt_ns/gas_etot) + &
             ihelp
     endif
     prt_nexsrc = prt_nexsrc + gas_nvolex(i,j,k)
     prt_nnew = prt_nnew + gas_nvolex(i,j,k)
  enddo
  enddo
  enddo

end subroutine sourcenumbers
