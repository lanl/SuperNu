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

  integer :: i
  integer :: ihelp
! gas_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

  dd_nvol = 0
  dd_nvolex = 0
  gas_esurf = 0d0
  dd_emit = 0d0
  dd_emitex = 0d0
  
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
  ! Calculating dd_emitex from analytic distribution
  call analytic_source

  ! Calculating fictitious emission energy per cell: loop
  do i = 1, dd_ncell
     
     dd_emit(i) =  tsp_dt*dd_fcoef(i)*dd_siggrey(i)*pc_c* &
          dd_ur(i)*dd_vol(i)
!
     dd_emit(i) = dd_emit(i) + &
          tsp_dt*dd_vol(i)*(1d0-dd_fcoef(i))*&
          dd_matsrc(i)
!-- accounting for total material source energy:
!-- dd_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-dd_fcoef is thermal radiation.
     if(tsp_it==1) then
        gas_eext = gas_eext+tsp_dt*dd_vol(i)*&
             dd_matsrc(i)
     else
!-- rtw: gas_eext is only broadcast for tsp_it==1
        gas_eext = gas_eext+tsp_dt*dd_vol(i)*&
             dd_matsrc(i)*dble(nmpi)
     endif
!
!-- accounting for thermal decay radiation source energy
     if(.not.gas_novolsrc .and. gas_srctype=='none') then
        dd_emit(i) = dd_emit(i) + dd_nisource(i)
        if(tsp_it==1) then
           gas_eext = gas_eext + dd_nisource(i)
        else
!-- rtw: gas_eext is only broadcast for tsp_it==1
           gas_eext = gas_eext + dble(nmpi)*dd_nisource(i)
        endif
     endif
     gas_etot = gas_etot + dd_emit(i)
  enddo

  ! Adding external energy to gas_etot
  gas_etot = gas_etot + sum(dd_emitex)
  
  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf = nint(gas_esurf*prt_ns/gas_etot)
  if(prt_nsurf>0) stop 'sourcenumber: nsurf unsupported'
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
  do i = 1, dd_ncell
     if(dd_emit(i)<=0d0) then
        dd_nvol(i)=0
     else
        dd_nvol(i)=nint(dd_emit(i)*prt_ns/gas_etot) + &
             ihelp
     endif
     prt_nnew = prt_nnew + dd_nvol(i)
     !external source volume numbers
     if(dd_emitex(i)<=0d0) then
        dd_nvolex(i)=0
     else
        dd_nvolex(i)=nint(dd_emitex(i)*prt_ns/gas_etot) + &
             ihelp
     endif
     prt_nexsrc = prt_nexsrc + dd_nvolex(i)
     prt_nnew = prt_nnew + dd_nvolex(i)
     !write(0,*) impi,prt_nnew,dd_nvol(i),dd_nvolex(i),dd_emitex(i),gas_etot
  enddo

end subroutine sourcenumbers
