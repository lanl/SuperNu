subroutine sourcenumbers

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

  integer :: ir, ig
  real*8 :: exsumg
  ! gas_esurf for any new prt_particles from a surface source
  ! prt_nsurf = number of surface prt_particles
  ! prt_nnew = total number of new prt_particles~=prt_ns

  gas_nvol = 0
  gas_nvolex = 0
  gas_esurf = 0d0
  gas_emit = 0d0
  gas_emitex = 0d0

  gas_esurf = 0.25*tsp_dt*pc_c*pc_acoef*(4.0*pc_pi*gas_rarr(1)**2)*gas_tempb(1)**4
  if(gas_isvelocity) then
     gas_etot = gas_esurf*tsp_texp**2
  else
     gas_etot = gas_esurf
  endif

  ! Calculating gas_exsource from analytic distribution
  if(gas_srctype/='none') then
     call analytic_source
  else
     gas_exsource=0d0
  endif

  ! Calculating fictitious emission energy per cell: loop
  do ir = 1, gas_nr
     
     gas_emit(ir) =  tsp_dt*gas_fcoef(ir)*gas_siggrey(ir)*pc_c* &
          gas_vals2(ir)%ur*gas_vals2(ir)%vol
     
     if(.not.gas_novolsrc .and. gas_srctype=='none') then
        gas_emit(ir) = gas_emit(ir) + gas_vals2(ir)%nisource
     endif

     gas_etot = gas_etot + gas_emit(ir)

  enddo
  
  ! Adding external energy to gas_etot
  do ir = 1, gas_nr
     do ig = 1, gas_ng
        
        gas_etot = gas_etot+tsp_dt*gas_exsource(ig,ir)* &
             gas_vals2(ir)%vol
        
     enddo
  enddo
  
  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf = nint(gas_esurf*prt_ns/gas_etot)+1
  prt_nnew = prt_nsurf

  ! Calculating number of particles per cell (gas_vals2%nvol): loop
  prt_nexsrc=0
  do ir = 1, gas_nr
     gas_nvol(ir)=nint(gas_emit(ir)*prt_ns/gas_etot)+1
     prt_nnew = prt_nnew + gas_nvol(ir)
     !external source volume numbers
     exsumg = 0d0
     
     do ig=1,gas_ng
        exsumg=exsumg+tsp_dt*gas_exsource(ig,ir)* &
             gas_vals2(ir)%vol
     enddo
     gas_emitex(ir) = exsumg
     gas_nvolex(ir)=nint(exsumg*prt_ns/gas_etot)
     prt_nexsrc = prt_nexsrc + gas_nvolex(ir)
     prt_nnew = prt_nnew + gas_nvolex(ir)
  enddo
  !write(*,*) gas_vals2(1)%nvolex

end subroutine sourcenumbers
