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

  integer :: ir
  real*8 :: sou
  ! gas_esurf for any new prt_particles from a surface source
  ! prt_nsurf = number of surface prt_particles
  ! prt_nnew = total number of new prt_particles~=prt_ns
  gas_esurf = 0.0d0

  !gas_esurf = 0.25*tsp_dt*pc_c*a_coef*(4.0*pc_pi*gas_rarr(1)**2)*gas_tempb(1)**4
  gas_etot = gas_esurf

  ! Volume external radiation source loops (gamma ray energy)
  do ir = 1, 39
     gas_vals2(ir)%nisource = (gas_vals2(ir)%rho/56.0)*pc_navo*gas_nidecay
  enddo
  !
  do ir = 40, gas_nr!9*gas_nr/10, gas_nr
     gas_vals2(ir)%nisource = 0.0d0
  enddo
  

  ! Calculating fictitious emission energy per cell: loop
  do ir = 1, gas_nr
     !gas_vals2(ir)%emit = (1.e35)*(4.0*pc_pi*gas_vals2(ir)%dr3_34pi/3.0)
     gas_vals2(ir)%emit =  tsp_dt*gas_fcoef(ir)*gas_sigmap(ir)*pc_c*gas_vals2(ir)%ur*(4.0*pc_pi*gas_vals2(ir)%dr3_34pi/3.0)
     !gas_vals2(ir)%emit = gas_vals2(ir)%emit*(gas_velno*1.0+gas_velyes*tsp_texp**3)
     sou = gas_vals2(ir)%nisource*(4.0*pc_pi*gas_vals2(ir)%dr3_34pi/3.0)*(gas_velno*1.0+gas_velyes*tsp_texp**3)*tsp_dt
     gas_vals2(ir)%emit = gas_vals2(ir)%emit+sou
     gas_etot = gas_etot+gas_vals2(ir)%emit
  enddo
  !write(*,*) gas_nisource(1), gas_nisource(2)

  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf=nint(gas_esurf*prt_ns/gas_etot)+1
  prt_nnew = prt_nsurf
  
  ! Calculating number of particles per cell (gas_vals2%nvol): loop 
  do ir = 1, gas_nr
     gas_vals2(ir)%nvol=nint(gas_vals2(ir)%emit*prt_ns/gas_etot)+1
     prt_nnew = prt_nnew+gas_vals2(ir)%nvol
  enddo
  
end subroutine sourcenumbers
