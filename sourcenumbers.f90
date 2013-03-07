SUBROUTINE sourcenumbers

  USE gasgridmod
  USE timestepmod
  USE particlemod
  USE physconstmod
  USE inputparmod
  IMPLICIT NONE

!##################################################
  !This subroutine computes the distribution of source particles each
  !time step.  A fraction of the source particle number prt_ns is given
  !to each cell based on the amount of energy emitted by the cell.
!##################################################

  INTEGER :: ir
  REAL*8 :: sou
  ! prt_esurf for any new prt_particles from a surface source
  ! prt_nsurf = number of surface prt_particles
  ! prt_nnew = total number of new prt_particles~=prt_ns
  prt_esurf = 0.0d0

  !prt_esurf = 0.25*tsp_dt*pc_c*a_coef*(4.0*pc_pi*gas_rarr(1)**2)*gas_tempb(1)**4
  prt_etot = prt_esurf

  ! Volume external radiation source loops (gamma ray energy)
  DO ir = 1, 39
     gas_nisource(ir) = (gas_rhoarr(ir)/56.0)*pc_navo*gas_nidecay
  ENDDO
  !
  DO ir = 40, gas_nr!9*gas_nr/10, gas_nr
     gas_nisource(ir) = 0.0d0
  ENDDO
  

  ! Calculating fictitious emission energy per cell: loop
  DO ir = 1, gas_nr
     !gas_emit(ir) = (1.e35)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)
     gas_emit(ir) =  tsp_dt*gas_fcoef(ir)*gas_sigmap(ir)*pc_c*gas_ur(ir)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)
     !gas_emit(ir) = gas_emit(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**3)
     sou = gas_nisource(ir)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)*(gas_velno*1.0+gas_velyes*tsp_texp**3)*tsp_dt
     gas_emit(ir) = gas_emit(ir)+sou
     prt_etot = prt_etot+gas_emit(ir)
  ENDDO
  !WRITE(*,*) gas_nisource(1), gas_nisource(2)

  ! Calculating number of domain inner boundary particles (if any)
  prt_nsurf=NINT(prt_esurf*prt_ns/prt_etot)+1
  prt_nnew = prt_nsurf
  
  ! Calculating number of particles per cell (gas_nvol): loop 
  DO ir = 1, gas_nr
     gas_nvol(ir)=NINT(gas_emit(ir)*prt_ns/prt_etot)+1
     prt_nnew = prt_nnew+gas_nvol(ir)
  ENDDO
  
END SUBROUTINE sourcenumbers
