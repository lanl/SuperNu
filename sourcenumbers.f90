SUBROUTINE sourcenumbers

  USE gasgridmod
  USE timestepmod
  USE particlemod
  USE physconstmod
  USE inputparmod
  IMPLICIT NONE

  INTEGER(iknd) :: ir
  REAL(rknd) :: sou
  ! prt_esurf for any new prt_particles from a surface source
  ! prt_nsurf = number of surface prt_particles
  ! prt_nnew = total number of new prt_particles~=prt_ns
  prt_esurf = 0.0_rknd

  !prt_esurf = 0.25*tsp_dt*pc_c*a_coef*(4.0*pc_pi*gas_rarr(1)**2)*gas_tempb(1)**4
  prt_etot = prt_esurf

  ! External source
  DO ir = 1, 39
     gas_nisource(ir) = (gas_rhoarr(ir)/56.0)*pc_navo*gas_nidecay
  ENDDO
  !
  DO ir = 40, gas_nr!9*gas_nr/10, gas_nr
     gas_nisource(ir) = 0.0_rknd
  ENDDO
  !ELSE
  !   DO ir = 1, gas_nr
  !      gas_nisource(ir) = 0.0_rknd
  !   ENDDO
  !ENDIF
  !WRITE(*,*) gas_sigmap(1), gas_sigmapg(1,1), gas_sigmapg(2,1)

  DO ir = 1, gas_nr
     !gas_emit(ir) = (1.e35)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)
     gas_emit(ir) =  tsp_dt*gas_fcoef(ir)*gas_sigmap(ir)*pc_c*gas_ur(ir)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)
     !gas_emit(ir) = gas_emit(ir)*(gas_velno*1.0+gas_velyes*tsp_texp**3)
     sou = gas_nisource(ir)*(4.0*pc_pi*gas_dr3arr(ir)/3.0)*(gas_velno*1.0+gas_velyes*tsp_texp**3)*tsp_dt
     gas_emit(ir) = gas_emit(ir)+sou
     prt_etot = prt_etot+gas_emit(ir)
  ENDDO
  !WRITE(*,*) gas_nisource(1), gas_nisource(2)

  prt_nsurf=NINT(prt_esurf*prt_ns/prt_etot)+1
  prt_nnew = prt_nsurf
  
  DO ir = 1, gas_nr
     gas_nvol(ir)=NINT(gas_emit(ir)*prt_ns/prt_etot)+1
     prt_nnew = prt_nnew+gas_nvol(ir)
  ENDDO
  
END SUBROUTINE sourcenumbers
