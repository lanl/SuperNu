SUBROUTINE sourcenumbers

  USE data_mod
  IMPLICIT NONE

  INTEGER(iknd) :: ir, nsr
  REAL(rknd) :: sou
  ! Esurf for any new particles from a surface source
  ! Nsurf = number of surface particles
  ! Nnew = total number of new particles~=Ns
  Esurf = 0.0_rknd

  !Esurf = 0.25*dt*lspeed*a_coef*(4.0*pi*rarr(1)**2)*Tempb(1)**4
  Etot = Esurf

  ! External source
  DO ir = 1, 39
     nisource(ir) = (rhoarr(ir)/56.0)*Nav*nidecay
  ENDDO
  !
  DO ir = 40, nr!9*nr/10, nr
     nisource(ir) = 0.0_rknd
  ENDDO
  !ELSE
  !   DO ir = 1, nr
  !      nisource(ir) = 0.0_rknd
  !   ENDDO
  !ENDIF
  !WRITE(*,*) sigmap(1), sigmapg(1,1), sigmapg(2,1)

  DO ir = 1, nr
     !Emit(ir) = (1.e35)*(4.0*pi*dr3arr(ir)/3.0)
     Emit(ir) =  dt*fcoef(ir)*sigmap(ir)*lspeed*Ur(ir)*(4.0*pi*dr3arr(ir)/3.0)
     !Emit(ir) = Emit(ir)*(velno*1.0+velyes*texp**3)
     sou = nisource(ir)*(4.0*pi*dr3arr(ir)/3.0)*(velno*1.0+velyes*texp**3)*dt
     Emit(ir) = Emit(ir)+sou
     Etot = Etot+Emit(ir)
  ENDDO
  !WRITE(*,*) nisource(1), nisource(2)

  Nsurf=NINT(Esurf*Ns/Etot)+1
  Nnew = Nsurf
  
  DO ir = 1, nr
     Nvol(ir)=NINT(Emit(ir)*Ns/Etot)+1
     Nnew = Nnew+Nvol(ir)
  ENDDO
  
END SUBROUTINE sourcenumbers
