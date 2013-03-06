SUBROUTINE interior_source

  USE gasgridmod
  USE timestepmod
  USE particlemod
  USE physconstmod
  USE inputparmod

  IMPLICIT NONE

  INTEGER(iknd) :: ir, ipart, ivac, ig, iig
  INTEGER(iknd), DIMENSION(gas_nr) :: irused
  REAL(rknd) :: r1, r2, r3, uul, uur, uumax, mu0, r0, Ep0
  REAL(rknd) :: denom2
  LOGICAL :: isnotvacnt

  ir = 1
  irused(1:gas_nr) = 0
  DO ipart = Nsurf+1, Nnew
     ivac = vacantarr(ipart)
     isnotvacnt = .FALSE.
     DO WHILE (isnotvacnt.EQV..FALSE.)
        IF (irused(ir)<Nvol(ir)) THEN
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0.0
           r1 = RAND()
           DO ig = 1, gas_ng
              iig = ig
              IF (r1>=denom2.AND.r1<denom2+EmitProbg(ig,ir)) EXIT
              denom2 = denom2+EmitProbg(ig,ir)
           ENDDO
           particles(ivac)%gsrc = iig
           !Calculating radial position
           r1 = 0.0
           r2 = 1.0
           uul = Tempb(ir)**4
           uur = Tempb(ir+1)**4
           uumax = MAX(uul,uur)
           DO WHILE (r2 > r1)
              r3 = RAND()
              r0 = (r3*rarr(ir+1)**3+(1.0-r3)*rarr(ir)**3)**(1.0/3.0)
              r3 = (r0-rarr(ir))/drarr(ir)
              r1 = (r3*uur+(1.0-r3)*uul)/uumax
              r2 = RAND()
           ENDDO
           particles(ivac)%rsrc = r0
           
           !Calculating direction cosine (comoving)
           r1 = RAND()
           mu0 = 1.0-2.0*r1
           !Calculating particle time
           r1 = RAND()
           particles(ivac)%tsrc = time+r1*dt
           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = Emit(ir)/REAL(Nvol(ir),rknd)
           IF ((sigmapg(iig,ir)*drarr(ir)*(velno*1.0+velyes*texp)<5.0_rknd).OR.(in_puretran.EQV..TRUE.)) THEN
              particles(ivac)%Esrc = Ep0*(1.0+velyes*r0*mu0/pc_c)
              particles(ivac)%Ebirth = Ep0*(1.0+velyes*r0*mu0/pc_c)
              particles(ivac)%musrc = (mu0+velyes*r0/pc_c)/(1.0+velyes*r0*mu0/pc_c)
              particles(ivac)%rtsrc = 1
           ELSE
              particles(ivac)%Esrc = Ep0
              particles(ivac)%Ebirth = Ep0
              particles(ivac)%musrc = mu0
              particles(ivac)%rtsrc = 2
           ENDIF
           !Setting ir = zone of particle
           particles(ivac)%zsrc = ir
           !Setting particle index to not vacant
           particles(ivac)%isvacant = .FALSE.
           
           isnotvacnt = .TRUE.
           
        ELSE
           ir = ir + 1
        ENDIF
     ENDDO
  ENDDO
  
  !DEALLOCATE(vacantarr)

END SUBROUTINE interior_source
