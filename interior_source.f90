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
  DO ipart = prt_nsurf+1, prt_nnew
     ivac = prt_vacantarr(ipart)
     isnotvacnt = .FALSE.
     DO WHILE (isnotvacnt.EQV..FALSE.)
        IF (irused(ir)<prt_nvol(ir)) THEN
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0.0
           r1 = RAND()
           DO ig = 1, gas_ng
              iig = ig
              IF (r1>=denom2.AND.r1<denom2+gas_emitprobg(ig,ir)) EXIT
              denom2 = denom2+gas_emitprobg(ig,ir)
           ENDDO
           prt_particles(ivac)%gsrc = iig
           !Calculating radial position
           r1 = 0.0
           r2 = 1.0
           uul = gas_tempb(ir)**4
           uur = gas_tempb(ir+1)**4
           uumax = MAX(uul,uur)
           DO WHILE (r2 > r1)
              r3 = RAND()
              r0 = (r3*gas_rarr(ir+1)**3+(1.0-r3)*gas_rarr(ir)**3)**(1.0/3.0)
              r3 = (r0-gas_rarr(ir))/gas_drarr(ir)
              r1 = (r3*uur+(1.0-r3)*uul)/uumax
              r2 = RAND()
           ENDDO
           prt_particles(ivac)%rsrc = r0
           
           !Calculating direction cosine (comoving)
           r1 = RAND()
           mu0 = 1.0-2.0*r1
           !Calculating particle tsp_time
           r1 = RAND()
           prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt
           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = gas_emit(ir)/REAL(prt_nvol(ir),rknd)
           IF ((gas_sigmapg(iig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0_rknd).OR.(puretran.EQV..TRUE.)) THEN
              prt_particles(ivac)%Esrc = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%Ebirth = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%rtsrc = 1
           ELSE
              prt_particles(ivac)%Esrc = Ep0
              prt_particles(ivac)%Ebirth = Ep0
              prt_particles(ivac)%musrc = mu0
              prt_particles(ivac)%rtsrc = 2
           ENDIF
           !Setting ir = zone of particle
           prt_particles(ivac)%zsrc = ir
           !Setting particle index to not vacant
           prt_particles(ivac)%isvacant = .FALSE.
           
           isnotvacnt = .TRUE.
           
        ELSE
           ir = ir + 1
        ENDIF
     ENDDO
  ENDDO
  
  !DEALLOCATE(prt_vacantarr)

END SUBROUTINE interior_source
