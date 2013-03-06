SUBROUTINE boundary_source

  USE particlemod
  USE timestepmod
  USE physconstmod
  USE gasgridmod
  USE inputparmod
  IMPLICIT NONE

  INTEGER(iknd) :: ipart, ivac, ig, z0
  REAL(rknd) :: r1, r2, P, mu0, r0, Esurfpart

  Esurfpart = Esurf/REAL(Nsurf,rknd)
  Einp = Einp+Esurf

  DO ipart = 1, Nsurf
     ivac = vacantarr(ipart)
     !Picket fence group sampling
     r1 = RAND()
     IF (r1 <= Ppick(1)) THEN
        particles(ivac)%gsrc = 1
     ELSE
        particles(ivac)%gsrc = 2
     ENDIF
     ig = particles(ivac)%gsrc
     
     r1 = RAND()
     r2 = RAND()
     particles(ivac)%musrc = 1.0*MAX(r1,r2)
     IF (ABS(particles(ivac)%musrc)<0.0000001) THEN
        particles(ivac)%musrc = 0.0000001
     ENDIF
     mu0 = particles(ivac)%musrc
     P = PPL(ig,1)*(1.0+1.5*particles(ivac)%musrc)

     r1 = RAND()
     particles(ivac)%tsrc = time+r1*dt

     particles(ivac)%zsrc = 1
     z0 = particles(ivac)%zsrc

     particles(ivac)%rsrc = rarr(1)
     r0 = particles(ivac)%rsrc

     IF ((sigmapg(ig,z0)*drarr(z0)*(velno*1.0+velyes*texp)<5.0_rknd).OR.(in_puretran.EQV..TRUE.)) THEN
        !transport => lab frame quantities
        particles(ivac)%Esrc = Esurfpart*(1.0+velyes*r0*mu0/lspeed)
        particles(ivac)%Ebirth = Esurfpart*(1.0+velyes*r0*mu0/lspeed)
        particles(ivac)%musrc = (mu0+velyes*r0/lspeed)/(1.0+velyes*r0*mu0/lspeed)
        particles(ivac)%rtsrc = 1
     ELSE
        !diffusion => comoving frame quantities (with diffuse reflection accounted)
        particles(ivac)%Esrc = P*Esurfpart
        particles(ivac)%Ebirth = P*Esurfpart
        particles(ivac)%rtsrc = 2
     ENDIF
     
     particles(ivac)%isvacant = .FALSE.

  ENDDO
  !DEALLOCATE(vacantarr)

END SUBROUTINE boundary_source
