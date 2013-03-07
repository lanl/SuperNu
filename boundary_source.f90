SUBROUTINE boundary_source

  USE particlemod
  USE timestepmod
  USE physconstmod
  USE gasgridmod
  USE inputparmod
  IMPLICIT NONE

  INTEGER :: ipart, ivac, ig, z0
  REAL*8 :: r1, r2, P, mu0, r0, Esurfpart

  Esurfpart = prt_esurf/REAL(prt_nsurf)
  prt_einp = prt_einp+prt_esurf

  DO ipart = 1, prt_nsurf
     ivac = prt_vacantarr(ipart)
     !Picket fence group sampling
     r1 = RAND()
     IF (r1 <= gas_ppick(1)) THEN
        prt_particles(ivac)%gsrc = 1
     ELSE
        prt_particles(ivac)%gsrc = 2
     ENDIF
     ig = prt_particles(ivac)%gsrc
     
     r1 = RAND()
     r2 = RAND()
     prt_particles(ivac)%musrc = 1.0*MAX(r1,r2)
     IF (ABS(prt_particles(ivac)%musrc)<0.0000001) THEN
        prt_particles(ivac)%musrc = 0.0000001
     ENDIF
     mu0 = prt_particles(ivac)%musrc
     P = gas_ppl(ig,1)*(1.0+1.5*prt_particles(ivac)%musrc)

     r1 = RAND()
     prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt

     prt_particles(ivac)%zsrc = 1
     z0 = prt_particles(ivac)%zsrc

     prt_particles(ivac)%rsrc = gas_rarr(1)
     r0 = prt_particles(ivac)%rsrc

     IF ((gas_sigmapg(ig,z0)*gas_drarr(z0)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0).OR.(in_puretran.EQV..TRUE.)) THEN
        !transport => lab frame quantities
        prt_particles(ivac)%Esrc = Esurfpart*(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%Ebirth = Esurfpart*(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%rtsrc = 1
     ELSE
        !diffusion => comoving frame quantities (with diffuse reflection accounted)
        prt_particles(ivac)%Esrc = P*Esurfpart
        prt_particles(ivac)%Ebirth = P*Esurfpart
        prt_particles(ivac)%rtsrc = 2
     ENDIF
     
     prt_particles(ivac)%isvacant = .FALSE.

  ENDDO
  !DEALLOCATE(prt_vacantarr)

END SUBROUTINE boundary_source
