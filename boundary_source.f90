subroutine boundary_source

  use particlemod
  use timestepmod
  use physconstmod
  use gasgridmod
  use inputparmod
  implicit none

  integer :: ipart, ivac, ig, z0
  real*8 :: r1, r2, P, mu0, r0, Esurfpart

  Esurfpart = gas_esurf/real(prt_nsurf)
  gas_einp = gas_einp+gas_esurf

  do ipart = 1, prt_nsurf
     ivac = prt_vacantarr(ipart)
     !Picket fence group sampling
     r1 = rand()
     if (r1 <= gas_ppick(1)) then
        prt_particles(ivac)%gsrc = 1
     else
        prt_particles(ivac)%gsrc = 2
     endif
     ig = prt_particles(ivac)%gsrc
     
     r1 = rand()
     r2 = rand()
     prt_particles(ivac)%musrc = 1.0*max(r1,r2)
     if (abs(prt_particles(ivac)%musrc)<0.0000001) then
        prt_particles(ivac)%musrc = 0.0000001
     endif
     mu0 = prt_particles(ivac)%musrc
     P = gas_ppl(ig,1)*(1.0+1.5*prt_particles(ivac)%musrc)

     r1 = rand()
     prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt

     prt_particles(ivac)%zsrc = 1
     z0 = prt_particles(ivac)%zsrc

     prt_particles(ivac)%rsrc = gas_rarr(1)
     r0 = prt_particles(ivac)%rsrc

     if ((gas_sigmapg(ig,z0)*gas_drarr(z0)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0).OR.(in_puretran.eqv..true.)) then
        !transport => lab frame quantities
        prt_particles(ivac)%Esrc = Esurfpart*(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%Ebirth = Esurfpart*(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
        prt_particles(ivac)%rtsrc = 1
     else
        !diffusion => comoving frame quantities (with diffuse reflection accounted)
        prt_particles(ivac)%Esrc = P*Esurfpart
        prt_particles(ivac)%Ebirth = P*Esurfpart
        prt_particles(ivac)%rtsrc = 2
     endif
     
     prt_particles(ivac)%isvacant = .false.

  enddo
  !deallocate(prt_vacantarr)

end subroutine boundary_source
