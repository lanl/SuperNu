subroutine boundary_source

  use particlemod
  use timestepmod
  use physconstmod
  use gasgridmod
  use inputparmod
  implicit none

  integer :: ipart, ivac, ig, iig, z0
  real*8 :: r1, r2, P, mu0, r0, Esurfpart
  real*8 :: denom2, x1, x2, specint
  real*8, dimension(gas_ng) :: emitsurfprobg  !surface emission probabilities 
  !, Ryan W.: size will=# of groups in first cell

  Esurfpart = gas_esurf/real(prt_nsurf)

  !Calculating grouped thermal emission probabilities
  if(gas_isanalgrp.and.gas_grptype=='pick') then
     emitsurfprobg(1) = gas_ppick(1)
     emitsurfprobg(2) = gas_ppick(2)
     do ig = 3, gas_ng
        emitsurfprobg(3) = 0d0
     enddo
  else
     do ig = 1, gas_ng
        x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(1))
        x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(1))
        emitsurfprobg(ig) = 15d0*specint(x1,x2,3)/pc_pi**4 
     enddo
  endif
  
  !write(6,*) prt_nsurf
  !Instantiating surface particles:
  do ipart = 1, prt_nsurf
     ivac = prt_vacantarr(ipart)
     !Picket fence group sampling
     denom2 = 0d0
     r1 = rand()
     do ig = 1, gas_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     !if(ipart==1.or.ipart==prt_nsurf) then
     !   write(6,*) iig
     !endif
     prt_particles(ivac)%gsrc = iig

     r1 = rand()
     r2 = rand()
     prt_particles(ivac)%musrc = 1.0*max(r1,r2)
     if (abs(prt_particles(ivac)%musrc)<0.0000001) then
        prt_particles(ivac)%musrc = 0.0000001
     endif
     mu0 = prt_particles(ivac)%musrc
     P = gas_ppl(iig,1)*(1.0+1.5*prt_particles(ivac)%musrc)

     r1 = rand()
     prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt

     prt_particles(ivac)%zsrc = 1
     z0 = prt_particles(ivac)%zsrc

     prt_particles(ivac)%rsrc = gas_rarr(1)
     r0 = prt_particles(ivac)%rsrc

     if ((gas_sigmapg(iig,z0)*gas_drarr(z0)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0).OR.(in_puretran.eqv..true.)) then
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
