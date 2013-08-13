subroutine boundary_source

  use particlemod
  use timestepmod
  use physconstmod
  use gasgridmod
  use inputparmod
  implicit none

  integer :: ipart, ivac, ig, iig, z0
  real*8 :: r1, r2, P, mu0, r0, Esurfpart, wl0
  real*8 :: denom2, x1, x2, specint, help
  real*8, dimension(gas_ng) :: emitsurfprobg  !surface emission probabilities 
  !, Ryan W.: size will=# of groups in first cell

  Esurfpart = gas_esurf/real(prt_nsurf)

  if(gas_isvelocity) then
     help = tsp_texp
  else
     help = 1d0
  endif
  !Calculating grouped thermal emission probabilities
  if(gas_opacanaltype=='pick') then
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
     !sampling group
     denom2 = 0d0
     r1 = rand()
     do ig = 1, gas_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     !particle group removed (rev 120)
     !prt_particles(ivac)%gsrc = iig
     !Calculating comoving particle wavelength uniformly in group
     if(gas_isvelocity.and.gas_srctype=='manu') then
        iig = 2
     endif
     r1 = rand()
     wl0 = (1d0-r1)*gas_wl(iig)+r1*gas_wl(iig+1)

     r1 = rand()
     r2 = rand()
     prt_particles(ivac)%musrc = 1d0*max(r1,r2)
     if (abs(prt_particles(ivac)%musrc)<0.0000001) then
        prt_particles(ivac)%musrc = 0.0000001
     endif
     mu0 = prt_particles(ivac)%musrc
     P = gas_ppl(iig,1)*(1.0+1.5*prt_particles(ivac)%musrc)
     !write(*,*) P, prt_particles(ivac)%musrc, iig
     r1 = rand()
     prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt

     if(gas_isvelocity.and.gas_srctype=='manu') then
        prt_particles(ivac)%zsrc = gas_nr
        mu0 = -mu0
        prt_particles(ivac)%musrc = mu0
     else
        prt_particles(ivac)%zsrc = 1
     endif
     z0 = prt_particles(ivac)%zsrc
     

     prt_particles(ivac)%rsrc = gas_rarr(1)
     r0 = prt_particles(ivac)%rsrc
     
     if (((gas_sig(z0)+gas_cap(iig,z0))*gas_drarr(z0)* &
          help < prt_tauddmc*gas_curvcent(z0)) &
          .or.(in_puretran.eqv..true.).or.P>1d0.or.P<0d0) then
        if(gas_isvelocity) then
        !transport => lab frame quantities
           prt_particles(ivac)%Esrc = Esurfpart*(1.0+r0*mu0/pc_c)
           prt_particles(ivac)%Ebirth = Esurfpart*(1.0+r0*mu0/pc_c)
        !(rev 120)
           prt_particles(ivac)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
        !
           prt_particles(ivac)%musrc = (mu0+r0/pc_c)/(1.0+r0*mu0/pc_c)
        else
           prt_particles(ivac)%Esrc = Esurfpart
           prt_particles(ivac)%Ebirth = Esurfpart
        !
           prt_particles(ivac)%wlsrc = wl0
        !
           prt_particles(ivac)%musrc = mu0
        endif
        prt_particles(ivac)%rtsrc = 1
     else
        !diffusion => comoving frame quantities (with diffuse reflection accounted)
        prt_particles(ivac)%Esrc = P*Esurfpart
        prt_particles(ivac)%Ebirth = P*Esurfpart
        !(rev 120)
        prt_particles(ivac)%wlsrc = wl0
        !
        prt_particles(ivac)%musrc = mu0
        !
        prt_particles(ivac)%rtsrc = 2
     endif

     prt_particles(ivac)%isvacant = .false.

     !source tally
     !if(prt_particles(ivac)%rtsrc==1) then
     !   gas_eraddens(iig,z0)=gas_eraddens(iig,z0)+Esurfpart
     !else
     !   gas_eraddens(iig,z0)=gas_eraddens(iig,z0)+P*Esurfpart
     !endif
  enddo
  !deallocate(prt_vacantarr)

end subroutine boundary_source
