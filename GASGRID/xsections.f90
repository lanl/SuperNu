subroutine xsections

  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes cross sections (opacities) used in
  !the particle advance phase of the program.  These opacities
  !include the grey Planck, grouped Planck, grouped Rosseland,
  !and DDMC grouped leakage opacities.
!##################################################

  integer :: ir, ig
  real*8 :: Um, beta, tt, gg, ggg, eps, bb, sigtot
  real*8 :: x1, x2, rrcenter
  real*8 :: specint
  ! Here: left=>toward r=0 and right=>outward

  !Interpolating cell boundary temperatures (in keV currently): loop
  if(gas_isshell) then
     gas_tempb(1)=gas_templ0
  else
     gas_tempb(1)=gas_vals2(1)%tempkev
  endif
  !gas_tempb(1) = 1.0
  do ir = 2, gas_nr
     gas_tempb(ir) = (gas_vals2(ir)%tempkev**4+gas_vals2(ir-1)%tempkev**4)/2.0
     gas_tempb(ir) = gas_tempb(ir)**0.25
  enddo
  gas_tempb(gas_nr+1)=gas_vals2(gas_nr)%tempkev
  !Interpolating cell boundary densities (in g/cm^3): loop
  gas_rhob(1)=gas_vals2(1)%rho
  do ir = 2, gas_nr
     !gas_rhob(ir)=(gas_vals2(ir)%rho*gas_vals2(ir)%vol+ &
     !     gas_vals2(ir-1)%rho*gas_vals2(ir-1)%vol)/ &
     !     (gas_vals2(ir)%vol+gas_vals2(ir-1)%vol)
     gas_rhob(ir)=(gas_vals2(ir)%rho*gas_vals2(ir-1)%rho)**0.5d0
  enddo
  gas_rhob(gas_nr+1) = gas_vals2(gas_nr)%rho

  !Calculating power law heat capacity
  do ir = 1, gas_nr
     gas_vals2(ir)%bcoef=gas_cvcoef*gas_vals2(ir)%tempkev**gas_cvtpwr*gas_vals2(ir)%rho**gas_cvrpwr
  enddo

  !Calculating simple physical group/grey opacities: Planck and Rosseland 
  !Ryan W.:(moved from supernu.f90 in rev 105)
  call analytic_opacity
  
  !Calculating Fleck factor: 
  do ir = 1, gas_nr
     Um = gas_vals2(ir)%bcoef*gas_vals2(ir)%tempkev
     beta = 4.0*gas_vals2(ir)%ur/Um
     gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_sigmap(ir))
  enddo

  !Calculating grouped volume emission probabilities:
  if(gas_isanalgrp.and.gas_grptype=='pick') then
     do ir = 1, gas_nr
        gas_emitprobg(1,ir) = gas_ppick(1)*gas_sigmapg(1,ir)/gas_sigmap(ir)
        gas_emitprobg(2,ir) = gas_ppick(2)*gas_sigmapg(2,ir)/gas_sigmap(ir)
        do ig = 3, gas_ng
           gas_emitprobg(ig,ir) = 0d0
        enddo
     enddo
  else
     if(gas_ng==1) then
        do ir = 1, gas_nr
           gas_emitprobg(1,ir) = 1d0
        enddo
     else
        do ir = 1, gas_nr
           do ig = 1, gas_ng
              x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
              x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
              gas_emitprobg(ig,ir) = 15d0*specint(x1,x2,3)*gas_sigmapg(ig,ir)/ &
                   (gas_sigmap(ir)*pc_pi**4)
           enddo
        enddo
     endif
  endif
  write(*,*) gas_fcoef(1), gas_fcoef(2)
  !Calculating IMC-to-DDMC leakage albedo coefficients (Densmore, 2007): loop
  !These quantities may not need to be stored directly (pending further analysis)
  do ir = 1, gas_nr
!-----------------------------------------------------------------
     !Ryan W.: emissivity albedo condition deprecated until curvature
     !conditions are added.
!----------------------------------------------------------------- 
     !do ig = 1, gas_ng
        !Calculating for leakage from left
        !total optical depth      ||   ||
     !   sigtot = gas_sigbl(ir)+gas_sigmargleft(ig,ir)
     !   gg = (3.0*(1d0-gas_sigmargleft(ig,ir)/sigtot))**0.5
     !   eps = (4.0/3.0)*gg/(1.0+0.7104*gg)
        !
     !   tt = (gas_sigmargleft(ig,ir)+gas_sigbl(ir)) &
     !        *gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        !
     !   ggg = (gg*tt)**2
     !   bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
     !   gas_ppl(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
        !Calculating for leakage from right
        !total optical depth      ||   ||
     !   sigtot = gas_sigbr(ir)+gas_sigmargright(ig,ir)
     !   gg = (3.0*(1d0-gas_sigmargright(ig,ir)/sigtot))**0.5
     !   eps = (4.0/3.0)*gg/(1.0+0.7104*gg)
        !
     !   tt = (gas_sigmargright(ig,ir)+gas_sigbr(ir)) &
     !        *gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
     !   ggg = (gg*tt)**2
     !   bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
     !   gas_ppr(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
     !enddo
     do ig = 1, gas_ng
        rrcenter=0.5*(gas_rarr(ir)+gas_rarr(ir+1))
        tt = (gas_sigmargleft(ig,ir)+gas_sigbl(ir)) &
             *gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp) !&
             !*gas_rarr(ir)**2/(rrcenter**2-gas_rarr(ir)**2)
        gas_ppl(ig,ir) = 4.0d0/(3d0*tt+6d0*0.7104d0)
        !
        tt = (gas_sigmargright(ig,ir)+gas_sigbr(ir)) &
             *gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp) !&
             !*gas_rarr(ir+1)**2/(gas_rarr(ir+1)**2-rrcenter**2)
        gas_ppr(ig,ir) = 4.0d0/(3d0*tt+6d0*0.7104d0)
     enddo
  enddo
  
  !Calculating DDMC(-to-IMC) leakage opacities (Densmore, 2007, 2012): loop
  do ir = 1, gas_nr
     do ig = 1, gas_ng
        !Computing left-leakage opacities
        if (ir==1) then
        !
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp))
        !   
        elseif((gas_sig(ir-1)+gas_sigmapg(ig,ir-1))*gas_drarr(ir-1) &
             *(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
        !   
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp))
        !   
        else
        !
           tt = (gas_sigbl(ir)+gas_sigmargleft(ig,ir))*gas_drarr(ir)+ &
                (gas_sigbr(ir-1)+gas_sigmargright(ig,ir-1))*gas_drarr(ir-1)
        !   
           gas_sigmal(ig,ir)=(2.0*gas_rarr(ir)**2)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmal(ig,ir) = gas_sigmal(ig,ir)/tt
        endif
        !Computing right-leakage opacities
        if (ir==gas_nr) then
        !
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp))
        !   
        elseif((gas_sig(ir+1)+gas_sigmapg(ig,ir+1))*gas_drarr(ir+1) &
             *(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
        !   
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp))
        !   
        else
        !
           tt = (gas_sigbr(ir)+gas_sigmargright(ig,ir))*gas_drarr(ir)+ &
                (gas_sigbl(ir+1)+gas_sigmargleft(ig,ir+1))*gas_drarr(ir+1)
           gas_sigmar(ig,ir) = (2.0*gas_rarr(ir+1)**2)/(gas_vals2(ir)%dr3_34pi &
                *(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmar(ig,ir) = gas_sigmar(ig,ir)/tt
        endif
     enddo
  enddo

end subroutine xsections
