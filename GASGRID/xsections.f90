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
  real*8 :: Um, beta, tt, gg, ggg, eps, bb
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

  !Calculating power law heat capacity
  do ir = 1, gas_nr
     gas_vals2(ir)%bcoef=gas_cvcoef*gas_vals2(ir)%tempkev**gas_cvtpwr*gas_vals2(ir)%rho**gas_cvrpwr
  enddo

  !Calculating simple physical group/grey opacities: Planck and Rosseland 
  !Ryan W.:(moved from supernu.f90 in rev 105)
  call analytic_opacity

  !Calculating grey Planck and gouped Planck opacities: loop
  do ir = 1, gas_nr
     gas_sigmapg(1,ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     gas_sigmapg(2,ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     do ig = 3, gas_ng
        gas_sigmapg(ig,ir) = 1.0
     enddo
     gas_sigmap(ir)=0.0
     do ig = 1, gas_ng
        gas_sigmap(ir) = gas_sigmap(ir)+gas_ppick(ig)*gas_sigmapg(ig,ir)
     enddo
     Um = gas_vals2(ir)%bcoef*gas_vals2(ir)%tempkev
     beta = 4.0*gas_vals2(ir)%ur/Um
     gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_sigmap(ir))
     do ig = 1, gas_ng
        gas_emitprobg(ig,ir) = gas_ppick(ig)*gas_sigmapg(ig,ir)/gas_sigmap(ir)
     enddo
  enddo

  !Calculating group Rosseland opacities: loop
  do ir = 1, gas_nr
     gas_sigmargleft(1,ir) = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     gas_sigmargleft(2,ir) = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     do ig = 3, gas_ng
        gas_sigmargleft(ig,ir) = 1.0
     enddo
     gas_sigmargright(1,ir) = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     gas_sigmargright(2,ir) = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     do ig = 3, gas_ng
        gas_sigmargright(ig,ir) = 1.0
     enddo
  enddo

  !Calculating IMC-to-DDMC leakage albedo coefficients (Densmore, 2007): loop
  !These quantities may not need to be stored directly (pending further analysis)
  do ir = 1, gas_nr
     gg = (3.0*gas_fcoef(ir))**0.5
     eps = (4.0/3.0)*gg/(1.0+0.7104*gg)
     do ig = 1, gas_ng
        !Calculating for leakage from left
        !tt = gas_sigmargleft(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        tt = gas_sigmapg(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        gas_ppl(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
        !Calculating for leakage from right
        !tt = gas_sigmargright(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        tt = gas_sigmapg(ig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)
        ggg = (gg*tt)**2
        bb = (3.0/4.0)*gas_fcoef(ir)*tt**2+(ggg+(ggg**2)/4.0)**0.5
        gas_ppr(ig,ir) = 0.5*eps*bb/(bb-(3.0/4.0)*eps*tt)
     enddo
  enddo

  !Calculating DDMC(-to-IMC) leakage opacities (Densmore, 2007, 2012): loop
  do ir = 1, gas_nr
     do ig = 1, gas_ng
        !Computing left-leakage opacities
        if (ir==1) then
           !gas_sigmal(ig,ir)=0.5*gas_ppl(ig,ir)/gas_drarr(ir)
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp))
        elseif(gas_sigmapg(ig,ir-1)*gas_drarr(ir-1)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
           !gas_sigmal(ig,ir)=0.5*gas_ppl(ig,ir)/gas_drarr(ir)
           gas_sigmal(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_sigmal(ig,ir)=gas_sigmal(ig,ir)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp))
        else
           tt = gas_sigmargleft(ig,ir)*gas_drarr(ir)+gas_sigmargright(ig,ir-1)*gas_drarr(ir-1)
           !gas_sigmal(ig,ir) = 2.0/(3.0*gas_drarr(ir))
           gas_sigmal(ig,ir) = (2.0*gas_rarr(ir)**2)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmal(ig,ir) = gas_sigmal(ig,ir)/tt
        endif
        !Computing right-leakage opacities
        if (ir==gas_nr) then
           !gas_sigmar(ig,ir)=0.5*gas_ppr(ig,ir)/gas_drarr(ir)
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp))
        elseif(gas_sigmapg(ig,ir+1)*gas_drarr(ir+1)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) then
           !gas_sigmar(ig,ir)=0.5*gas_ppr(ig,ir)/gas_drarr(ir)
           gas_sigmar(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_sigmar(ig,ir)=gas_sigmar(ig,ir)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp))
        else
           tt = gas_sigmargright(ig,ir)*gas_drarr(ir)+gas_sigmargleft(ig,ir+1)*gas_drarr(ir+1)
           !gas_sigmar(ig,ir) = 2.0/(3.0*gas_drarr(ir))
           gas_sigmar(ig,ir) = (2.0*gas_rarr(ir+1)**2)/(gas_vals2(ir)%dr3_34pi*(gas_velno*1.0+gas_velyes*tsp_texp**2))
           gas_sigmar(ig,ir) = gas_sigmar(ig,ir)/tt
        endif
     enddo
  enddo

end subroutine xsections
