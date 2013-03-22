subroutine analytic_opacity

  use gasgridmod
  use physconstmod
  implicit none

!#####################################
  !This subroutine computes Planck and Rosseland
  !opacities for several simple test group structures.
  !The power law coefficients are used in any structure
  !selection to either fully or partially determine
  !opacity dependence on temperature and density.
!#####################################

  integer :: ir, ig
  real*8 :: fgren, fgrren, fglren !renormalization factors
  real*8 :: sigll, sigrr    !dummy variables
  real*8 :: x1, x2  !unitless energy group bounds
  real*8 :: specint !debye type function integrator

  if(gas_grptype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        sigll = gas_sigmap(ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        sigrr = gas_sigmap(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        do ig = 1, gas_ng
           gas_sigmapg(ig,ir) = gas_sigmap(ir)
           gas_sigmargleft(ig,ir) = sigll
           gas_sigmargright(ig,ir) = sigrr
        enddo
     enddo
  elseif(gas_grptype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to integral_g(1/nu^3)
     !
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        sigll = gas_sigmap(ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        sigrr = gas_sigmap(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        do ig = 1, gas_ng
           !
           !group Planck opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           gas_sigmapg(ig,ir) = gas_sigmap(ir)*log((1.0-exp(-x2))/(1.0-exp(-x1)))/specint(x1,x2,3)
           !
           !group left Rosseland opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir))
           gas_sigmargleft(ig,ir) = sigll*specint(x1,x2,3)/(specint(x1,x2,6)*gas_tempb(ir)**3)
           !
           !group right Rosseland opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir+1))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir+1))
           gas_sigmargright(ig,ir) = sigrr*specint(x1,x2,3)/(specint(x1,x2,6)*gas_tempb(ir+1)**3)
           !
        enddo
        x1 = (pc_h*pc_c/(pc_ev*gas_wl(gas_ng+1)))/(1d3*gas_vals2(ir)%tempkev)
        x2 = (pc_h*pc_c/(pc_ev*gas_wl(1)))/(1d3*gas_vals2(ir)%tempkev)
        gas_sigmap(ir) = gas_sigmap(ir)*log((1.0-exp(-x2))/(1.0-exp(-x1)))/specint(x1,x2,3)
     enddo
  elseif(gas_grptype=='pick') then
     ! sigmaP = sigmaR = constant = A (gas_sigcoef set in input.par)
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     enddo
     if(gas_suol=='tsta') then    !Case: A
        do ir = 1, gas_nr
           gas_sigmapg(1,ir) = gas_sigmap(ir) 
           gas_sigmapg(2,ir) = gas_sigmap(ir)
           gas_sigmargleft(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargleft(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           do ig = 3, gas_ng
              gas_sigmapg(ig,ir) = 0d0
              gas_sigmargleft(ig,ir) = 0d0
              gas_sigmargright(ig,ir) = 0d0
           enddo
        enddo
     elseif(gas_suol=='tstb') then  !Case: B
        do ir = 1, gas_nr
           gas_sigmapg(1,ir) = 2d0*gas_sigmap(ir)/11d0
           gas_sigmapg(2,ir) = 20d0*gas_sigmap(ir)/11d0
           gas_sigmargleft(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargleft(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           do ig = 3, gas_ng
              gas_sigmapg(ig,ir) = 0d0
              gas_sigmargleft(ig,ir) = 0d0
              gas_sigmargright(ig,ir) = 0d0
           enddo
        enddo
     elseif(gas_suol=='tstc') then  !Case: C
        do ir = 1, gas_nr
           gas_sigmapg(1,ir) = 2d0*gas_sigmap(ir)/101d0
           gas_sigmapg(2,ir) = 200d0*gas_sigmap(ir)/101d0
           gas_sigmargleft(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargleft(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(1,ir) = gas_sigmapg(1,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           gas_sigmargright(2,ir) = gas_sigmapg(2,ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
           do ig = 3, gas_ng
              gas_sigmapg(ig,ir) = 0d0
              gas_sigmargleft(ig,ir) = 0d0
              gas_sigmargright(ig,ir) = 0d0
           enddo
        enddo
     else
        stop 'analytic_opacity: gas_suol invalid'
     endif
  elseif(gas_grptype=='line') then
     ! Highly structured line test: group opacities alternate 10^7 in magnitude
     ! sigmaP = A*T^B*rho^C
     ! sigmaP_g = sigmaP*func_P(g), sigmaR_g=sigmaP
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        sigll = gas_sigmap(ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        sigrr = gas_sigmap(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        !
        fgren = 0d0
        fglren= 0d0
        fgrren= 0d0
        do ig = 1, gas_ng, 2
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           fgren = fgren+1d-3*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir))
           fglren = fglren+1d-3*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir+1))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir+1))
           fgrren = fgrren+1d-3*15d0*specint(x1,x2,3)/pc_pi**4
        enddo
        do ig = 2, gas_ng, 2
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           fgren = fgren+1d4*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir))
           fglren = fglren+1d4*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir+1))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir+1))
           fgrren = fgrren+1d4*15d0*specint(x1,x2,3)/pc_pi**4
        enddo
        !set odd group magnitudes (low)
        do ig = 1, gas_ng, 2
           gas_sigmapg(ig,ir) = gas_sigmap(ir)*1d-3/fgren
           gas_sigmargleft(ig,ir) = sigll*1d-3/fglren
           gas_sigmargright(ig,ir) = sigrr*1d-3/fgrren
        enddo
        !set even group magnitudes (high)
        do ig = 2, gas_ng, 2
           gas_sigmapg(ig,ir) = gas_sigmap(ir)*1d4/fgren
           gas_sigmargleft(ig,ir) = sigll*1d4/fglren
           gas_sigmargright(ig,ir) = sigrr*1d4/fgrren
        enddo
     enddo
  else
    stop 'analytic_opacity: gas_grptype invalid'
  endif
  

end subroutine analytic_opacity
