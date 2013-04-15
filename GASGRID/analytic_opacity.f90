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
  !
  !Since revision 121, calculates grey scattering opacity, gas_sig
!#####################################

  integer :: ir, ig
  real*8 :: fgren, fgrren, fglren !renormalization factors
  real*8 :: sigll, sigrr    !dummy variables
  real*8 :: x1, x2  !unitless energy group bounds
  real*8 :: specint !debye type function integrator

  !Calculating grey scattering opacity
  do ir = 1, gas_nr
     gas_sig(ir) = gas_sigcoefs*gas_vals2(ir)%tempkev**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
     gas_sigbl(ir) = gas_sigcoefs*gas_tempb(ir)**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
     gas_sigbr(ir) = gas_sigcoefs*gas_tempb(ir+1)**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
  enddo
       
  !Calculating grouped Planck and Rosseland opacities
  if(gas_grptype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
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
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_sigmap(ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        sigrr = gas_sigmap(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        do ig = 1, gas_ng
           !
           !group (Planck) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_sigmapg(ig,ir) = 0.5d0*gas_sigmap(ir)*(x1+x2)/(x1*x2)**2
           !
           !group left (Rosseland) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_sigmargleft(ig,ir) = 0.5d0*sigll*(x1+x2)/(x1*x2)**2
           !
           !group right (Rosseland) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_sigmargright(ig,ir) = 0.5d0*sigrr*(x1+x2)/(x1*x2)**2
           !
        enddo
        x1 = (pc_h*pc_c/(pc_ev*gas_wl(gas_ng+1)))/(1d3)
        x2 = (pc_h*pc_c/(pc_ev*gas_wl(1)))/(1d3)
        !gas_sigmap(ir) = 0.5d0*gas_sigmap(ir)*gas_vals2(ir)%tempkev**(-3)*(x1+x2)/(x1*x2)**2
        gas_sigmap(ir)=0d0
        do ig = 1, gas_ng
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           gas_sigmap(ir)=gas_sigmap(ir)+15d0*gas_sigmapg(ig,ir)*specint(x1,x2,3)/pc_pi**4
        enddo
     enddo
  elseif(gas_grptype=='pick') then
     ! sigmaP = sigmaR = constant = A (gas_sigcoef set in input.par)
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     enddo
     if(gas_suol=='tsta') then    !Case: A
        do ir = 1, gas_nr
           gas_sigmapg(1,ir) = gas_sigmap(ir) 
           gas_sigmapg(2,ir) = gas_sigmap(ir)
           gas_sigmargleft(1,ir) = sigll
           gas_sigmargleft(2,ir) = sigll
           gas_sigmargright(1,ir) = sigrr
           gas_sigmargright(2,ir) = sigrr
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
           gas_sigmargleft(1,ir) = 2d0*sigll/11d0
           gas_sigmargleft(2,ir) = 20d0*sigll/11d0
           gas_sigmargright(1,ir) = 2d0*sigrr/11d0
           gas_sigmargright(2,ir) = 20d0*sigrr/11d0
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
           gas_sigmargleft(1,ir) = 2d0*sigll/101d0
           gas_sigmargleft(2,ir) = 200d0*sigll/101d0
           gas_sigmargright(1,ir) = 2d0*sigrr/101d0
           gas_sigmargright(2,ir) = 200d0*sigrr/101d0
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
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_sigmap(ir)*(gas_tempb(ir)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        sigrr = gas_sigmap(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%tempkev)**gas_sigtpwr
        !
        fgren = 0d0
        fglren= 0d0
        fgrren= 0d0
        do ig = 1, gas_ng, 2
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           fgren = fgren+15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir))
           fglren = fglren+15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir+1))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir+1))
           fgrren = fgrren+15d0*specint(x1,x2,3)/pc_pi**4
        enddo
        do ig = 2, gas_ng, 2
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           fgren = fgren+gas_ldisp*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir))
           fglren = fglren+gas_ldisp*15d0*specint(x1,x2,3)/pc_pi**4
           !
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_tempb(ir+1))
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_tempb(ir+1))
           fgrren = fgrren+gas_ldisp*15d0*specint(x1,x2,3)/pc_pi**4
        enddo
        !set odd group magnitudes (low)
        do ig = 1, gas_ng, 2
           gas_sigmapg(ig,ir) = gas_sigmap(ir)/fgren
           gas_sigmargleft(ig,ir) = sigll/fglren
           gas_sigmargright(ig,ir) = sigrr/fgrren
        enddo
        !set even group magnitudes (high)
        do ig = 2, gas_ng, 2
           gas_sigmapg(ig,ir) = gas_sigmap(ir)*gas_ldisp/fgren
           gas_sigmargleft(ig,ir) = sigll*gas_ldisp/fglren
           gas_sigmargright(ig,ir) = sigrr*gas_ldisp/fgrren
        enddo
     enddo
  else
    stop 'analytic_opacity: gas_grptype invalid'
  endif
  

end subroutine analytic_opacity
