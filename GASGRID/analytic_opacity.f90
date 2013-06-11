subroutine analytic_opacity

  use gasgridmod
  use physconstmod
  !use timestepmod
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
  real*8 :: sigll, sigrr    !dummy variables
  real*8 :: x1, x2  !unitless energy group bounds
  real*8 :: specint !debye type function integrator

  gas_siggrey = 0d0
  gas_cap = 0.
!-- opacity of zero should contribute infinitely little to the rosseland
!   opacity sum (calculated in convert_cap2capros).
  gas_caprosl = huge(gas_cap)
  gas_caprosr = huge(gas_cap)

  gas_sig = 0d0
  gas_sigbl = 0d0
  gas_sigbr = 0d0

  !Calculating grouped Planck and Rosseland opacities
  if(gas_opacanaltype=='none') then
     return
  elseif(gas_opacanaltype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)!{{{
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_siggrey(ir) = gas_sigcoef*gas_vals2(ir)%temp**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_siggrey(ir)*(gas_tempb(ir)/gas_vals2(ir)%temp)**gas_sigtpwr
        sigrr = gas_siggrey(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%temp)**gas_sigtpwr
        do ig = 1, gas_ng
           gas_cap(ig,ir) = gas_siggrey(ir)
           gas_caprosl(ig,ir) = sigll
           gas_caprosr(ig,ir) = sigrr
        enddo
     enddo!}}}
  elseif(gas_opacanaltype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)!{{{
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to integral_g(1/nu^3)
     !
     do ir = 1, gas_nr
        gas_siggrey(ir) = gas_sigcoef*gas_vals2(ir)%temp**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_siggrey(ir)*(gas_tempb(ir)/gas_vals2(ir)%temp)**gas_sigtpwr
        sigrr = gas_siggrey(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%temp)**gas_sigtpwr
        do ig = 1, gas_ng
           !
           !group (Planck) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_cap(ig,ir) = 0.5d0*gas_siggrey(ir)*(x1+x2)/(x1*x2)**2
           !
           !group left (Rosseland) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_caprosl(ig,ir) = 0.5d0*sigll*(x1+x2)/(x1*x2)**2
           !
           !group right (Rosseland) opacities:
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3)
           gas_caprosr(ig,ir) = 0.5d0*sigrr*(x1+x2)/(x1*x2)**2
           !
        enddo
        x1 = (pc_h*pc_c/(pc_ev*gas_wl(gas_ng+1)))/(1d3)
        x2 = (pc_h*pc_c/(pc_ev*gas_wl(1)))/(1d3)
        !gas_siggrey(ir) = 0.5d0*gas_siggrey(ir)*gas_vals2(ir)%tempkev**(-3)*(x1+x2)/(x1*x2)**2
        gas_siggrey(ir) = 0d0
        do ig = 1, gas_ng
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           gas_siggrey(ir) = gas_siggrey(ir)+15d0*gas_cap(ig,ir)*specint(x1,x2,3)/pc_pi**4
        enddo
     enddo!}}}
  elseif(gas_opacanaltype=='pick') then
     ! sigmaP = sigmaR = constant = A (gas_sigcoef set in input.par)!{{{
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     do ir = 1, gas_nr
        gas_siggrey(ir) = gas_sigcoef*gas_vals2(ir)%temp**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
     enddo
     if(gas_suol=='tsta') then    !Case: A
        do ir = 1, gas_nr
           gas_cap(1,ir) = gas_siggrey(ir)
           gas_cap(2,ir) = gas_siggrey(ir)
           gas_caprosl(1,ir) = sigll
           gas_caprosl(2,ir) = sigll
           gas_caprosr(1,ir) = sigrr
           gas_caprosr(2,ir) = sigrr
           do ig = 3, gas_ng
              gas_cap(ig,ir) = 0.
              gas_caprosl(ig,ir) = 0d0
              gas_caprosr(ig,ir) = 0d0
           enddo
        enddo
     elseif(gas_suol=='tstb') then  !Case: B
        do ir = 1, gas_nr
           gas_cap(1,ir) = 2d0*gas_siggrey(ir)/11d0
           gas_cap(2,ir) = 20d0*gas_siggrey(ir)/11d0
           gas_caprosl(1,ir) = 2d0*sigll/11d0
           gas_caprosl(2,ir) = 20d0*sigll/11d0
           gas_caprosr(1,ir) = 2d0*sigrr/11d0
           gas_caprosr(2,ir) = 20d0*sigrr/11d0
           do ig = 3, gas_ng
              gas_cap(ig,ir) = 0.
              gas_caprosl(ig,ir) = 0d0
              gas_caprosr(ig,ir) = 0d0
           enddo
        enddo
     elseif(gas_suol=='tstc') then  !Case: C
        do ir = 1, gas_nr
           gas_cap(1,ir) = 2d0*gas_siggrey(ir)/101d0
           gas_cap(2,ir) = 200d0*gas_siggrey(ir)/101d0
           gas_caprosl(1,ir) = 2d0*sigll/101d0
           gas_caprosl(2,ir) = 200d0*sigll/101d0
           gas_caprosr(1,ir) = 2d0*sigrr/101d0
           gas_caprosr(2,ir) = 200d0*sigrr/101d0
           do ig = 3, gas_ng
              gas_cap(ig,ir) = 0.
              gas_caprosl(ig,ir) = 0d0
              gas_caprosr(ig,ir) = 0d0
           enddo
        enddo
     else
        stop 'analytic_opacity: gas_suol invalid'
     endif!}}}
  elseif(gas_opacanaltype=='line') then
     ! Highly structured line test: group opacities alternate in magnitude!{{{
     ! sigmaP = A*T^B*rho^C
     ! sigmaP_g = sigmaP*func_P(g), sigmaR_g = sigmaP
     do ir = 1, gas_nr
        gas_siggrey(ir) = gas_sigcoef*gas_vals2(ir)%temp**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        !sigll = gas_sigcoef*gas_tempb(ir)**gas_sigtpwr*gas_rhob(ir)**gas_sigrpwr
        !sigrr = gas_sigcoef*gas_tempb(ir+1)**gas_sigtpwr*gas_rhob(ir+1)**gas_sigrpwr
        sigll = gas_siggrey(ir)*(gas_tempb(ir)/gas_vals2(ir)%temp)**gas_sigtpwr
        sigrr = gas_siggrey(ir)*(gas_tempb(ir+1)/gas_vals2(ir)%temp)**gas_sigtpwr
        !
        !set odd group magnitudes (low)
        do ig = 1, gas_ng, 2
           gas_cap(ig,ir) = gas_siggrey(ir)*gas_ldisp1
           gas_caprosl(ig,ir) = sigll*gas_ldisp1
           gas_caprosr(ig,ir) = sigrr*gas_ldisp1
        enddo
        !set even group magnitudes (high)
        do ig = 2, gas_ng, 2
           gas_cap(ig,ir) = gas_siggrey(ir)*gas_ldisp2
           gas_caprosl(ig,ir) = sigll*gas_ldisp2
           gas_caprosr(ig,ir) = sigrr*gas_ldisp2
        enddo
        !
        !calculate Planck, Rosseland opacities
        gas_siggrey(ir) = 0d0
        do ig = 1, gas_ng
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           gas_siggrey(ir) = gas_siggrey(ir)+15d0*gas_cap(ig,ir)* &
                specint(x1,x2,3)/pc_pi**4
        enddo
        !
     enddo!}}}
  else
    stop 'analytic_opacity: gas_opacanaltype invalid'
  endif

  !Calculating grey scattering opacity
  do ir = 1, gas_nr
     gas_sig(ir) = gas_sigcoefs*gas_vals2(ir)%temp**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
     gas_sigbl(ir) = gas_sigcoefs*gas_tempb(ir)**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
     gas_sigbr(ir) = gas_sigcoefs*gas_tempb(ir+1)**gas_sigtpwrs*gas_vals2(ir)%rho**gas_sigrpwrs
  enddo

end subroutine analytic_opacity
