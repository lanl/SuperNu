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
  real*8 :: x1, x2  !unitless energy group bounds
  real*8 :: specint !debye type function integrator

  if(gas_grptype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        do ig = 1, gas_ng
           gas_sigmapg(ig,ir) = gas_sigmap(ir)
           gas_sigmargleft(ig,ir) = gas_sigmap(ir)
           gas_sigmargright(ig,ir) = gas_sigmap(ir)
        enddo
     enddo
  elseif(gas_grptype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to integral_g(1/nu^3)
     do ir = 1, gas_nr
        gas_sigmap(ir) = gas_sigcoef*gas_vals2(ir)%tempkev**gas_sigtpwr*gas_vals2(ir)%rho**gas_sigrpwr
        do ig = 1, gas_ng
           x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
           x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
           gas_sigmapg(ig,ir) = gas_sigmap(ir)*log((1.0-exp(-x2))/(1.0-exp(-x1)))/specint(x1,x2,3)
        enddo
        x1 = (pc_h*pc_c/(pc_ev*gas_wl(gas_ng+1)))/(1d3*gas_vals2(ir)%tempkev)
        x2 = (pc_h*pc_c/(pc_ev*gas_wl(1)))/(1d3*gas_vals2(ir)%tempkev)
        gas_sigmap(ir) = gas_sigmap(ir)*log((1.0-exp(-x2))/(1.0-exp(-x1)))/specint(x1,x2,3)
     enddo
  elseif(gas_grptype=='pick') then
  
  elseif(gas_grptype=='line') then

  else

  endif
  

end subroutine analytic_opacity
