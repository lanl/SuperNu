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
  real*8 :: x1, x2  !unitless energy group bounds
  real*8,external :: specint !debye type function integrator

  gas_siggrey = 0d0
  gas_cap = 0d0
  gas_sig = 0d0

  !Calculating grouped Planck and Rosseland opacities
  if(gas_opacanaltype=='none') then
     return
  elseif(gas_opacanaltype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)!{{{
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     do ir = 1, gas_nx
        gas_siggrey(ir,1,1) = gas_sigcoef*gas_temp(ir,1,1)**gas_sigtpwr* &
             gas_vals2(ir,1,1)%rho**gas_sigrpwr
        gas_cap(:,ir,1,1) = gas_siggrey(ir,1,1)
     enddo!}}}
  elseif(gas_opacanaltype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)!{{{
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to
     ! integral_g(1/nu^3)
     !
     do ir = 1, gas_nx
        gas_siggrey(ir,1,1) = gas_sigcoef*gas_temp(ir,1,1)**gas_sigtpwr* &
             gas_vals2(ir,1,1)%rho**gas_sigrpwr
        do ig = 1, gas_ng
           !
           !group opacities:

            x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb)
            x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb)

            gas_cap(ig,ir,1,1) = 0.5d0*gas_siggrey(ir,1,1)*(x1+x2)/(x1*x2)**2
           !
           !
        enddo
        gas_siggrey(ir,1,1) = 0d0
        do ig = 1, gas_ng
           x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_temp(ir,1,1))
           x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(ir,1,1))
           gas_siggrey(ir,1,1) = gas_siggrey(ir,1,1)+15d0*gas_cap(ig,ir,1,1)* &
                specint(x1,x2,3)/pc_pi**4
        enddo
     enddo!}}}
  elseif(gas_opacanaltype=='pick') then
     if(gas_ng/=2) stop 'analytic_opacity: invalid gas_ng'
     ! sigmaP = sigmaR = constant = A (gas_sigcoef set in input.par)!{{{
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     do ir = 1, gas_nx
        gas_siggrey(ir,1,1) = gas_sigcoef*gas_temp(ir,1,1)**gas_sigtpwr* &
             gas_vals2(ir,1,1)%rho**gas_sigrpwr     
        if(gas_suol=='tsta') then    !Case: A
           gas_cap(1,ir,1,1) = gas_siggrey(ir,1,1)
           gas_cap(2,ir,1,1) = gas_siggrey(ir,1,1)
        elseif(gas_suol=='tstb') then  !Case: B
           gas_cap(1,ir,1,1) = 2d0*gas_siggrey(ir,1,1)/11d0
           gas_cap(2,ir,1,1) = 20d0*gas_siggrey(ir,1,1)/11d0
        elseif(gas_suol=='tstc') then  !Case: C
           gas_cap(1,ir,1,1) = 2d0*gas_siggrey(ir,1,1)/101d0
           gas_cap(2,ir,1,1) = 200d0*gas_siggrey(ir,1,1)/101d0
!-- added evacuated picket test
        elseif(gas_suol=='tstd') then !Case: D (not in SuOlson)
           gas_cap(1,ir,1,1) = 0d0
           gas_cap(2,ir,1,1) = 2d0*gas_siggrey(ir,1,1)
        else
           stop 'analytic_opacity: gas_suol invalid'
        endif!}}}
     enddo
  elseif(gas_opacanaltype=='line') then
     ! Highly structured line test: group opacities alternate in magnitude!{{{
     ! sigmaP = A*T^B*rho^C
     ! sigmaP_g = sigmaP*func_P(g), sigmaR_g = sigmaP
     do ir = 1, gas_nx
        gas_siggrey(ir,1,1) = gas_sigcoef*gas_temp(ir,1,1)**gas_sigtpwr * &
           gas_vals2(ir,1,1)%rho**gas_sigrpwr
        !
        !set odd group magnitudes (low)
        do ig = 1, gas_ng, 2
           gas_cap(ig,ir,1,1) = gas_siggrey(ir,1,1)*gas_ldisp1
        enddo
        !set even group magnitudes (high)
        do ig = 2, gas_ng, 2
           gas_cap(ig,ir,1,1) = gas_siggrey(ir,1,1)*gas_ldisp2
        enddo
        !
        !calculate Planck, Rosseland opacities
        gas_siggrey(ir,1,1) = 0d0
        do ig = 1, gas_ng
           x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_temp(ir,1,1))
           x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(ir,1,1))
           gas_siggrey(ir,1,1) = gas_siggrey(ir,1,1)+15d0*gas_cap(ig,ir,1,1)* &
                specint(x1,x2,3)/pc_pi**4
        enddo
        !
     enddo!}}}
  else
    stop 'analytic_opacity: gas_opacanaltype invalid'
  endif

  !Calculating grey scattering opacity
  gas_sig = gas_sigcoefs*gas_temp**gas_sigtpwrs* &
       gas_vals2%rho**gas_sigrpwrs

end subroutine analytic_opacity
