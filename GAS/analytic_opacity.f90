subroutine analytic_opacity

  use inputparmod
  use groupmod
  use gridmod
  use gasmod
  use physconstmod
  use miscmod, only:specint
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

  integer :: i, ig
  real*8 :: x1, x2  !unitless energy group bounds

  gas_capgrey = 0d0
  gas_cap = 0.
  gas_sig = 0d0

  !Calculating grey scattering opacity
  gas_sig = in_gas_sigcoef*gas_temp**in_gas_sigtpwr* &
       gas_rho**in_gas_sigrpwr

  !Calculating grouped Planck and Rosseland opacities
  if(in_opacanaltype=='none') then
     return
  elseif(in_opacanaltype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     gas_capgrey = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        gas_cap(:,i) = sngl(gas_capgrey(i))
     enddo

  elseif(in_opacanaltype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)!{{{
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to
     ! integral_g(1/nu^3)
     !
     gas_capgrey = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        do ig = 1, grp_ng
            x1 = pc_h*pc_c*grp_wlinv(ig+1)/pc_kb
            x2 = pc_h*pc_c*grp_wlinv(ig)/pc_kb
            gas_cap(ig,i) = sngl(0.5d0*gas_capgrey(i)*(x1+x2)/(x1*x2)**2)
        enddo
        gas_capgrey(i) = 0d0
        do ig = 1, grp_ng
           x1 = pc_h*pc_c*grp_wlinv(ig+1)/(pc_kb*gas_temp(i))
           x2 = pc_h*pc_c*grp_wlinv(ig)/(pc_kb*gas_temp(i))
           gas_capgrey(i) = gas_capgrey(i)+15d0*gas_cap(ig,i)* &
                specint(x1,x2,3,10)/pc_pi**4
        enddo
     enddo!}}}
  elseif(in_opacanaltype=='pick') then
     if(grp_ng/=2) stop 'analytic_opacity: invalid grp_ng'
     ! sigmaP = sigmaR = constant = A (dd_sigcoef set in input.par)!{{{
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     if(grd_ny>1) stop 'analytic_opacity: no 2D for opacanaltyp=pick'
     do i = 1, gas_ncell
        gas_capgrey(i) = in_gas_capcoef*gas_temp(i)**in_gas_captpwr* &
             gas_rho(i)**in_gas_caprpwr     
        if(in_suol=='tsta') then    !Case: A
           gas_cap(1,i) = sngl(gas_capgrey(i))
           gas_cap(2,i) = sngl(gas_capgrey(i))
        elseif(in_suol=='tstb') then  !Case: B
           gas_cap(1,i) = sngl(2d0*gas_capgrey(i)/11d0)
           gas_cap(2,i) = sngl(20d0*gas_capgrey(i)/11d0)
        elseif(in_suol=='tstc') then  !Case: C
           gas_cap(1,i) = sngl(2d0*gas_capgrey(i)/101d0)
           gas_cap(2,i) = sngl(200d0*gas_capgrey(i)/101d0)
!-- added evacuated picket test
        elseif(in_suol=='tstd') then !Case: D (not in SuOlson)
           gas_cap(1,i) = sngl(0d0)
           gas_cap(2,i) = sngl(2d0*gas_capgrey(i))
        else
           stop 'analytic_opacity: in_suol invalid'
        endif!}}}
     enddo
  elseif(in_opacanaltype=='line') then
     ! Highly structured line test: group opacities alternate in magnitude!{{{
     ! sigmaP = A*T^B*rho^C
     ! sigmaP_g = sigmaP*func_P(g), sigmaR_g = sigmaP
     gas_capgrey = in_gas_capcoef*gas_temp**in_gas_captpwr * &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        !
        !set odd group magnitudes (low)
        do ig = 1, grp_ng, 2
           gas_cap(ig,i) = sngl(gas_capgrey(i)*in_ldisp1)
        enddo
        !set even group magnitudes (high)
        do ig = 2, grp_ng, 2
           gas_cap(ig,i) = sngl(gas_capgrey(i)*in_ldisp2)
        enddo
        !
        !calculate Planck, Rosseland opacities
        gas_capgrey(i) = 0d0
        do ig = 1, grp_ng
           x1 = pc_h*pc_c*grp_wlinv(ig+1)/(pc_kb*gas_temp(i))
           x2 = pc_h*pc_c*grp_wlinv(ig)/(pc_kb*gas_temp(i))
           gas_capgrey(i) = gas_capgrey(i)+15d0*gas_cap(ig,i)* &
                specint(x1,x2,3,10)/pc_pi**4
        enddo
        !
     enddo!}}}
  else
    stop 'analytic_opacity: in_opacanaltype invalid'
  endif

end subroutine analytic_opacity
