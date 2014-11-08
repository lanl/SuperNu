subroutine analytic_opacity

  use gridmod
  use gasgridmod
  use physconstmod
  use miscmod, only:specint
  !use timestepmod
  implicit none

!#####################################
  !This subroutine computes Planck and Rosseland
  !opacities for several simple test group structures.
  !The power law coefficients are used in any structure
  !selection to either fully or partially determine
  !opacity dependence on temperature and density.
  !
  !Since revision 121, calculates grey scattering opacity, dd_sig
!#####################################

  integer :: i, ig
  real*8 :: x1, x2  !unitless energy group bounds

  dd_siggrey = 0d0
  dd_cap = 0d0
  dd_sig = 0d0

  !Calculating grey scattering opacity
  dd_sig = gas_sigcoefs*dd_temp**gas_sigtpwrs* &
       dd_rho**gas_sigrpwrs

  !Calculating grouped Planck and Rosseland opacities
  if(gas_opacanaltype=='none') then
     return
  elseif(gas_opacanaltype=='grey') then
     ! sigmaP = A*T^B*rho^C (A,B,C set in input.par)
     ! sigmaP_g, sigmaR_g = sigmaP for all g 
     ! Input wavelength grid not used
     dd_siggrey = gas_sigcoef*dd_temp**gas_sigtpwr* &
          dd_rho**gas_sigrpwr
     do i = 1, dd_ncell
        dd_cap(:,i) = dd_siggrey(i)
     enddo

  elseif(gas_opacanaltype=='mono') then
     ! sigmaP = A*T^B*rho^C*f(T)!{{{
     ! sigmaP_g = sigmaP*func_P(T,g), sigmaR_g=sigmaP*func_R(T,g)
     ! func_P(T,g) and func_R(T,g) are functions proportional to
     ! integral_g(1/nu^3)
     !
     dd_siggrey = gas_sigcoef*dd_temp**gas_sigtpwr* &
          dd_rho**gas_sigrpwr
     do i = 1, dd_ncell
        do ig = 1, gas_ng
            x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb)
            x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb)
            dd_cap(ig,i) = 0.5d0*dd_siggrey(i)*(x1+x2)/(x1*x2)**2
        enddo
        dd_siggrey(i) = 0d0
        do ig = 1, gas_ng
           x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*dd_temp(i))
           x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*dd_temp(i))
           dd_siggrey(i) = dd_siggrey(i)+15d0*dd_cap(ig,i)* &
                specint(x1,x2,3)/pc_pi**4
        enddo
     enddo!}}}
  elseif(gas_opacanaltype=='pick') then
     if(gas_ng/=2) stop 'analytic_opacity: invalid gas_ng'
     ! sigmaP = sigmaR = constant = A (dd_sigcoef set in input.par)!{{{
     ! Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     ! Input wavelength grid not used
     if(grd_ny>1) stop 'analytic_opacity: no 2D for opacanaltyp=pick'
     do i = 1, dd_ncell
        dd_siggrey(i) = gas_sigcoef*dd_temp(i)**gas_sigtpwr* &
             dd_rho(i)**gas_sigrpwr     
        if(gas_suol=='tsta') then    !Case: A
           dd_cap(1,i) = dd_siggrey(i)
           dd_cap(2,i) = dd_siggrey(i)
        elseif(gas_suol=='tstb') then  !Case: B
           dd_cap(1,i) = 2d0*dd_siggrey(i)/11d0
           dd_cap(2,i) = 20d0*dd_siggrey(i)/11d0
        elseif(gas_suol=='tstc') then  !Case: C
           dd_cap(1,i) = 2d0*dd_siggrey(i)/101d0
           dd_cap(2,i) = 200d0*dd_siggrey(i)/101d0
!-- added evacuated picket test
        elseif(gas_suol=='tstd') then !Case: D (not in SuOlson)
           dd_cap(1,i) = 0d0
           dd_cap(2,i) = 2d0*dd_siggrey(i)
        else
           stop 'analytic_opacity: gas_suol invalid'
        endif!}}}
     enddo
  elseif(gas_opacanaltype=='line') then
     ! Highly structured line test: group opacities alternate in magnitude!{{{
     ! sigmaP = A*T^B*rho^C
     ! sigmaP_g = sigmaP*func_P(g), sigmaR_g = sigmaP
     dd_siggrey = gas_sigcoef*dd_temp**gas_sigtpwr * &
          dd_rho**gas_sigrpwr
     do i = 1, dd_ncell
        !
        !set odd group magnitudes (low)
        do ig = 1, gas_ng, 2
           dd_cap(ig,i) = dd_siggrey(i)*gas_ldisp1
        enddo
        !set even group magnitudes (high)
        do ig = 2, gas_ng, 2
           dd_cap(ig,i) = dd_siggrey(i)*gas_ldisp2
        enddo
        !
        !calculate Planck, Rosseland opacities
        dd_siggrey(i) = 0d0
        do ig = 1, gas_ng
           x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*dd_temp(i))
           x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*dd_temp(i))
           dd_siggrey(i) = dd_siggrey(i)+15d0*dd_cap(ig,i)* &
                specint(x1,x2,3)/pc_pi**4
        enddo
        !
     enddo!}}}
  else
    stop 'analytic_opacity: gas_opacanaltype invalid'
  endif

end subroutine analytic_opacity
