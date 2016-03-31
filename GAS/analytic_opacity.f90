!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2015 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
  real*8 :: x1(grp_ng),x2(grp_ng)  !unitless energy group bounds
  real*8 :: capcoef(gas_ncell)

  capcoef =0d0
  gas_cap = 0.
  gas_sig = 0d0

  !Calculating grey scattering opacity
  gas_sig = in_gas_sigcoef*gas_temp**in_gas_sigtpwr* &
       gas_rho**in_gas_sigrpwr

  !Calculating grouped Planck and Rosseland opacities
  if(in_opacanaltype=='none') then
     return
  elseif(in_opacanaltype=='grey') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        gas_cap(:,i) = sngl(capcoef(i))
     enddo

  elseif(in_opacanaltype=='mono') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
!-- monotonic group dependence
     do i = 1, gas_ncell
        x1 = pc_h*pc_c*grp_wlinv(2:)/pc_kb
        x2 = pc_h*pc_c*grp_wlinv(:grp_ng)/pc_kb
        gas_cap(:,i) = sngl(0.5d0*capcoef(i)*(x1+x2)/(x1*x2)**2)
     enddo!}}}
  elseif(in_opacanaltype=='pick') then
!-- sanity check
     if(grp_ng/=2) stop 'analytic_opacity: invalid grp_ng'
     if(grd_ny>1) stop 'analytic_opacity: no 2D for opacanaltyp=pick'
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr* &
          gas_rho**in_gas_caprpwr
!-- Su&Olson picket-fence distributions (tests: A,B,C (Su and Olson 1999))
     if(in_suol=='tsta') then    !Case: A
        gas_cap(1,:) = sngl(capcoef)
        gas_cap(2,:) = sngl(capcoef)
     elseif(in_suol=='tstb') then  !Case: B
        gas_cap(1,:) = sngl(2d0*capcoef/11d0)
        gas_cap(2,:) = sngl(20d0*capcoef/11d0)
     elseif(in_suol=='tstc') then  !Case: C
        gas_cap(1,:) = sngl(2d0*capcoef/101d0)
        gas_cap(2,:) = sngl(200d0*capcoef/101d0)
!-- added evacuated picket test
     elseif(in_suol=='tstd') then !Case: D (not in SuOlson)
        gas_cap(1,:) = sngl(0d0)
        gas_cap(2,:) = sngl(2d0*capcoef)
     else
        stop 'analytic_opacity: in_suol invalid'
     endif!}}}
  elseif(in_opacanaltype=='line') then
!-- using power law to make opacity coefficient
     capcoef = in_gas_capcoef*gas_temp**in_gas_captpwr * &
          gas_rho**in_gas_caprpwr
     do i = 1, gas_ncell
        !
        !set odd group magnitudes (low)
        do ig = 1, grp_ng, 2
           gas_cap(ig,i) = sngl(capcoef(i)*in_ldisp1)
        enddo
        !set even group magnitudes (high)
        do ig = 2, grp_ng, 2
           gas_cap(ig,i) = sngl(capcoef(i)*in_ldisp2)
        enddo
        !
     enddo!}}}
  else
    stop 'analytic_opacity: in_opacanaltype invalid'
  endif

end subroutine analytic_opacity
! vim: fdm=marker
