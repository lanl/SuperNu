! © 2024. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
! Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
! Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
! copies to the public, perform publicly and display publicly, and to permit others to do so.
!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.

subroutine rprocess_fitfrm_source(it,nt,tcenter)
  use physconstmod, only:pc_msun,pc_pi,pc_c
  use gasmod, only:gas_ye0,gas_mass,gas_rho,gas_matsrc,gas_ncell
  use rprocmod, only:heating_rate,v_grid
  use gridmod, only:grd_vol
  implicit none
  integer,intent(in) :: it, nt
  real*8,intent(in) :: tcenter
!--------------------------------------------------
!- Calculate heating rate from analytic fit formula (interpolated).
!- See Rosswog & Korobkin (2023) for details
!--------------------------------------------------
  integer :: icol, i
  real*8 :: helparr(gas_ncell),vexp(gas_ncell)
  real*8 :: therm_frac(gas_ncell)
  real*8, parameter :: M_ref = 0.05d0*pc_msun  !< reference mass when computing vexp

!---------------------------------------------------
!- Radioactive species partitioning (hardcoded for now!)
!  See Wollaeger et al. (2018), Sec. 2.2 
!---------------------------------------------------
  real*8, parameter :: A_alpha = 1.3d-11 ! [g cm^-3 s]
  real*8, parameter :: A_beta  = 1.3d-11 ! [g cm^-3 s]
  real*8, parameter :: A_ff    = 0.2d-11 ! [g cm^-3 s]
  real*8, parameter :: X_alpha = 0.05d0  ! [g cm^-3 s]
  real*8, parameter :: X_beta  = 0.20d0  ! [g cm^-3 s]
  real*8, parameter :: X_ff    = 0.00d0  ! [g cm^-3 s]

!-- reset material source
  vexp = 0.1d0 ! default value of the expansion velocity
  if (tcenter > 0) &
     where(gas_mass > 0) &
        vexp = (3*M_ref/(4*pc_pi) * grd_vol/gas_mass)**(1d0/3d0)/tcenter/pc_c

  ! truncate expansion velocity to the limits of the grid
  vexp = max(v_grid(1), vexp)
  vexp = min(v_grid(size(v_grid)), vexp)

  where(gas_rho > 0d0) 
     gas_matsrc = gas_rho*heating_rate(vexp, gas_ye0, tcenter)
     helparr = 2d0/(tcenter*gas_rho)
     therm_frac = X_alpha*dlog(1d0 + helparr*A_alpha)/(helparr*A_alpha) &
                + X_beta *dlog(1d0 + helparr*A_beta) /(helparr*A_beta) &
                + X_ff   *dlog(1d0 + helparr*A_ff)   /(helparr*A_ff)
     gas_matsrc = gas_matsrc*therm_frac
  endwhere

  ! compute themralization fractions

end subroutine rprocess_fitfrm_source
