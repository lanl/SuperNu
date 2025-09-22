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
  use gasmod, only:gas_vol,gas_ye0,gas_mass,gas_rho,gas_matsrc,gas_ncell,gas_decaygamma,gas_dynfr
  use rprocmod, only:heating_rate,v_grid
  implicit none
  integer,intent(in) :: it, nt
  real*8,intent(in) :: tcenter
!--------------------------------------------------
!- Calculate heating rate from analytic fit formula (interpolated).
!- See Rosswog & Korobkin (2023) for details
!--------------------------------------------------
  integer :: icol, i
  real*8 :: helparr(gas_ncell),vexp(gas_ncell),gas_matsrc_raw(gas_ncell), gas_matsrc_dyn(gas_ncell), gas_matsrc_wnd(gas_ncell)
  real*8 :: therm_frac(gas_ncell)
  real*8, parameter :: M_ref = 0.05d0*pc_msun  !< reference mass when computing vexp

!---------------------------------------------------
!- Radioactive species partitioning (hardcoded for now!)
!  See Wollaeger et al. (2018), Sec. 2.2
!---------------------------------------------------
  real*8, parameter :: A_alpha = 1.2d-11 ! [g cm^-3 s]
  real*8, parameter :: A_beta  = 1.3d-11 ! [g cm^-3 s]
  real*8, parameter :: A_ff    = 0.2d-11 ! [g cm^-3 s]
  real*8, parameter :: X_alpha = 0.05d0
  real*8, parameter :: X_beta  = 0.20d0
  real*8, parameter :: X_gamma = 0.35d0
  real*8, parameter :: X_ff    = 0.00d0
  real*8, parameter :: ye_w    = 0.37d0
  real*8, parameter :: ye_d    = 0.05d0 ! electron fractions of wind and dyn ejecta (hardcoded for now)

!-- reset material source
  vexp = 0.1d0 ! default value of the expansion velocity
  if (tcenter > 0) &
     where(gas_mass > 0) &
        vexp = (3*M_ref/(4*pc_pi) * gas_vol/gas_mass)**(1d0/3d0)/tcenter/pc_c

  ! truncate expansion velocity to the limits of the grid
  vexp = max(v_grid(1), vexp)
  vexp = min(v_grid(size(v_grid)), vexp)
  
!   write(*,*)                                                                                                                                                                                                                                                                                                                                                                                                
!   write(*,'(A,ES12.3)') 'vexp = ', vexp
!  write(*,'(A,ES12.3)') 'Ye = ', gas_ye0  
!  print *, "time = ", tcenter                                                                                                                                                                                                                                                                                                                                                                              
!  print *, "gas_ye0 = ", gas_ye0                                                                                                                                                                                                                                                                                                                                                                           
 ! print *, "vexp = ", vexp                                                                                                                                                                                                                                                                                                                                                                                 
  !print *, "hrate_rproc = ", heating_rate(vexp, gas_ye0, tcenter)
  
  ! kill external gamma source
  gas_decaygamma = 0d0
  helparr = 0d0
  therm_frac = 0d0
  gas_matsrc = 0d0
   
  !gas_matsrc_raw = gas_rho*heating_rate(vexp, gas_ye0, tcenter) !heating rate of mass_frac weighted Ye's of each component

  
  gas_matsrc_raw = gas_dynfr*heating_rate(vexp, ye_d, tcenter) + (1-gas_dynfr)*heating_rate(vexp, ye_w, tcenter) !mass frac weighted heating rates of each component summed
 ! gas_matsrc_wnd = (1-gas_dynfr)*heating_rate(vexp, ye_w, tcenter)
 ! gas_matsrc_dyn = gas_dynfr*heating_rate(vexp, ye_d, tcenter)
  !gas_matsrc_dyn = gas_dynfr*heating_rate(vexp, ye_d, tcenter)
  where(gas_mass<=0d0) gas_matsrc_raw = 0d0
  where(gas_rho>0d0) helparr = 2d0/(tcenter*gas_rho)
  gas_matsrc = gas_rho*gas_matsrc_raw
  !-- create external gamma source for grey mc transport
  gas_decaygamma = X_gamma*gas_matsrc

  !-- create thermalized material source
  !helparr = 2d0/(tcenter*gas_rho)
  where(helparr>0d0)
     therm_frac = X_alpha*dlog(1d0 + helparr*A_alpha)/(helparr*A_alpha) &
                 + X_beta *dlog(1d0 + helparr*A_beta) /(helparr*A_beta) &
                 + X_ff   *dlog(1d0 + helparr*A_ff)   /(helparr*A_ff)
  endwhere
  gas_matsrc = gas_matsrc*therm_frac
     
     
  
  !print *, "thermfrac_rproc = ", therm_frac
  !print '(G10.3)', "rproc_gas_matsrc = ", gas_matsrc_raw
 ! write(*,*)
 ! write(*,'(A,ES12.3)') 'rproc_gasmat_src = ', gas_matsrc
 ! write(*,*)
 ! write(*,'(A,ES12.3)') 'rproc_gasmat_src_raw = ', gas_matsrc_raw
 ! write(*,*)
 ! write(*,'(A,ES12.3)') 'rproc_gasmat_src_dyn = ', gas_matsrc_dyn
 ! write(*,*)
 ! write(*,'(A,ES12.3)') 'rproc_gasmat_src_wnd = ', gas_matsrc_wnd
  

!  stop
  ! compute themralization fractions

end subroutine rprocess_fitfrm_source
