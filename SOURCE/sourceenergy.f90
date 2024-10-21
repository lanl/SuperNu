! © 2023. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
! Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
! Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
! copies to the public, perform publicly and display publicly, and to permit others to do so.
!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine sourceenergy(lgamma)

  use gasmod
  use sourcemod
  use totalsmod
  use timestepmod
  use physconstmod
  use inputparmod
  use manufacmod
  implicit none
  logical, intent(in) :: lgamma
!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################
  real*8 :: q1,q2,q3

!-- prepare manufactured solution temperature source
  if(in_srctype=='manu') then
     call generate_manutempsrc(in_str_totmass,in_gas_capcoef,tsp_t,tsp_dt)
  elseif(in_srctype=='tabl') then
     call tabular_source(tsp_it,tsp_nt,tsp_tcenter,in_srctype)
  elseif(in_srctype=='rpro') then
     call rprocess_fitfrm_source(tsp_it,tsp_nt,tsp_tcenter)
  elseif(in_gas_srccoef>0d0) then
!-- short-cuts
     q1=in_gas_srcrpwr
     q2=in_gas_srctpwr
     q3=in_gas_srctimepwr
!-- power-law material energy source
     gas_matsrc=in_gas_srccoef*gas_rho**q1 * &
        gas_temp**q2 * (tsp_t+.5d0*tsp_dt)**q3
  endif

! Calculating fictitious emission energy per cell
!-- thermal source
  gas_emit =  tsp_dt*gas_vol*gas_fcoef*gas_capgrey*pc_c*gas_ur
  tot_sthermal = sum(gas_emit)
!
!-- manufactured solution material source
  gas_emit = gas_emit + tsp_dt*gas_vol*(1d0-gas_fcoef)*gas_matsrc
  tot_smanufac = sum(tsp_dt*gas_vol*(1d0-gas_fcoef)*gas_matsrc)
!
!-- accounting for total material source energy:
!-- gas_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-gas_fcoef is thermal radiation.
  tot_eext = tot_eext + tsp_dt*sum(gas_vol*gas_matsrc)

!-- non-thermal decay radiation source energy
  if(lgamma) then
     if(in_srctype=='tabl' .or. in_srctype=='rpro') then
!-- gamma source calculated by tabular_source in erg/cc/s
        gas_decaygamma = gas_decaygamma*tsp_dt*gas_vol
        gas_decaybeta = 0d0
     elseif(in_sgamcoef>0d0) then
!-- power-law
!-- short-cuts
        q1=in_sgamrpwr
        q2=in_sgamtpwr
        q3=in_sgamtimepwr
        gas_decaygamma = tsp_dt*gas_vol * &
             in_sgamcoef*gas_rho**q1 * &
             gas_temp**q2 * (tsp_t+.5d0*tsp_dt)**q3
        gas_decaybeta = 0d0
     endif
!-- gamma source is considered external energy
     gas_emitex = gas_decaygamma  !grey transport
     gas_emit = gas_emit + gas_decaybeta  !local deposition
!-- totals
     tot_sdecaygamma = sum(gas_decaygamma)
     tot_sdecaybeta = sum(gas_decaybeta)
  else
     gas_emitex = 0d0
     tot_sdecaygamma = 0d0
     tot_sdecaybeta = 0d0
  endif

end subroutine sourceenergy
! vim: fdm=marker
