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
subroutine analytic_initial

  use gridmod
  use gasmod
  use timestepmod
  use inputparmod
  use physconstmod
  use manufacmod
  implicit none
!###############################################
! This subroutines attributes radiation energy to
! each cell and group depeding on user specification
! of in_srctype
!###############################################
!
!-- map radiation temperature to grd_evolinit
  call grid_volume(grd_igeom,grd_isvelocity,tsp_t)
  grd_evolinit = grd_evolinit * grd_vol
!--
!
!-- source specific initial conditions (overrides gas_inittyp)
!-- currently only supplying nonzero for in_srctype=manu
  if(any(['none','tabl','surf']==in_srctype)) then
     if(in_opacanaltype=='pick') then
!-- tstd initial energy profile currently approximation
        stop 'analytic_initial: in_opacanaltype==pick not implemented'
     else
        return
     endif
  elseif(in_srctype=='heav') then
     return
  elseif(in_srctype=='strt') then
     return
  elseif(in_srctype=='manu') then
     call init_manuprofile
  else
     stop 'analytic_initial: invalid in_srctype'
  endif

end subroutine analytic_initial
! vim: fdm=marker
