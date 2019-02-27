!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
  if(in_srctype=='none'.or.in_srctype=='surf') then
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
