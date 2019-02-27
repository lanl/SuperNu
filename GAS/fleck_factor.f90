!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine fleck_factor(tempalt,capgreyalt)

  use inputparmod
  use gasmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------
  real*8, intent(in) :: tempalt(gas_ncell)
  real*8, intent(in) :: capgreyalt(gas_ncell)
  integer :: i
  real*8 :: Um, beta, beta2, dlogsig

!-- init (necessary for domain decomposition
  gas_fcoef = 0d0

!-- calculating modified Fleck factor
  do i=1,gas_ncell
     if(gas_mass(i)<=0d0) cycle
!
     Um = gas_bcoef(i)*gas_temp(i)
     if(gas_temp(i)<=0d0.or.gas_bcoef(i)==0d0) then
        beta = 0d0
     elseif(gas_capgrey(i)<=0d0.or.capgreyalt(i)<=0d0) then
        beta = 4.0*gas_ur(i)/Um
     else
        if(.not.in_ismodimc .or. gas_temp(i)==tempalt(i) .or. gas_eraddens(i)==0d0) then
           beta2 = 0d0
        else
           dlogsig = log(gas_capgrey(i)/(capgreyalt(i)*gas_rho(i)))/& !-- convert from per gram
              (gas_temp(i)-tempalt(i))
           beta2 = min(0d0,(gas_eraddens(i)-gas_ur(i))*&
              dlogsig/gas_bcoef(i))
        endif
        beta = 4.0*gas_ur(i)/Um - beta2

     endif
     gas_fcoef(i) = 1.0/(1.0+in_alpha*beta*pc_c*tsp_dt*gas_capgrey(i))
!-- sanity check
     if(gas_fcoef(i)<=0d0) stop 'fleck_factor: fcoef<=0d0'
     if(gas_fcoef(i)/=gas_fcoef(i)) stop 'fleck_factor: fcoef nan'
  enddo !i

end subroutine fleck_factor
! vim: fdm=marker
