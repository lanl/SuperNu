!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine vacancies

  use particlemod
  use sourcemod
  implicit none

!##################################################
  !This subroutine creates an array of vacant index locations
  !in the particle array each time step.
!##################################################

  integer :: ipart, ivac

!-- Initializing index counters and full array checking boolean
  ipart = 0
  ivac = 0

!-- Filling src_ivacant with particle index of vacant particles: loop
  if(src_nnew > prt_npartmax) stop 'vacancies: nnew>npartmax'
  do ipart=1,prt_npartmax
     if(.not.prt_isvacant(ipart)) cycle
     ivac = ivac+1
     src_ivacant(ivac) = ipart
     if (ivac == src_nnew) exit
  enddo
  if(ipart > prt_npartmax) stop 'Maximum number of prt_particles reached'

end subroutine vacancies
! vim: fdm=marker
