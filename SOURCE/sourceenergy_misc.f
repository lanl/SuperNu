*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine sourceenergy_misc(lmpi0)
c     ----------------------------------
      use gridmod
      use totalsmod
!     use timestepmod
      implicit none
      logical,intent(in) :: lmpi0
************************************************************************
* Add the energy deposition from gamma absorption and amplification
* factors to the energy source for optical particles.
************************************************************************
!     integer :: i,l
c
c-- dump whole profile (1D only)
!      do i=grd_nx,1,-1
!       l = grd_icell(i,1,1)
!       write(6,*) 65-i,grd_emitex(l)/tsp_dt,grd_tally(1,l)/tsp_dt,
!     &   grd_tally(1,l)/grd_emitex(l)
!      enddo
c
c-- sanity check energy deposition
      if(any(grd_tally(1,:)<0d0)) stop 'sourceenergy_misc: energy<0'
c
c-- gamma deposition is energy source
      grd_emit = grd_emit + grd_tally(1,:)
      if(lmpi0) tot_sdeposgamma = sum(grd_tally(1,:))
c-- clear eamp in the dummy cell
      if(grd_ivoid>0) grd_eamp(grd_ncell) = 0d0
c-- 'particle-amplification' factor
      grd_emit = grd_emit + grd_eamp
      if(lmpi0) tot_samp = sum(grd_eamp)
c-- verify zero emission energy in dummy cell
      if(grd_ivoid>0 .and. grd_emit(grd_ncell)/=0d0)
     &  stop 'soureceenergy_misc: emission energy in dummy cell'
c
c-- add gamma radiation source tot total
      if(lmpi0) tot_eext = tot_eext + sum(grd_tally(1,:))
c
      end subroutine sourceenergy_misc
c vim: fdm=marker
