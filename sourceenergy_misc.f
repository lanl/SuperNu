      subroutine sourceenergy_misc
c     ----------------------------
      use mpimod
      use gridmod
      use totalsmod
!     use timestepmod
      implicit none
************************************************************************
* Add the energy deposition from gamma absorption and amplification
* factors to the energy source for optical particles.
************************************************************************
!     integer :: i,l
c
c-- dump whole profile (1D only)
!      do i=grd_nx,1,-1
!       l = grd_icell(i,1,1)
!       write(6,*) 65-i,grd_emitex(l)/tsp_dt,grd_edep(l)/tsp_dt,
!     &   grd_edep(l)/grd_emitex(l)
!      enddo
c
c-- sanity check energy deposition
      if(any(grd_edep<0d0)) stop 'sourceenergy_misc: negative energy'
c
c-- gamma deposition is energy source
      grd_emit = grd_emit + grd_edep
c-- clear eamp in the dummy cell
      if(grd_lpad) grd_eamp(grd_nc) = 0d0
c-- 'particle-amplification' factor
      grd_emit = grd_emit + grd_eamp
c-- verify zero emission energy in dummy cell
      if(grd_lpad .and. grd_emit(grd_nc)/=0d0)
     &  stop 'soureceenergy_misc: emission energy in dummy cell'
c
c-- add gamma radiation source tot total
      if(impi==impi0) tot_eext = tot_eext + sum(grd_edep)
c
      end subroutine sourceenergy_misc
