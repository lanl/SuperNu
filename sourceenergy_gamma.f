      subroutine sourceenergy_gamma
c     -----------------------------
      use gridmod
!     use timestepmod
      implicit none
************************************************************************
* Add the energy deposition from gamma absorption to the energy source
* for optical particles.
************************************************************************
!     integer :: i
c
c-- dump whole profile (1D only)
!      do i=grd_nx,1,-1
!       write(6,*) 65-i,grd_emitex(i,1,1)/tsp_dt,grd_edep(i,1,1)/tsp_dt,
!     &   grd_edep(i,1,1)/grd_emitex(i,1,1)
!      enddo
c
c-- gamma deposition is energy source
      grd_emit = grd_emit + grd_edep
c
      end subroutine sourceenergy_gamma
