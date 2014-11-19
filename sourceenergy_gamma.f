      subroutine sourceenergy_gamma
c     -----------------------------
      use gridmod
      use timestepmod
      use profiledatamod
      implicit none
************************************************************************
* Add the energy deposition from gamma absorption to the energy source
* for optical particles.
************************************************************************
!     integer :: i
!     real*8 :: hlparr(grd_nx)
!     real*8 :: help
!c
!c-- dump integral numbers
!      help = sum(grd_emitex)
!      hlparr = gamma_profile(tsp_t)
!      write(6,*) sum(grd_emitex),sum(grd_edep),sum(hlparr)*help,
!     &  sum(grd_edep)/sum(grd_emitex),
!     &  sum(hlparr)*help/sum(grd_emitex)
!      do i=grd_nx,1,-1
!       write(6,*) 65-i,grd_emitex(i,1,1)/tsp_dt,grd_edep(i,1,1)/tsp_dt,
!     &   hlparr(i)*help/tsp_dt,
!     &   grd_edep(i,1,1)/grd_emitex(i,1,1),
!     &   hlparr(i)*help/grd_emitex(i,1,1)
!      enddo
c
c-- gamma deposition is energy source
      grd_emit = grd_emit + grd_edep
c
      end subroutine sourceenergy_gamma
