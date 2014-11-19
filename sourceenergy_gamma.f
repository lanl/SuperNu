      subroutine sourceenergy_gamma
c     -----------------------------
      use gridmod
      implicit none
cc
cc-- dump integral numbers
c      write(6,*) sum(grd_edep),sum(grd_emitex),
c     &  sum(grd_edep)/sum(grd_emitex)
c
c-- gamma deposition is energy source
      grd_emit = grd_emit + grd_edep
c
      end subroutine sourceenergy_gamma
