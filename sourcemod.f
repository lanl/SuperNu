      module sourcemod
c     ----------------
      implicit none
      integer :: src_ns,src_ninit
      integer :: src_nsurf,src_nnonth,src_nnew,src_ninitnew
      integer :: src_nvacantmin=huge(src_ns)
      integer :: src_nflux=0
c
      integer,allocatable :: src_ivacant(:) !array of vacant particle array locations
      integer,allocatable :: src_nvacantall(:) !(nmpi) number of vacancies on each of the ranks
c
      save
c
      contains

      subroutine sourcemod_init(nmpi)
c     ----------------------------------------
      implicit none
      integer,intent(in) :: nmpi
************************************************************************
* init particle module
************************************************************************
      allocate(src_nvacantall(nmpi))
      end subroutine sourcemod_init
c
      end module sourcemod
