      module sourcemod
c     ----------------
      implicit none
      integer :: src_ns,src_ninit
      integer :: src_nnew
      integer :: src_nnonth
      integer :: src_nsurf
      integer :: src_ninitnew=0
      integer*8 :: src_nvacantmin=huge(src_ns)
      integer :: src_nflux=0
      integer :: src_nsurftot
c
      integer,allocatable :: src_ivacant(:) !array of vacant particle array locations
      integer*8,allocatable :: src_nvacantall(:) !(nmpi) number of vacancies on each of the ranks
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
