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

      subroutine sourcemod_init(nmpi,ns,ninit)
c     ----------------------------------------!{{{
      implicit none
      integer,intent(in) :: nmpi,ns,ninit
************************************************************************
* init particle module
************************************************************************
!-- adopt input values in module internal storage
      src_ns = ns
      src_ninit = ninit
      allocate(src_nvacantall(nmpi))
!}}}
      end subroutine sourcemod_init
c
      end module sourcemod
