      module sourcemod
c     ----------------
      implicit none
      integer :: src_ns,src_ninit
      integer :: src_nsurf,src_nnonth,src_nnew,src_ninitnew

      integer,allocatable :: src_ivacant(:) !array of vacant particle array locations
c
      save
c
      contains

      subroutine source_init(ns,ninit)
c     ------------------------------------------!{{{
      implicit none
      integer,intent(in) :: ns, ninit
************************************************************************
* init particle module
************************************************************************
!-- adopt input values in module internal storage
      src_ns = ns
      src_ninit = ninit
!}}}
      end subroutine source_init
c
      end module sourcemod
