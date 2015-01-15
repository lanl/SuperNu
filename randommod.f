      module randommod
      implicit none
************************************************************************
* Random number generator based on mzran algorithm of Marsaglia and
* Zaman (1993)
************************************************************************
      integer*4 :: rnd_imax = 2147483579
c
c-- a data type for storing the state of the random number generator
      type :: rnd_t
       integer*4 :: part(4)
      end type rnd_t
c
      type(rnd_t),save :: rnd_state
c
c
      contains
c
c
      integer*4 function rnd_i(state)
c     -------------------------------!{{{
      implicit none
      type(rnd_t),intent(inout) :: state
************************************************************************
* Draws a uniform real number on [0,rnd_imax].
************************************************************************
      integer*4 :: imz
c
      imz = state%part(1) - state%part(3)
      if(imz<0) imz = imz + 2147483579
c
      state%part(1) = state%part(2)
      state%part(2) = state%part(3)
      state%part(3) = imz
      state%part(4) = 69069*state%part(4) + 1013904243
      imz = imz + state%part(4)
      rnd_i = imz!}}}
      end function rnd_i
c
c
      real*8 function rnd_r(state)
c     ----------------------------!{{{
      implicit none
      type(rnd_t),intent(inout) :: state
************************************************************************
* Draws a uniform real number on [0,1].
************************************************************************
      integer*4 :: imz
c
      imz = state%part(1) - state%part(3)
      if(imz<0) imz = imz + 2147483579
c
      state%part(1) = state%part(2)
      state%part(2) = state%part(3)
      state%part(3) = imz
      state%part(4) = 69069*state%part(4) + 1013904243
      imz = imz + state%part(4)
      rnd_r = 0.5d0 + 0.23283064d-9*imz !(0,1)!}}}
      end function rnd_r
c
c
      subroutine rnd_seed(state,i)
c     ----------------------------!{{{
      implicit none
      integer,intent(in) :: i
      type(rnd_t),intent(out) :: state
************************************************************************
* Return a random initial state of the random number generator.
************************************************************************
      integer :: j
      type(rnd_t) :: st
c-- init
      state%part = [521288629, 362436069, 16163801, 1131199299]
c-- advance to the selected offset
      call rnd_advance(state,i*4)!}}}
c-- draw four random numbers
      do j=1,4
       st%part(j) = rnd_i(state)
      enddo
c-- save
      state = st
      end subroutine rnd_seed
c
c
      subroutine rnd_advance(state,n)
c     -------------------------------!{{{
      implicit none
      type(rnd_t),intent(inout) :: state
      integer,intent(in) :: n
************************************************************************
* advance the random number generator n steps
************************************************************************
      integer :: i
      integer*4 :: imz
c
      do i=1,n
       imz = state%part(1) - state%part(3)
       if(imz<0) imz = imz + 2147483579
c
       state%part(1) = state%part(2)
       state%part(2) = state%part(3)
       state%part(3) = imz
       state%part(4) = 69069*state%part(4) + 1013904243
       imz = imz + state%part(4)
      enddo!}}}
      end subroutine rnd_advance
c
c
      end module randommod
