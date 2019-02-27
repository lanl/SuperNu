*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
      integer :: rnd_nstate
      type(rnd_t) :: rnd_state
      type(rnd_t),allocatable :: rnd_states(:)
c
      save
c
      contains
c
c
      subroutine rnd_init(n,ioffset)
c     ------------------------------!{{{
      implicit none
      integer,intent(in) :: n,ioffset
************************************************************************
* Initialize n states of the random number generator.
************************************************************************
      integer :: i,j
      type(rnd_t) :: state
c-- alloc
      rnd_nstate = n
      allocate(rnd_states(rnd_nstate))
c-- init
      state%part = [521288629, 362436069, 16163801, 1131199299]
c-- advance to the requested offset
      call rnd_advance(state,n*ioffset*4)
c-- draw four random numbers per state
      do i=1,n
       do j=1,4
        call rnd_i(rnd_states(i)%part(j),state)
       enddo
      enddo
!}}}
      end subroutine rnd_init
c
c
      pure subroutine rnd_i(i,state)
c     -------------------------------!{{{
      implicit none
      integer*4,intent(out) :: i
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
      i = imz!}}}
      end subroutine rnd_i
c
c
      pure subroutine rnd_r(x,state)
c     ----------------------------!{{{
      implicit none
      real*8,intent(out) :: x
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
      x = 0.5d0 + 0.23283064d-9*imz !(0,1)!}}}
      end subroutine rnd_r
c
c
      pure subroutine rnd_advance(state,n)
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
c vim: fdm=marker
