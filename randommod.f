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
c
      contains
c
c
      integer*4 function irand(state)
c     --------------------------------------------!{{{
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
      imz = imz + state%part(4)!}}}
      end function irand
c
c
      real*8 function rrand(state)
c     ------------------------------------!{{{
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
      rrand = 0.5d0 + 0.23283064d-9*imz !(0,1)!}}}
      end function rrand
c
      end module randommod
