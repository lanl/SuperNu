*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module countersmod
c     ------------------
      implicit none
************************************************************************
* Collection of runtime counters for different parts of the code.
* This is an integer clone of timingmod.
************************************************************************
      integer,private,parameter :: mreg = 11
      integer*8,private,target :: iregisters(4,mreg)
c
      integer*8,pointer :: ct_nnonvacant(:) !non-vacant particle slots
      integer*8,pointer :: ct_npcreate(:)   !newly created particles
      integer*8,pointer :: ct_npactive(:)   !transported particles
      integer*8,pointer :: ct_npflux(:)     !flux particles
      integer*8,pointer :: ct_npfluxbuf(:)  !flux particles remaining in buffer
      integer*8,pointer :: ct_npcensimc(:)  !censussed IMC particles
      integer*8,pointer :: ct_npcensddmc(:) !censussed DDMC particles
      integer*8,pointer :: ct_npstepimc(:)  !IMC interactions (steps)
      integer*8,pointer :: ct_npstepddmc(:) !DDMC interactions (steps)
      integer*8,pointer :: ct_npstepmax(:)  !interactions (steps)
c-- transport method swaps
      integer*8,pointer :: ct_npmethswap(:)
c
      save
c
      contains
c
      subroutine countersmod_init
c     ---------------------------
      implicit none
      ct_nnonvacant  => iregisters(:,1)
      ct_npcreate    => iregisters(:,2)
      ct_npactive    => iregisters(:,3)
      ct_npflux      => iregisters(:,4)
      ct_npfluxbuf   => iregisters(:,5)
      ct_npcensimc   => iregisters(:,6)
      ct_npcensddmc  => iregisters(:,7)
      ct_npstepimc   => iregisters(:,8)
      ct_npstepddmc  => iregisters(:,9)
      ct_npstepmax   => iregisters(:,10)
      ct_npmethswap  => iregisters(:,11)
      end subroutine countersmod_init
c
c
      subroutine counterreg(ireg,n)
c     -----------------------------!{{{
      implicit none
      integer*8,intent(inout) :: ireg(4)
      integer,intent(in) :: n
************************************************************************
* Put the counter c in the register reg. The first position in reg stores
* the last value of n, the second position stores the sum.
************************************************************************
      ireg(1) = n
      ireg(2) = ireg(2) + n!}}}
      end subroutine counterreg
c
c
      subroutine counters_cycle(impi,ldummystep)
c     --------------------------------------------
      implicit none
      integer,intent(in) :: impi
      logical,intent(in) :: ldummystep
************************************************************************
* reset timestep timers and dump timing output (on master rank only).
************************************************************************
      logical :: lexist
      integer :: istat,i
c
c-- add to total
      iregisters(3,:) = iregisters(3,:) + iregisters(2,:)
c-- update max
      iregisters(4,:) = max(iregisters(4,:),iregisters(2,:))
c
c-- write output on master rank only
      if(impi==0) then
       inquire(file='output.counters',exist=lexist)
       open(4,file='output.counters',position='append',iostat=istat)
       if(istat/=0) stop 'counters_timestep: file open error'
c-- header
       if(.not.lexist) then
         write(4,'("#",30a12)') 'nonvacant','create','transport',
     &     'flux','fluxbuffer','censusimc','censusddmc','stepimc',
     &     'stepddmc','stepmax','methswap'
       endif
c-- body
       if(ldummystep) then
        write(4,'("#",1p,30g12.4)') (iregisters(2,i),i=1,mreg)
       else
        write(4,'(1x,1p,30g12.4)') (iregisters(2,i),i=1,mreg)
       endif
       close(4)
      endif
c
c-- reset timers
      iregisters(:2,:) = 0
      end subroutine counters_cycle
c
c
      subroutine print_counters
c     -------------------------
      implicit none
************************************************************************
* Print the timing totals
************************************************************************
      integer,parameter :: i=3 !total runtime timing
c
      write(6,*)
      write(6,*) 'counters (on master rank):'
      write(6,*) '============================'
      write(6,1) 'pckt transport    :',ct_npactive(i)/1000
      write(6,1) ' steps (imc|ddmc) :',ct_npstepimc(i)/1000,
     &  ct_npstepddmc(i)/1000
      write(6,2) ' max steps        :',ct_npstepmax(4)
      write(6,1) ' method swaps     :',ct_npmethswap(i)/1000
1     format(1x,a,10(i10,"k"))
2     format(1x,a,10i10)
      end subroutine print_counters
c
      end module countersmod
c vim: fdm=marker
