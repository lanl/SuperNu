      module countersmod
c     ------------------
      implicit none
************************************************************************
* Collection of runtime counters for different parts of the code.
* This is an integer clone of timingmod.
************************************************************************
      integer,private,parameter :: mreg = 9
      integer,private,target :: iregisters(4,mreg)
c
c-- transported particles
      integer,pointer :: ct_npcreate(:)
      integer,pointer :: ct_nptransport(:)
c-- flux particles
      integer,pointer :: ct_npflux(:)
c-- censussed particles
      integer,pointer :: ct_npcensimc(:)
      integer,pointer :: ct_npcensddmc(:)
c-- transport interactions (steps)
      integer,pointer :: ct_npstepimc(:)
      integer,pointer :: ct_npstepddmc(:)
      integer,pointer :: ct_npstepmax(:)
c-- transport method swaps
      integer,pointer :: ct_npmethswap(:)
c
      save
c
      contains
c
      subroutine countersmod_init
c     ---------------------------
      implicit none
      ct_npcreate    => iregisters(:,1)
      ct_nptransport => iregisters(:,2)
      ct_npflux      => iregisters(:,3)
      ct_npcensimc   => iregisters(:,4)
      ct_npcensddmc  => iregisters(:,5)
      ct_npstepimc   => iregisters(:,6)
      ct_npstepddmc  => iregisters(:,7)
      ct_npstepmax   => iregisters(:,8)
      ct_npmethswap  => iregisters(:,9)
      end subroutine countersmod_init
c
c
      subroutine counterreg(ireg,n)
c     -----------------------------!{{{
      implicit none
      integer,intent(inout) :: ireg(4)
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
      subroutine counters_timestep(impi)
c     ----------------------------------
      implicit none
      integer,intent(in) :: impi
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
         write(4,'("#",30a12)') 'npcreate','nptransport',
     &     'npflux','npcensimc','npcensddmc','npstepimc',
     &     'npstepddmc','nstepmax','npmethswap'
       endif
c-- body
       write(4,'(1x,1p,30g12.4)') (iregisters(2,i),i=1,mreg)
       close(4)
      endif
c
c-- reset timers
      iregisters(2,:) = 0
      end subroutine counters_timestep
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
      write(6,1) 'pckt transport    :',ct_nptransport(i)/1000
      write(6,1) ' steps (imc|ddmc) :',ct_npstepimc(i)/1000,
     &  ct_npstepddmc(i)/1000
      write(6,2) ' max steps        :',ct_npstepmax(4)
      write(6,1) ' method swaps     :',ct_npmethswap(i)/1000
1     format(1x,a,10(i10,"k"))
2     format(1x,a,10i10)
      end subroutine print_counters
c
      end module countersmod
