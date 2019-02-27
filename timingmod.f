*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module timingmod
c     ----------------
      implicit none
************************************************************************
* Collection of runtime timers for different parts of the code.
************************************************************************
c-- system_clock helper variables
      integer,private :: icount_prev=-1,imax
      real*8,private :: tick
c
c-- one-time events
      real*8 :: t_setup
      real*8 :: t_all
c
c-- timeline
      integer,parameter,private :: ntimeline=7
      real*8 :: t_timelin(ntimeline+1)
      real*8 :: t_timeline(ntimeline)
c
      integer,private,parameter :: mreg = 17
      real*8,private,target :: registers(4,mreg)
c
c-- global-flow time registers:
      real*8,pointer :: t_gasupd(:)   !update gas
      real*8,pointer :: t_eos(:)      !equation of state
      real*8,pointer :: t_emitp(:)    !emission probability
      real*8,pointer :: t_opac(:)     !all opacity
      real*8,pointer :: t_opacleak(:) !leakage opacity
      real*8,pointer :: t_bb(:)       !bound-bound opacity
      real*8,pointer :: t_bf(:)       !bound-free opacity
      real*8,pointer :: t_ff(:)       !free-free opacity
c-- communication
      real*8,pointer :: t_mpibcast(:)
      real*8,pointer :: t_mpimisc(:)
      real*8,pointer :: t_mpireduc(:)
c-- packet transport
      real*8,pointer :: t_pcktmin(:) !collect the max runtimes across all ranks
      real*8,pointer :: t_pcktmea(:)
      real*8,pointer :: t_pcktmax(:)
      real*8,pointer :: t_pcktgam(:)  !gamma transport
c-- flux
      real*8,pointer :: t_fluxtally(:)
c-- output
      real*8,pointer :: t_output(:)
c
c-- parallel statistics packet timer
      real*8 :: t_pckt_stat(3)  !min,mean,max
c
      save
c
      contains
c
      subroutine timingmod_init
c     ----------------------
      implicit none
      t_gasupd =>   registers(:,1)
      t_eos =>      registers(:,2)
      t_emitp =>    registers(:,3)
      t_opacleak=>  registers(:,4)
      t_opac =>     registers(:,5)
      t_bb =>       registers(:,6)
      t_bf =>       registers(:,7)
      t_ff =>       registers(:,8)
      t_mpibcast => registers(:,9)
      t_mpimisc  => registers(:,10)
      t_mpireduc => registers(:,11)
      t_pcktgam =>  registers(:,12)
      t_pcktmin =>  registers(:,13)  !collect the max runtimes across all ranks
      t_pcktmea =>  registers(:,14)  !collect the mean runtimes across all ranks
      t_pcktmax =>  registers(:,15)  !collect the min runtimes across all ranks
      t_fluxtally =>registers(:,16)
      t_output =>   registers(:,17)
      end subroutine timingmod_init
c
c
      subroutine timereg(reg,t)
c     -------------------------!{{{
      implicit none
      real*8,intent(inout) :: reg(4)
      real*8,intent(in) :: t
************************************************************************
* Put the time t in the register reg. The first position in reg stores
* the last value of t, the second position stores the sum.
* t values.
*
* Note that t is single precision (in correspondence with cpu_time)!
************************************************************************
      reg(1) = dble(t)
      reg(2) = reg(2) + dble(t)!}}}
      end subroutine timereg
c
c
      subroutine timing_cycle(impi,ldummystep)
c     -------------------------------------------
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
      registers(3,:) = registers(3,:) + registers(2,:)
c-- update max
      registers(4,:) = max(registers(4,:),registers(2,:))
c
c-- write output on master rank only
      if(impi==0) then
       inquire(file='output.timing',exist=lexist)
       open(4,file='output.timing',position='append',iostat=istat)
       if(istat/=0) stop 'timing_timestep: file open error'
c-- header
       if(.not.lexist) then
         write(4,'("#",30a12)') 't_gasupd','t_eos','t_emitp',
     &   't_opacleak','t_opac','t_bb','t_bf','t_ff',
     &   't_mpibcast','t_mpimisc','t_mpireduc',
     &   't_pgam','t_pmin','t_pmean','t_pmax','t_fluxtally','t_output'
       endif
c-- body
       if(ldummystep) then
        write(4,'("#",1p,30g12.4)') (registers(2,i),i=1,mreg)
       else
        write(4,'(1x,1p,30g12.4)') (registers(2,i),i=1,mreg)
       endif
       close(4)
      endif
c
c-- reset timers
      registers(:2,:) = 0d0
      end subroutine timing_cycle
c
c
      subroutine print_timing
c     ------------------------
      implicit none
************************************************************************
* Print the timing totals
************************************************************************
      integer,parameter :: i=3 !total runtime timing
      real*8 :: tmpi,taccounted
c
      tmpi = t_mpibcast(i)+t_mpimisc(i)+t_mpireduc(i)
      taccounted = tmpi+t_setup+t_gasupd(i)+t_opacleak(i)+t_pcktmax(i)+
     &  t_pcktgam(i)+t_output(i)
c
      write(6,*)
      write(6,*) 'timing results:'
      write(6,*) '============================'
      write(6,1) 'EOS               :',t_eos(i)
      write(6,1) 'opacity (bb|bf|ff):',t_opac(i),t_bb(i),t_bf(i),t_ff(i)
      write(6,*) '----------------------------'
      write(6,1) 'timeline          :',t_timeline
      write(6,*) '----------------------------'
      write(6,1) 'setup             :',t_setup
      write(6,1) 'gas update        :',t_gasupd(i)
      write(6,1) 'gas opacleak      :',t_opacleak(i)
      write(6,1) 'mpi (bc|misc|red) :',tmpi,
     &  t_mpibcast(i),t_mpimisc(i),t_mpireduc(i)
      write(6,1) 'transport min|max :',t_pcktmin(i),t_pcktmea(i),
     &  t_pcktmax(i)
      write(6,1) 'gamma transport   :',t_pcktgam(i)
      write(6,1) 'flux tally        :',t_fluxtally(i)
      write(6,1) 'output            :',t_output(i)
      write(6,1) 'unaccounted       :',t_all - taccounted
      write(6,*) '----------------------------'
      write(6,1) 'all               :',t_all
1     format(1x,a,10f9.1)
      end subroutine print_timing
c
c
c
      function t_time()
c     ------------------!{{{
      implicit none
      real*8 :: t_time
************************************************************************
* Determine system time in seconds since some point in history.
* Note that the result is a single precision real, corresponding to the
* cpu_t_time() FORTRAN 95 intrinsic.
************************************************************************
      integer :: icount,irate
c
      if(icount_prev<0) then
       call system_clock(count=icount, count_rate=irate, count_max=imax)
       tick = 1.d0/irate
      else
       call system_clock(icount)
      endif
c
      t_time = icount*tick
      if(icount<icount_prev) t_time = t_time + tick*imax
      icount_prev = icount
c!}}}
      end function t_time
c
      end module timingmod
c vim: fdm=marker
