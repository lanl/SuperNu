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
      integer,private,parameter :: mreg = 11
      real*8,private,target :: registers(3,mreg)
c
c-- global-flow time registers:
      real*8,pointer :: t_gasupd(:) !update gas grid
      real*8,pointer :: t_eos(:) !equation of state
      real*8,pointer :: t_opac(:) !bound-bound opacity
      real*8,pointer :: t_bb(:) !bound-bound opacity
      real*8,pointer :: t_bf(:) !bound-free opacity
      real*8,pointer :: t_ff(:) !bound-free opacity
c-- packet transport
      real*8,pointer :: t_pckt_allrank(:) !collect the max runtimes across all ranks
      real*8,pointer :: t_pckt(:)
      real*8,pointer :: t_pcktnpckt(:)
      real*8,pointer :: t_pcktnddmc(:)
      real*8,pointer :: t_pcktnimc(:)
c
c-- parallel statistics packet timer
      real*8 :: t_pckt_stat(3)  !min,mean,max
c
      save
c
      contains
c
      subroutine timing_init
c     ----------------------
      implicit none
      t_gasupd => registers(:,1) !update gas grid
      t_eos =>    registers(:,2)   !equation of state
      t_opac =>   registers(:,3)    !bound-bound opacity
      t_bb =>     registers(:,4)    !bound-bound opacity
      t_bf =>     registers(:,5)    !bound-free opacity
      t_ff =>     registers(:,6)    !bound-free opacity
c--
      t_pckt_allrank => registers(:, 7)  !collect the max runtimes across all ranks
      t_pckt =>         registers(:, 8)
      t_pcktnpckt =>    registers(:, 9)
      t_pcktnddmc =>    registers(:,10)
      t_pcktnimc =>     registers(:,11)
      end subroutine timing_init
c
c
      subroutine timereg(reg,t)
c     -------------------------!{{{
      implicit none
      real*8,intent(inout) :: reg(3)
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
      subroutine timing_timestep(impi)
c     -----------------------------------
      implicit none
      integer,intent(in) :: impi
************************************************************************
* reset timestep timers and dump timing output (on master rank only).
************************************************************************
      logical :: lexist
      integer :: istat
c
c-- add to total
      registers(3,:) = registers(3,:) + registers(2,:)
c
c-- write output on master rank only
      if(impi==0) then
       inquire(file='output.timing',exist=lexist)
       open(4,file='output.timing',position='append',iostat=istat)
       if(istat/=0) stop 'timing_timestep: file open error'
c-- header
       if(.not.lexist) then
         write(4,'("#",30a12)') 't_gasupd','t_eos',
     &   't_opac','t_bb','t_bf','t_ff',
     &   't_p_allrank','t_pckt',
     &   't_pcktnpckt','t_pcktnddmc','t_pcktnimc'
       endif
c-- body
       write(4,'(1x,30g12.2)') registers(2,:)
       close(4)
      endif
c
c-- reset timers
      registers(2,:) = 0d0
      end subroutine timing_timestep
c
c
      subroutine print_timing
c     ------------------------
      implicit none
************************************************************************
* Print the timing totals
************************************************************************
      integer,parameter :: i=3 !total runtime timing
      write(6,*)
      write(6,*) 'timing results:'
      write(6,*) '============================'
      write(6,1) 'EOS               :',t_eos(i)
      write(6,1) 'opacity (bb|bf|ff):',t_opac(i),t_bb(i),t_bf(i),t_ff(i)
      write(6,*) '----------------------------'
      write(6,1) 'setup             :',t_setup
      write(6,1) 'gas-grid update   :',t_gasupd(i)
      write(6,1) 'packet transport  :',t_pckt_allrank(i)
      write(6,*) '----------------------------'
      write(6,1) 'all               :',t_all
1     format(1x,a,4f9.1)
      end subroutine print_timing
c
c
c
      subroutine time(t)
c     ------------------!{{{
      implicit none
      real*8,intent(out) :: t
************************************************************************
* Determine system time in seconds since some point in history.
* Note that the result is a single precision real, corresponding to the
* cpu_time() FORTRAN 95 intrinsic.
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
      t = icount*tick
      if(icount<icount_prev) t = t + tick*imax
      icount_prev = icount
c!}}}
      end subroutine time
c
      end module timingmod
