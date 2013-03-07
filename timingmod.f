      module timingmod
c     ----------------
      IMPLICIT NONE
************************************************************************
* Collection of runtime timers for different parts of the code.
************************************************************************
c-- system_clock helper variables
      integer,private :: icount_prev=-1,imax
      REAL*8,private :: tick
c
c-- one-time events
      REAL*8 :: t_setup
      REAL*8 :: t_all
c-- global-flow time registers:
      REAL*8 :: t_gasupd(2)=0d0  !update gas grid
      REAL*8 :: t_pckt(2)=0d0   !packet transport
c-- task-based time registers:
      REAL*8 :: t_eos(2)=0d0    !equation of state
      REAL*8 :: t_bb(2)=0d0     !bound-bound opacity
      REAL*8 :: t_bf(2)=0d0     !bound-free opacity
      REAL*8 :: t_ff(2)=0d0     !bound-free opacity
c
c-- parallel statistics packet timer
      REAL*8 :: t_pckt_stat(3)!min,mean,max
c
      save
c
      contains
c
      subroutine timereg(reg,t)
c     -------------------------
      IMPLICIT NONE
      REAL*8,intent(inout) :: reg(2)
      real,intent(in) :: t
************************************************************************
* Put the time t in the register reg. The first position in reg stores
* the last value of t, the second position stores the sum.
* t values.
*
* Note that t is single precision (in correspondence with cpu_time)!
************************************************************************
      reg(1) = dble(t)
      reg(2) = reg(2) + dble(t)
      end subroutine timereg
c
c
      subroutine print_timing
c     -----------------------
      IMPLICIT NONE
************************************************************************
* Print the timing totals
************************************************************************
      write(6,*)
      write(6,*) 'timing results:'
      write(6,*) '============================'
      write(6,1) 'EOS               :',t_eos(2)
      write(6,1) 'opacity (bb|bf|ff):',t_bb(2)+t_bf(2)+t_ff(2),
     &                                 t_bb(2),t_bf(2),t_ff(2)
      write(6,*) '----------------------------'
      write(6,1) 'setup             :',t_setup
      write(6,1) 'gas-grid update   :',t_gasupd(2)
      write(6,1) 'packet transport  :',t_pckt(2)
      write(6,*) '----------------------------'
      write(6,1) 'all               :',t_all
1     format(1x,a,4f9.1)
      end subroutine print_timing
c
c
c
      subroutine time(t)
c     ------------------
      IMPLICIT NONE
      real,intent(out) :: t
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
c
      end subroutine time
c
      end module timingmod
