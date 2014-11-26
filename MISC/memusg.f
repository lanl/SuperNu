      function memusg() result(mbsize)
c     -------------------------
      use miscmod, only:warn
      implicit none
      integer :: mbsize(2)
************************************************************************
* Read memory statistics from /proc file system.
************************************************************************
      integer,parameter :: commwidth=10
c
c-- /proc/self/stat file data
      type statdata
       integer :: pid        !%d !The process ID.
       character(commwidth) :: comm !%s !The filename of the executable, in parentheses. This is visible whether or not the executable is swapped out.
       character :: state    !%c !One character from the string "RSDZTW" where R is running, S is sleeping in an interruptible wait, D is waiting in uninterruptible disk sleep, Z is zombie, T is traced or stopped (on a signal), and W is paging.
       integer :: ppid       !%d !The PID of the parent.
       integer :: pgrp       !%d !The process group ID of the process.
       integer :: session    !%d !The session ID of the process.
       integer :: tty_nr     !%d !The tty the process uses.
       integer :: tpgid      !%d !The process group ID of the process which currently owns the tty that the process is connected to.
       integer :: flags      !%lu !The kernel flags word of the process. For bit meanings, see the PF_* defines in <linux/sched.h>. Details depend on the kernel version.
       integer :: minflt     !%lu !The number of minor faults the process has made which have not required loading a memory page from disk.
       integer :: cminflt    !%lu !The number of minor faults that the process's waited-for children have made.
       integer :: majflt     !%lu !The number of major faults the process has made which have required loading a memory page from disk.
       integer :: cmajflt    !%lu !The number of major faults that the process's waited-for children have made.
       integer :: utime      !%lu !The number of jiffies that this process has been scheduled in user mode.
       integer :: stime      !%lu !The number of jiffies that this process has been scheduled in kernel mode.
       integer :: cutime     !%ld !The number of jiffies that this process's waited-for children have been scheduled in user mode. (See also times(2).)
       integer :: cstime     !%ld !The number of jiffies that this process's waited-for children have been scheduled in kernel mode.
       integer :: priority   !%ld !The standard nice value, plus fifteen. The value is never negative in the kernel.
       integer :: nice       !%ld !The nice value ranges from 19 (nicest) to -19 (not nice to others).
       integer :: nul        !%ld !This value is hard coded to 0 as a placeholder for a removed field.
       integer :: itrealvalue!%ld !The time in jiffies before the next SIGALRM is sent to the process due to an interval timer.
       integer :: starttime  !%lu !The time in jiffies the process started after system boot.
       integer :: vsize      !%lu !Virtual memory size in bytes.
       integer :: rss        !%ld !Resident Set Size: number of pages the process has in real memory, minus 3 for administrative purposes. This is just the pages which count towards text, data, or stack space. This does not include pages which have not been demand-loaded in, or which are swapped out.
      end type statdata
      type(statdata) :: stat
c
      integer :: istat
c
c-- /proc/self/statm file data
      integer,save :: pagesize=0
      integer :: vss_pages,rss_pages
c
*     integer :: i,j,n,m
*     character(300) :: line
c
*     character(30) :: words(24)
*     character(30),allocatable :: words(:)
      logical,save :: fail=.false.
c
c-- invalid by default
      mbsize = -1  !-- return invalid number if error
c
c-- quick exit
      if(fail) return
c
c-- read statm file
      if(pagesize==0) then
       open(4,file='/proc/self/statm',status='old',iostat=istat)
       if(istat==0) read(4,*,iostat=istat) vss_pages,rss_pages
       close(4)
c
c-- check status
       if(istat/=0) then
        fail = .true.
        call warn('memusg','read error: /proc/self/statm')
        return
       endif
      endif
c
c-- read stat file
      open(4,file='/proc/self/stat',status='old',iostat=istat)
      if(istat==0) read(4,*,iostat=istat) stat
      close(4)
c
c-- check status
      if(istat/=0) then
       fail = .true.
       call warn('memusg','read error: /proc/self/stat')
       return
      endif
c
c-- pagesize
      if(pagesize==0) pagesize = nint(stat%vsize/dble(vss_pages))
c
*!     write(6,*) 'line'
*!     write(6,*) line
*c-- count the number of words
*      i = 0
*      do n=1,len(line)
*       j = scan(trim(line(i+1:)),' ')
*       if(j==0) exit
*       i = i + j
*      enddo
*!     write(6,*) '#words:',n
*c-- split line into words
*      allocate(words(n))
*      i = 0
*      do m=1,n-1
*       j = scan(trim(line(i+1:)),' ')
*       words(m) = line(i+1:i+j)
*       i = i + j
*      enddo
*      words(n) = line(i+1:)
*!     write(6,*) 'words'
*!     write(6,'(a)') (words(i),i=1,n)
*c-- read stat from words
*      read(words,*) stat
*      deallocate(words)
*!     write(6,*) 'stat'
*!     write(6,*) stat
c
      mbsize(1) = (stat%rss*pagesize)/1024**2
      mbsize(2) = (stat%vsize)/1024**2
c
      end function memusg
