      function memusg() result(mbrss)
c     -------------------------
      use miscmod
      implicit none
      integer :: mbrss
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
c-- /proc/self/statm file data
      integer :: vss_pages,rss_pages
      integer,save :: pagesize=0
c
*     integer :: i,j,n,m
*     character(300) :: line
c
      character(30) :: words(24)
*     character(30),allocatable :: words(:)
      logical :: exists
      logical,save :: skipwarn=.false.
c
c-- init
      mbrss = -1  !-- return invalid number if error
c
c-- read stat file
      inquire(file='/proc/self/stat',exist=exists)
      if(.not.exists .and. .not.skipwarn) then
       call warn('memusg','file missing: /proc/self/stat,no mem info')
       skipwarn = .true.
       return
      endif
c-- read
      open(4,file='/proc/self/stat',action='read',status='old')
      read(4,*,end=6) words !not all compilers read words shorter than len(words(i))from the file
*      read(4,'(a)',end=6) line
      close(4)
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
c-- read statm file
      if(pagesize==0) then
       inquire(file='/proc/self/statm',exist=exists)
       if(.not.exists .and. .not.skipwarn) then
        call warn('memusg','file missing:/proc/self/stat,no mem info')
        skipwarn = .true.
        return
       endif
c-- read
       open(4,file='/proc/self/statm',action='read',status='old')
       read(4,*) vss_pages,rss_pages
       close(4)
c
c-- page size
       pagesize = nint(stat%vsize/dble(vss_pages))
      endif !pagesize
c
      mbrss = stat%rss*pagesize/1024**2
c
      return
c
 6    call warn('memusg','read error: eof: /proc/self/stat')
      close(4)
      skipwarn = .true.
c
      end function memusg
