      subroutine open_logfiles
c     ------------------------
      use inputparmod
      implicit none
************************************************************************
* Open output file handles.
************************************************************************
      integer :: istat
c
      if(in_io_grabstdout) then
       write(6,*) 'write stdout to output.log'
       open(6,file='output.log',action='write',status='replace',
     &   recl=5000,iostat=istat) !write stdout to file
       if(istat/=0) stop 'parse_inputpars: open output.log error'
       call banner
      endif
c
      open(8,file='output.logdata',action='write',status='replace',
     &  recl=5000,iostat=istat) !write stdout to file
      if(istat/=0) stop 'parse_inputpars: open output.data error'
c
      end subroutine open_logfiles
