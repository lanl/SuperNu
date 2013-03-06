      subroutine warn(source,mesg,sunt)
c     ---------------------------------
      implicit none
      character(*),intent(in) :: source
      character(*),intent(in) :: mesg
      character(*),intent(in),optional :: sunt
************************************************************************
* write warning message to multiple file units
************************************************************************
      integer :: i,iunt
      integer,save :: nwarn=0 !warning counter
      character(27),parameter :: fmt='(" WARNING #",i0,2(": ",a))'
c
      nwarn = nwarn + 1 !update counter
c
      if(.not. present(sunt)) then
       write(0,fmt) nwarn,source,mesg
       write(6,fmt) nwarn,source,mesg
      else
       do i=1,len(sunt)
        read(sunt(i:i),*) iunt
        write(iunt,fmt) nwarn,source,mesg
       enddo
      endif
      end subroutine warn
