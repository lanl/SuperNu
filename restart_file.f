      subroutine read_restart_file
c     ----------------------------
      use gridmod
      use timestepmod
      implicit none
************************************************************************
* read restart file
************************************************************************
      character(13) :: fname = 'input.restart'
      integer :: istat
      integer :: it
c
      open(unit=4,file=fname,status='old',iostat=istat)
      if(istat/=0) stop 'read_restart: no input.restart file'
c-- assumes no header
      do it = 1, tsp_ntres-2
         read(4,*)
      enddo
      read(4,*) grd_tempinv
      close(4)
      end subroutine read_restart_file
c
c
c
      subroutine write_restart_file
c     -----------------------------
      use gridmod
      implicit none
************************************************************************
* write restart file
************************************************************************
      character(14) :: fname = 'output.restart'
      open(unit=4,file=fname,status='unknown',position='append')
      write(4,*) grd_tempinv
      close(4)
      end subroutine write_restart_file
