      subroutine read_restart_file
c     ----------------------------
      use gasgridmod
************************************************************************
* read restart file
************************************************************************
      character(13) :: fname = 'input.restart'
      logical :: lexist
c
      inquire(file=fname,exist=lexist)
      if(.not.lexist) return
c
      open(unit=4,file=fname,status='unknown')
      read(4,*) gas_vals2%temp
      close(4)
      end subroutine read_restart_file
c
c
c
      subroutine write_restart_file
c     -----------------------------
      use gasgridmod
************************************************************************
* write restart file
************************************************************************
      character(14) :: fname = 'output.restart'
      open(unit=4,file=fname,status='unknown',position='append')
      write(4,*) gas_vals2%temp
      close(4)
      end subroutine write_restart_file
