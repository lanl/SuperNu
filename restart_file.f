      subroutine read_restart_file
c     ----------------------------
      use gasgridmod
************************************************************************
* read restart file
************************************************************************
      !open(unit=4,file='Tinit.dat',status='unknown')
      open(unit=4,file='input.restart',status='unknown')
      read(4,*) gas_vals2%tempkev
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
      open(unit=4,file='output.restart',status='unknown')
      write(4,*) gas_vals2%tempkev
      close(4)
      end subroutine write_restart_file
