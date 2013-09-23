      subroutine read_restart_file
c     ----------------------------
      use gasgridmod
************************************************************************
* read restart file
************************************************************************
      character(13) :: fname = 'input.restart'
      integer :: istat
c
      open(unit=4,file=fname,status='old',iostat=istat)
      if(istat/=0) stop 'read_restart: no input.restart file'
      read(4,*) gas_temp
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
      write(4,*) gas_temp
      close(4)
      end subroutine write_restart_file
c
c
c
      subroutine read_restart_randcount
c     ---------------------------------
      use particlemod
************************************************************************
* read number of rand calls from each rank from restart file
* at some time step.
************************************************************************
      character(15) :: fname = 'input.retlyrand'
c
      subroutine read_restart_randcount
c
c
c
      subroutine write_restart_randcount
c     ----------------------------------
      use particlemod
************************************************************************
* write number of rand calls from each rank to restart file
* at each time step.
************************************************************************
      character(16) :: fname = 'output.retlyrand'
      open(unit=4,file=fname,status='unknown',position='append')
      write(4,*) prt_tlyrandarr
      close(4)
      subroutine write_restart_randcount
