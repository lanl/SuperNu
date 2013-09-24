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
      use timestepmod
************************************************************************
* read number of rand calls from each rank from restart file
* at some time step.
* output.retlyrand must be copied to input.retlyrand
************************************************************************
      character(15) :: fname = 'input.retlyrand'
      integer :: istat
      integer :: it
      integer, allocatable :: helprandarr(:)
c
      allocate(helprandarr(size(prt_tlyrandarr,2)))
      helprandarr = 0
c
      if(tsp_ntres>1) then
         open(unit=4,file=fname,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.retlyrand file'
         do it = 1:tsp_ntres
            read(4,*) prt_tlyrandarr
            helprandarr=helprandarr+prt_tlyrandarr
         enddo
         close(4)
      else
         if(tsp_ntres/=1) stop 'read_restar: tsp_ntres<1'
      endif
      prt_tlyrandarr = helprandarr
c
      end subroutine read_restart_randcount
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
      end subroutine write_restart_randcount
