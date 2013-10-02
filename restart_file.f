      subroutine read_restart_file
c     ----------------------------
      use gasgridmod
      use timestepmod
************************************************************************
* read restart file
************************************************************************
      character(13) :: fname = 'input.restart'
      integer :: istat
      integer :: it
c
      open(unit=4,file=fname,status='old',iostat=istat)
      if(istat/=0) stop 'read_restart: no input.restart file'
      if(tsp_ntres<=1) then
         read(4,*) gas_temp
      else
c
c-- assumes no header
         do it = 1, tsp_ntres-2
            read(4,*)
         enddo
         read(4,*) gas_temp
      endif
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
      character(13) :: fname = 'input.tlyrand'
      integer :: istat
      integer :: it
      integer, allocatable :: helprandarr(:)
c
      allocate(helprandarr(size(prt_tlyrandarr)))
      helprandarr = 0
c
      if(tsp_ntres>1) then
         open(unit=4,file=fname,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.retlyrand file'
         do it = 1, tsp_ntres
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
      character(14) :: fname = 'output.tlyrand'
      open(unit=4,file=fname,status='unknown',position='append')
      write(4,*) prt_tlyrandarr
      close(4)
      end subroutine write_restart_randcount
c
c
c
      subroutine read_restart_particles
c     ----------------------------------
      use particlemod
************************************************************************
* read particle from each rank to restart file
* at some time step
************************************************************************
c

      end subroutine read_restart_particles
c
c
c
      subroutine write_restart_particles
c     ----------------------------------
      use particlemod
************************************************************************
* write particle from each rank to restart file
* at each time step
************************************************************************
c
C$$$      character(14) :: fnamez = 'output.prtzsrc'
C$$$      character(15) :: fnamert = 'output.prtrtsrc'
C$$$      character(14) :: fnamer = 'output.prtrsrc'
C$$$      character(15) :: fnamemu = 'output.prtmusrc'
C$$$      character(14) :: fnamet = 'output.prttsrc'
C$$$      character(14) :: fnamee = 'output.prtesrc'
C$$$      character(16) :: fnameeb = 'output.prtebirth'
C$$$      character(15) :: fnamewl = 'output.prtwlsrc'
      character(18) :: fnamevac = 'output.prtisvacant'
c
c-- write particle array vacancies
      open(unit=4,file=fnamevac,status='unknown',position='rewind')
      write(4,*) prt_tlyvacant(1:10,1)
      write(4,*) prt_tlyvacant(1:10,2)
      close(4)
c
      end subroutine write_restart_particles
