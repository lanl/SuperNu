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
         read(4,*) grd_temp
      else
c
c-- assumes no header
         do it = 1, tsp_ntres-2
            read(4,*)
         enddo
         read(4,*) grd_temp
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
      write(4,*) grd_temp
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
         if(istat/=0) stop 'read_restart: no input.tlyrand file'
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
      character(13) :: fnamez = 'input.prtzsrc'
      character(14) :: fnamert = 'input.prtrtsrc'
      character(13) :: fnamer = 'input.prtrsrc'
      character(14) :: fnamemu = 'input.prtmusrc'
      character(13) :: fnamet = 'input.prttsrc'
      character(13) :: fnamee = 'input.prtesrc'
      character(15) :: fnameeb = 'input.prtebirth'
      character(14) :: fnamewl = 'input.prtwlsrc'
      character(17) :: fnamevac = 'input.prtisvacant'
c
c-- reading particle vacancies
         open(unit=4,file=fnamevac,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyvacant file'
         read(4,*) prt_tlyvacant
         close(4)
c
c-- reading particle zones
         open(unit=4,file=fnamez,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyzsrc file'
         read(4,*) prt_tlyzsrc
         close(4)
c
c-- reading particle transport index
         open(unit=4,file=fnamert,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyrtsrc file'
         read(4,*) prt_tlyrtsrc
         close(4)
c
c-- reading particle positions
         open(unit=4,file=fnamer,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyrsrc file'
         read(4,*) prt_tlyrsrc
         close(4)
c
c-- reading particle directions
         open(unit=4,file=fnamemu,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlymusrc file'
         read(4,*) prt_tlymusrc
         close(4)
c
c-- reading particle times
         open(unit=4,file=fnamet,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlytsrc file'
         read(4,*) prt_tlytsrc
         close(4)
c
c-- reading particle energies
         open(unit=4,file=fnamee,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyesrc file'
         read(4,*) prt_tlyesrc
         close(4)
c
c-- reading particle birth energies
         open(unit=4,file=fnameeb,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlyebirth file'
         read(4,*) prt_tlyebirth
         close(4)
c
c-- reading particle wavelengths
         open(unit=4,file=fnamewl,status='old',iostat=istat)
         if(istat/=0) stop 'read_restart: no input.tlywlsrc file'
         read(4,*) prt_tlywlsrc
         close(4)
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
      character(14) :: fnamez = 'output.prtzsrc'
      character(15) :: fnamert = 'output.prtrtsrc'
      character(14) :: fnamer = 'output.prtrsrc'
      character(15) :: fnamemu = 'output.prtmusrc'
      character(14) :: fnamet = 'output.prttsrc'
      character(14) :: fnamee = 'output.prtesrc'
      character(16) :: fnameeb = 'output.prtebirth'
      character(15) :: fnamewl = 'output.prtwlsrc'
      character(18) :: fnamevac = 'output.prtisvacant'
c
      call compress_particle_output
c-- writing particle vacancies
      open(unit=4,file=fnamevac,status='unknown',position='rewind')
      write(4,*) prt_tlyvacant
      close(4)
c
c-- writing particle zones
      open(unit=4,file=fnamez,status='unknown',position='rewind')
      write(4,*) prt_tlyzsrc
      close(4)
c
c-- writing particle transport index
      open(unit=4,file=fnamert,status='unknown',position='rewind')
      write(4,*) prt_tlyrtsrc
      close(4)
c
c-- writing particle positions
      open(unit=4,file=fnamer,status='unknown',position='rewind')
      write(4,*) prt_tlyrsrc
      close(4)
c
c-- writing particle directions
      open(unit=4,file=fnamemu,status='unknown',position='rewind')
      write(4,*) prt_tlymusrc
      close(4)
c
c-- writing particle times
      open(unit=4,file=fnamet,status='unknown',position='rewind')
      write(4,*) prt_tlytsrc
      close(4)
c
c-- writing particle energies
      open(unit=4,file=fnamee,status='unknown',position='rewind')
      write(4,*) prt_tlyesrc
      close(4)
c
c-- writing particle birth energies
      open(unit=4,file=fnameeb,status='unknown',position='rewind')
      write(4,*) prt_tlyebirth
      close(4)
c
c-- writing particle wavelength
      open(unit=4,file=fnamewl,status='unknown',position='rewind')
      write(4,*) prt_tlywlsrc
      close(4)
c
c
c

c
c
      contains
c
      subroutine compress_particle_output
************************************************************************
* Resizing arrays to only write non-vacant particle histories.
* Should reduce read/write work.
* May be advantageous for many ranks or few source particles per rank.
************************************************************************
c
c
      end subroutine compress_particle_output
c
      end subroutine write_restart_particles
