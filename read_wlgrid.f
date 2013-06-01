      subroutine read_wlgrid
c     ----------------------
      use gasgridmod
      use inputparmod
      implicit none
************************************************************************
* read wavelength grid from file
*
* EXAMPLE FILE LAYOUT:
* --------------------
* #wavelength bin boundaries. units: [cm]
* #ncell ngroupmax
* #icell ngroup wlbound
* 10 5
*  1 1 1e-5 32e-5
*  2 1 1e-5 32e-5
*  3 1 1e-5 32e-5
*  4 1 1e-5 32e-5
*  5 1 1e-5 32e-5
*  6 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
*  7 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
*  8 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
*  9 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
* 10 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
* 
************************************************************************
c
      real*8, allocatable :: wlstore(:) 
      integer :: ng,ngm,nrm,irr,ir
c
      open(4,file='input.wlgrid',status='old')
c
c-- strip header
      do ir = 1, 3
         read(4,*)
      enddo
c
c-- read dimensions
      read(4,*) nrm, ngm
      write(6,*) 'maximume cell number, group number: ',nrm, ngm
c
      do ir = 1,in_wldex-1
         read(4,*)
      enddo
      read(4,*) irr, ng
c
      allocate(gas_wl(ng+1))
      allocate(wlstore(ng+3))
      rewind(4)
      do ir = 1, in_wldex+3
         read(4,*)
      enddo
      read(4,*) wlstore(:)
      close(4)
      gas_wl = wlstore(3:)
      deallocate(wlstore)
c
      write(6,*)
      write(6,*) 'wavelength grid in [cm]'
      write(6,*) gas_wl
      !Resetting in_ng to value of particular group entry from file
      in_ng = ng
c
      end subroutine read_wlgrid
