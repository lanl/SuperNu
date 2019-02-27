*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine output_grid
c     ----------------------
      use inputparmod
      use timingmod
      use timestepmod
      use gridmod
      use totalsmod
      implicit none
************************************************************************
* Write grid values to file
************************************************************************
      integer :: i,j,k
      integer :: np,npbot
      integer,save :: nrow=0 !number of rows for printing grid variables
      integer,save :: ncpr=0 !number of cells per rank
      integer :: reclen
      logical,save :: lfirst=.true.
      character(16),save :: fstat='replace'
c
      integer,allocatable :: iarr(:)
      real*8,allocatable :: arr(:)
      real*8 :: t0,t1
c
      t0 = t_time()
c
      if(lfirst) then
c
c-- shape of volume quantities
       if(grd_ivoid>0) then
c-- no void cells
        ncpr = grd_nx
       else
c-- find optimal row size
        npbot = grd_ncell !init with worst case
        do i=grd_nx,max(grd_nx,grd_ny,grd_nz,grd_ncell/2)+1
           nrow = ceiling(float(grd_ncell)/i)
           np = nrow*i - grd_ncell !number of pad cells
           if(np<npbot) then
              npbot = np
              ncpr = i !optimal row size
           endif
        enddo
       endif
       nrow = ceiling(float(grd_ncell)/ncpr)
c
c-- write once
c=============
       reclen = max(4,grd_nx+1,grd_ny+1,grd_nz+1,grd_ny*grd_nz)*12
       open(unit=4,file='output.grd_grid',status='replace',recl=reclen)
c-- header: dimension
       write(4,*) "#",grd_igeom
       write(4,*) "#",grd_nx,grd_ny,grd_nz
       write(4,*) "#",nrow*ncpr,nrow,ncpr
c-- body
       write(4,'(1p,10000e12.4)') grd_xarr(:)
       write(4,'(1p,10000e12.4)') grd_yarr(:)
       write(4,'(1p,10000e12.4)') grd_zarr(:)
c-- cell indices
       do k=1,grd_nz
       do j=1,grd_ny
          write(4,'(10000i12)') grd_icell(:,j,k)
       enddo
       enddo
       close(4)
      endif

c
c-- scalars
c==========
      open(unit=4,file='output.tsp_time',status=fstat,position='append')
      if(lfirst) write(4,*) "#",tsp_nt
      if(lfirst) write(4,*) tsp_t  !beginning of first time step
      write(4,*) tsp_t + tsp_dt  !end of time step
      close(4)

      open(unit=4,file='output.tot_energy',status=fstat,
     &  position='append',recl=17*12)
      if(lfirst) then
       write(4,'("#",i5)') 16 !number of columns
       write(4,'("#",16a12)')
     &   'eerror','erad','emat','eext','eout','evelo',
     &   'sfluxgamma','sflux',
     &   'sthermal','smanufac','sanalvol','sanalsurf','samp',
     &   'sdecaygamma','sdecaybeta','sdeposgamma'
      endif
      write(4,'(1x,1p,16e12.4)')
     &  tot_eerror,tot_erad,tot_emat,tot_eext,tot_eout,tot_evelo,
     &  tot_sfluxgamma,tot_sflux,
     &  tot_sthermal,tot_smanufac,tot_sanalvol,tot_sanalsurf,tot_samp,
     &  tot_sdecaygamma,tot_sdecaybeta,tot_sdeposgamma
      close(4)

c
c-- grid arrays
c==============
      if(.not.in_io_nogriddump) then
       reclen = ncpr*12
c-- alloc
       allocate(iarr(ncpr*nrow),arr(ncpr*nrow))
       iarr(grd_ncell+1:) = 0
       arr(grd_ncell+1:) = 0d0
c
c-- grdinput values
       iarr(:grd_ncell) = grd_nvol
       open(unit=4,file='output.grd_nvol',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = 1d0/grd_tempinv
       open(unit=4,file='output.grd_temp',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = grd_fcoef
       open(unit=4,file='output.grd_fcoef',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = grd_capgrey
       open(unit=4,file='output.grd_capgrey',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = grd_sig
       open(unit=4,file='output.grd_sig',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
       arr(:grd_ncell) = grd_tally(2,:)/grd_vol
       open(unit=4,file='output.grd_eraddens',status=fstat,
     &   position='append',recl=reclen)
       do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
       enddo
       close(4)
c
c-- grdtally values
       if(in_io_dogrdtally) then
        iarr(:grd_ncell) = grd_methodswap
        open(unit=4,file='output.grd_methodswap',status=fstat,
     &   position='append',recl=reclen)
        do i=1,nrow
         write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
c
        iarr(:grd_ncell) = grd_numcensimc
        open(unit=4,file='output.grd_numcensimc',status=fstat,
     &   position='append',recl=reclen)
        do i=1,nrow
         write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
c
        iarr(:grd_ncell) = grd_numcensddmc
        open(unit=4,file='output.grd_numcensddmc',status=fstat,
     &   position='append',recl=reclen)
        do i=1,nrow
         write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
       endif !in_io_dogrdtally
c
       deallocate(iarr,arr)
      endif

c
c-- after the first iteration open files in append mode
      lfirst = .false.
      fstat = 'old'
c
c-- timing
      t1 = t_time()
      call timereg(t_output, t1-t0)
c
      end subroutine output_grid
c vim: fdm=marker
