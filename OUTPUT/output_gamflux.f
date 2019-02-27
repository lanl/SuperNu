*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine output_gamflux
c     -------------------------
      use timestepmod
      use fluxmod
      implicit none
************************************************************************
* Write flux arrays to file.
************************************************************************
      integer :: j,k
      logical,save :: lfirst=.true.
      character(16),save :: fstat='replace'
      integer :: reclen
c
c-- gamflux arrays
c==============
      reclen = flx_ng*12
      open(unit=4,file='output.flx_gamluminos',status=fstat,
     &  position='append')
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(1p,e12.4)') merge(flx_gamluminos(j,k),0d0,
     &   mask=flx_gamluminos(j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
      enddo
      enddo
      close(4)

      open(unit=4,file='output.flx_gamlumnum',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(i12)') flx_gamlumnum(j,k)
      enddo
      enddo
      close(4)

      open(unit=4,file='output.flx_gamlumtime',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(e12.4)') merge(flx_gamlumtime(j,k),0d0,
     &   mask=flx_gamlumtime(j,k)>1d-99)
      enddo
      enddo
      close(4)

      open(unit=4,file='output.flx_gamlumdev',status=fstat,
     &  position='append',recl=reclen)
      do k=1,flx_nom
      do j=1,flx_nmu
       write(4,'(1p,e12.4)') merge(flx_gamlumdev(j,k),0d0,
     &   mask=flx_gamlumdev(j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
      enddo
      enddo
      close(4)

c
c-- after the first iteration open files in append mode
      lfirst = .false.
      fstat = 'old'
c
      end subroutine output_gamflux
c vim: fdm=marker
