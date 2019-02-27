*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module ffxsmod
c     --------------
      implicit none
************************************************************************
* free-free cross sections based on the free-free gaunt factor
* calculations of Sutherland, 1998, MNRAS, 300, 321.
*
* Cloned from the Chianti atomic database.
*
* Note that Sutherland's Eq.(15) has units of erg/cm^3/s. Comparing
* with Rybicki & Lightman's Eq.5.14(a) (in their book 'Radiative
* Processes in Astrophysics'), suggests that Sutherland's units
* should be erg/cm^3/s/sr/Hz. We are assuming the latter to be
* correct in this routine.
************************************************************************
      integer,parameter :: ff_nu=81, ff_ngg=41
c
      real*8 :: ff_gff(ff_nu,ff_ngg)
c
      save
c
      contains
c
c
c
      subroutine ffxs_read_data
c     -------------------------
      implicit none
************************************************************************
* Read sutherland data
************************************************************************
      character(18),parameter :: fname='data.ff_sutherland'
      integer :: istat
c
      real*8 :: gff_raw(3,ff_nu,ff_ngg)
c
c-- read
      open(4,file=fname,action='read',status='old',
     &  iostat=istat)
      if(istat/=0) stop 'ffxs_read_data: cannot read ff_sutherland.dat'
      read(4,'(4/)') !comments
      read(4,*,err=99) gff_raw
      close(4)
c
c-- save in permanent array
      ff_gff = gff_raw(3,:,:)
c
      return
c
 99   stop 'ffxs_read_data: read error'
c
      end subroutine ffxs_read_data
c
      end module ffxsmod
c vim: fdm=marker
