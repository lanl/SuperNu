*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module bbxsmod
c     --------------
      implicit none
************************************************************************
* read bound-bound cross sections (oscillator strengths)
************************************************************************
      integer :: bb_nlevel=0
      integer :: bb_nline=0
c
c-- raw data read from file
c-- level data type
      type bb_raw_level_data
       real*4 :: chi     !in cm^-1
       integer :: id,g
      end type bb_raw_level_data
      type(bb_raw_level_data),allocatable :: bbxs_level(:) !bb_nlevel
c-- line data type
      type bb_raw_line_data
       integer :: lev1,lev2
       real*4 :: f
      end type bb_raw_line_data
      type(bb_raw_line_data),allocatable :: bbxs_line(:) !bb_nline
c
c-- permanent data
      type bb_xs_data
       real*4 :: wl0   !in ang
       real*4 :: gxs   !g*xs = g*f*(pi*e**2/m_e/c)
       real*4 :: chilw !exp(h*c*chi_low/(k*1e4))  ![chi] = cm^-1
       integer*2 :: iz
       integer*2 :: ii
      end type bb_xs_data
      type(bb_xs_data),allocatable :: bb_xs(:) !bb_nline
c
      save
c
      contains
c
c
c
      subroutine read_atom(iz,ii,istat,get_data)
c     ------------------------------------------!{{{
      use miscmod, only:lcase
      use elemdatamod, only:elem_data
      implicit none
      integer,intent(in) :: iz,ii
      integer,intent(out) :: istat
      logical,intent(in) :: get_data
************************************************************************
* Read a single .atom file, or just poll how many lines it contains.
************************************************************************
      character(80) :: word
      character(13) :: fname
c-- level id
      integer :: l,lidmax,istat2
      integer,allocatable :: lid(:)
      integer(1) :: byte
c
c-- filename
      write(fname,'("data.atom.",a,i1)')
     &  lcase(trim(elem_data(iz)%sym)),ii
      open(4,file='Atoms/'//trim(fname),status='old',action='read',
     &  iostat=istat)
      if(istat/=0) goto 66
c
c-- read data size
      read(4,*)
      read(4,*)
      read(4,*)
      read(4,*,iostat=istat) bb_nlevel,bb_nline,word
      if(bb_nlevel<=0 .or. bb_nline<=0) istat = 1
      if(istat/=0) goto 66
c-- poll ready, return
      if(.not. get_data) goto 67
c
c-- allocate
      allocate(bbxs_level(bb_nlevel),bbxs_line(bb_nline))
c
c-- read data
c-- level data
      read(4,'(f11.3,13x,2i5)',iostat=istat) bbxs_level
      if(istat/=0) then
       write(6,*) fname
       stop 'read_atom: level data read error'
      endif
c-- line data
      read(4,'(2i5,f7.3)',iostat=istat) bbxs_line
      if(istat/=0) then
       write(6,*) fname
       stop 'read_atom: level data read error'
      endif
c-- verify eof
      read(4,*,iostat=istat2) byte
      if(istat2>=0) then
       write(6,*) fname,byte
       stop 'read_atom: data remaining on input file'
      endif
c
c-- construct reverse level pointer
      lidmax = maxval(bbxs_level(:)%id)
      allocate(lid(lidmax))
      lid = 0
      do l=1,bb_nlevel
       lid(bbxs_level(l)%id) = l
      enddo !l
c-- fix level id
      do l=1,bb_nline
       bbxs_line(l)%lev1 = lid(bbxs_line(l)%lev1)
       bbxs_line(l)%lev2 = lid(bbxs_line(l)%lev2)
      enddo !l
      deallocate(lid)
c
67    continue
      close(4)
      return
c
66    continue
      if(.not. get_data) write(6,*) 'read_atom failed:',iz,ii
      close(4)
c!}}}
      end subroutine read_atom
c
c
c
      subroutine sort_lines
c     ---------------------!{{{
      implicit none
************************************************************************
* Sort the bound-bound transitions for wl in order to speed up the
* the insertion into the big opacity array.
*
* So far, this doesn't seem to speed-up filling the cap array
* significantly.
************************************************************************
      integer :: l,is,it
      real*8,allocatable :: wl(:)  !too big for the stack
      type(bb_xs_data) :: xs_src,xs_trg
      integer,allocatable :: indx(:),indx_inv(:)  !too big for the stack
c
c-- initialize arrays
      allocate(wl(bb_nline),indx(bb_nline),indx_inv(bb_nline))
      wl = dble(bb_xs%wl0)
      forall(l=1:bb_nline) indx(l) = l
c
c-- index sort
      call sorti(bb_nline,wl,indx)
      deallocate(wl)
c-- reverse indx
      forall(l=1:bb_nline) indx_inv(indx(l)) = l
c
c-- sort the big structure: move an element to its final position,
c-- shifting the value there to the temporary
      l = 2
      is = 1
      xs_src = bb_xs(is) !save original data in target
      do
c-- proceed to a new unsorted position if a closed subset was sorted
       if(indx_inv(is)==0) then
        do l=l,bb_nline
         if(indx_inv(l)/=0) exit
        enddo
        if(l>bb_nline) exit
        is = l             !source location
        xs_src = bb_xs(is) !save original data in target
       endif
c-- step 1: save target
       it = indx_inv(is)       !target location
       xs_trg = bb_xs(it)  !save original data in target
c-- step 2: move old source to sorted position
       bb_xs(it) = xs_src  !move source to target
       indx_inv(is) = 0        !mark source position ready
c-- step 3: rotate old target to new source
       is = it             !next source location
       xs_src = xs_trg     !rotate old target to new source
      enddo
c!}}}
      end subroutine sort_lines
c
      end module bbxsmod
c vim: fdm=marker
