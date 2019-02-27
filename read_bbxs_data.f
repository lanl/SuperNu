*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine read_bbxs_data(nelem)
c     --------------------------------
      use physconstmod
      use ionsmod
      use bbxsmod
      use timingmod
      use miscmod, only:warn
      implicit none
      integer,intent(in) :: nelem
************************************************************************
* read all bound-bound cross sections
************************************************************************
      real,parameter :: fconst = sngl(pc_pi*pc_e**2/(pc_me*pc_c))
      integer :: nlinall,ilinall
      integer :: l,iz,ii,istat,llw,lhg,n
      real*8 :: t0,t1
c
c-- quick exit
      if(nelem==0) then
       call warn('read_bbxs_data','nelem==0: no lines read')
       return
      endif
c
c-- determine total number of lines
      nlinall = 0
      do iz=1,nelem
       do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
        call read_atom(iz,ii,istat,get_data=.false.)
c
c-- test if succesfull
        if(istat/=0) cycle
        write(8,*) 'read_atom successful:',iz,ii,bb_nlevel,bb_nline
        nlinall = nlinall + bb_nline
       enddo !ii
      enddo !iz
      write(6,'(a,i8)') ' total number of lines read:',nlinall
      write(8,'(a,i8)') ' total number of lines read:',nlinall
c
c-- allocate permanent storage space for line data
      if(nlinall<=0) stop 'rd_bbxs_data: no sigle line read in'
      allocate(bb_xs(nlinall))
      n = int(sizeof(bb_xs)/1024) !kB
      write(6,*) 'ALLOC bb_xs    :',n,"kB",n/1024,"MB",n/1024**2,"GB"
c
      t0 = t_time()
      ilinall = 0
      do iz=1,nelem
       do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
        call read_atom(iz,ii,istat,get_data=.true.)
        if(istat/=0) cycle
c
c-- store data in permanent array
        do l=1,bb_nline
         ilinall = ilinall + 1
         llw = bbxs_line(l)%lev1
         lhg = bbxs_line(l)%lev2
c-- line center wavelength
         bb_xs(ilinall)%wl0 = 1e8/(abs(bbxs_level(lhg)%chi) - !in ang
     &     abs(bbxs_level(llw)%chi))
c-- flip low<->high levels
         if(bb_xs(ilinall)%wl0 < 0.) then
          llw = lhg
          lhg = bbxs_line(l)%lev1
          bb_xs(ilinall)%wl0 = -bb_xs(ilinall)%wl0
         endif
c-- g*xs
         bb_xs(ilinall)%gxs = fconst*bbxs_level(llw)%g*
     &     10.**bbxs_line(l)%f          !fconst = pi*e**2/(m_e*c)
c-- exp(chi)
         bb_xs(ilinall)%chilw = abs(bbxs_level(llw)%chi)
c-- ion code
         bb_xs(ilinall)%iz = int(iz,2)
         bb_xs(ilinall)%ii = int(ii,2)
        enddo !l
c-- ready with raw data
        deallocate(bbxs_level,bbxs_line)
       enddo !ii
      enddo !iz
c-- verify counter
      if(ilinall/=nlinall) stop 'rd_bbxs_data: line counter error'
c
c-- store counter globally
      bb_nline = nlinall
c
c-- sort lines - doesn't speed-up bb opacity.
      call sort_lines
c
      t1 = t_time()
!     write(6,'(a,f8.2,a)') ' time used for bbxs reading:',t1-t0,'s'

      end subroutine read_bbxs_data
c vim: fdm=marker
