      module bbxsmod
c     --------------
      implicit none
************************************************************************
* read bound-bound cross sections (oscillator strengths)
************************************************************************
      integer :: bb_nlevel,bb_nline
c
c-- raw data read from file
c-- level data type
      type bb_raw_level_data
       real :: chi     !in cm^-1
       integer :: id,g
      end type bb_raw_level_data
      type(bb_raw_level_data),allocatable :: bbxs_level(:) !bb_nlevel
c-- line data type
      type bb_raw_line_data
       integer :: lev1,lev2
       real :: f
      end type bb_raw_line_data
      type(bb_raw_line_data),allocatable :: bbxs_line(:) !bb_nline
c
c-- permanent data
      type bb_xs_data
       sequence
       integer*2 :: iz
       integer*2 :: ii
       real :: wl0    !in ang
       real :: gxs    !g*xs = g*f*(pi*e**2/m_e/c)
       real :: chilw !exp(h*c*chi_low/(k*1e4))  ![chi] = cm^-1
      end type bb_xs_data
      type(bb_xs_data),allocatable :: bb_xs(:)
c
      contains
c
c
c
      subroutine read_atom(iz,ii,istat,get_data)
c     ------------------------------------------
c     use physconstmod
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
      character(8) :: fname
c-- level id
      integer :: i,lidmax,istat2
      integer,allocatable :: lid(:)
      integer(1) :: byte
c
c-- filename
      write(fname,'(a2,i1,".atom")') lcase(trim(elem_data(iz)%sym)),ii
      open(4,file='Atom/'//adjustl(fname),status='old',action='read',
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
      do i=1,bb_nlevel
       lid(bbxs_level(i)%id) = i
      enddo !i
c-- fix level id
      do i=1,bb_nline
       bbxs_line(i)%lev1 = lid(bbxs_line(i)%lev1)
       bbxs_line(i)%lev2 = lid(bbxs_line(i)%lev2)
      enddo !i
      deallocate(lid)
c
67    continue
      close(4)
      return
c
66    continue
      if(.not. get_data) write(6,*) 'read_atom failed:',iz,ii
      close(4)
c
      end subroutine read_atom
c
c
c
      subroutine sort_lines
c     ---------------------
      implicit none
************************************************************************
* Sort the bound-bound transitions for wl in order to speed up the
* the insertion into the big opacity array.
*
* So far, this doesn't seem to speed-up filling the gas_cap array
* significantly.
************************************************************************
      integer :: i,is,it
      type(bb_xs_data) :: xs_src,xs_trg
      real*8 :: wl(bb_nline)
      integer :: indx(bb_nline)
c
c-- initialize arrays
      wl = dble(bb_xs%wl0)
      forall(i=1:bb_nline) indx(i) = i
c
c-- index sort
      call sorti(bb_nline,wl,indx)
c-- reverse indx
      forall(i=1:bb_nline) indx(indx(i)) = i
c
c-- sort the big structure: move an element to its final position,
c-- shifting the value there to the temporary
      i = 2
      is = 1
      xs_src = bb_xs(is) !save original data in target
      do
c-- proceed to a new unsorted position if a closed subset was sorted
       if(indx(is)==0) then
        do i=i,bb_nline
         if(indx(i)/=0) exit
        enddo
        if(i>bb_nline) exit
        is = i             !source location
        xs_src = bb_xs(is) !save original data in target
       endif
c-- step 1: save target
       it = indx(is)       !target location
       xs_trg = bb_xs(it)  !save original data in target
c-- step 2: move old source to sorted position
       bb_xs(it) = xs_src  !move source to target
       indx(is) = 0        !mark source position ready
c-- step 3: rotate old target to new source
       is = it             !next source location
       xs_src = xs_trg     !rotate old target to new source
      enddo
c
      end subroutine sort_lines
      end module bbxsmod
