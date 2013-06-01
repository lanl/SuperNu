      module inputstrmod
c     ------------------
      implicit none
************************************************************************
* Supernova atmospheric stratification
************************************************************************
      character(9),private :: fname='input.str'
      integer :: str_nabund=0
      integer,allocatable :: str_iabund(:) !(nabund)
c
      real*8,allocatable :: str_velright(:) !(nr)
      real*8,allocatable :: str_mass(:) !(nr)
      real*8,allocatable :: str_massfr(:,:) !(nabund,nr)
      character(8),allocatable :: str_abundlabl(:) !(nabund)
c
      save
      public
c
      contains
c
      subroutine read_inputstr(nr,velout)
c     -----------------------------------
      use physconstmod
      use miscmod
      implicit none
      integer,intent(in) :: nr
      real*8,intent(out) :: velout
************************************************************************
* Read the input structure file
************************************************************************
      integer :: ierr,nr_r
      character(2) :: dmy
      character(8) :: labl(2)
      real*8,allocatable :: raw(:,:)
      real*8 :: help
c
c-- open file
      open(4,file=fname,status='old',iostat=ierr)
      if(ierr/=0) stop 'read_inputstr: file missing: input.str'
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy,nr_r,str_nabund
      if(ierr/=0) stop 'read_inputstr: input.str format err: dimensions'
c-- verify dimension
      if(nr_r/=nr) stop 'read_inputstr: incompatible nr dimension'
c
c-- allocate arrays
      allocate(str_velright(nr))
      allocate(str_mass(nr))
      allocate(str_massfr(str_nabund,nr))
      allocate(str_abundlabl(str_nabund))
      allocate(raw(str_nabund+2,nr))
c
c-- read labels
      read(4,*,iostat=ierr) dmy, labl, str_abundlabl
      if(ierr/=0) stop 'read_inputstr: input.str format err: col labels'
c
c-- read body
      read(4,*,iostat=ierr) raw
      if(ierr/=0) stop 'read_inputstr: input.str format err: body'
c
c-- transer data to final arrays
      str_velright = raw(1,:)
      str_mass = raw(2,:)
      str_massfr = raw(3:,:)
c
c-- close file
      close(4)
c
c-- result
      velout = str_velright(nr)
c
c-- convert abundlabl to element codes
      call elnam2elcode
c
c-- output
      write(6,*)
      write(6,*) 'input structure:'
      write(6,*) '================'
      write(6,*) 'mass  :', sum(str_mass)/pc_msun
      help = sum(str_massfr(str_nabund,:)*str_mass)
      write(6,*) 'm_ni56:', help/pc_msun
c
      end subroutine read_inputstr
c
c
c
      subroutine elnam2elcode
c     ------------------------------
      use miscmod, only:lcase
      use gasgridmod, only:gas_ini56,gas_ico56
      use elemdatamod
      implicit none
************************************************************************
* convert the abundlabl labels to element codes (atomic z number), which
* also serve as indexing pointers to the gas_vals2%mass0fr array.
************************************************************************
      integer :: i,j,iabund
      character(4) :: elname
c
c-- allocate element code array (pointer to gas_vals2%mass0fr)
      allocate(str_iabund(str_nabund))
c
c-- determine atomic z number
      do i=1,str_nabund
       iabund = 0
       elname = lcase(trim(str_abundlabl(i)))
       if(elname=='ni56') then
c-- special case
        iabund = gas_ini56
       elseif(elname=='co56') then
c-- special case
        iabund = gas_ico56
       else
c
c-- normal case, search elem_data for corresponding entry
        do j=1,elem_neldata
         if(lcase(trim(elem_data(j)%sym))==elname) exit
        enddo
c-- verify hit
        if(j>elem_neldata) then
         write(*,*) 'unknown chemical element name:',elname
         stop 'elnam2elcode: no such element found in elemdata'
        endif
        iabund = j
       endif
c
c-- store element code (pointer to gas_vals2%mass0fr)
       str_iabund(i) = iabund
      enddo
c
      end subroutine elnam2elcode
c
      end module inputstrmod
