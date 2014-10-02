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
      real*8,allocatable :: str_velleft(:) !(nr+1)
      real*8,allocatable :: str_mass(:) !(nr)
      real*8,allocatable :: str_massfr(:,:) !(nabund,nr)
c
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
      use gasgridmod, only:gas_ini56,gas_ico56
      use miscmod
      implicit none
      integer,intent(in) :: nr
      real*8,intent(out) :: velout
************************************************************************
* Read the input structure file
************************************************************************
      integer :: ierr,nr_r,ini56
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
      allocate(str_velleft(nr+1))
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
c-- translate to velleft
      str_velleft(2:) = str_velright
      str_velleft(1) = 0d0
c
c-- close file
      close(4)
      deallocate(raw)
c
c-- result
      velout = str_velright(nr)
c
c-- convert abundlabl to element codes
      call elnam2elcode(ini56)
c
c-- output
      write(6,*)
      write(6,*) 'input structure:'
      write(6,*) '================'
      write(6,*) 'mass  :', sum(str_mass)/pc_msun, 'Msun'
c-- ni56 mass
      if(ini56>0) then
       help = sum(str_massfr(ini56,:)*str_mass)
      else
       help = 0d0
      endif
      write(6,*) 'm_ni56:', help/pc_msun, 'Msun'
c
      end subroutine read_inputstr
c
c
c
      subroutine elnam2elcode(ini56)
c     ------------------------------
      use miscmod, only:lcase
      use gasgridmod, only:gas_ini56,gas_ico56
      use elemdatamod
      implicit none
      integer,intent(out) :: ini56
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
c-- default
      ini56 = 0
c
c-- determine atomic z number
      do i=1,str_nabund
       iabund = 0
       elname = lcase(trim(str_abundlabl(i)))
       if(elname=='ni56') then
c-- special case
        iabund = gas_ini56
        ini56 = i
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
c
c
c-- use input.par variables to populate mass and grid arrays
c-- WARNING: size (nr+1) allocated for str_velright in generate_inputstr
      subroutine generate_inputstr(lr,velout)
      use inputparmod
      implicit none
      real*8,intent(out) :: lr,velout
************************************************************************
* generate stratification from input.par variables
* if in_noreadstruct==.true.
************************************************************************
      real*8,allocatable :: rout(:) !(nr+1)
      integer :: ir
      real*8 :: help, help2, dr
c
c-- verifications (input.par)
      if(in_velout<=0d0.and.in_isvelocity)
     &  stop 'generate_inputstr: invalid in_velout'
      if(in_lr<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr: invalid in_lr'
      if(in_totmass<0d0)
     &  stop 'generate_inputstr: invalid in_totmass'
c
c-- allocate arrays
      allocate(rout(in_nr+1))
      allocate(str_velright(in_nr))
      allocate(str_velleft(in_nr+1))
      allocate(str_mass(in_nr))
c
c-- local copies
      velout = in_velout
      lr = in_lr
c
c-- create unit sphere radii rout
      if(in_isvelocity) then
       help2 = in_velout
      else
       help2 = in_lr
      endif
      dr = 1d0/real(in_nr)
      forall(ir=1:in_nr+1) rout(ir) = help+(ir-1)*dr
c
c-- outer shells
      str_velleft = help2*rout
      str_velright = str_velleft(2:)
c
c-- mass
      if(in_dentype=='unif') then
       str_mass = in_totmass*(rout(2:)**3 - rout(:in_nr)**3)
       str_mass = str_mass/(1d0 - rout(1)**3)
      elseif(in_dentype=='mass') then
       forall(ir=1:in_nr)str_mass(ir) = in_totmass/real(in_nr)
      else
       stop 'generate_inputstr: invalid in_dentype'
      endif
c
      end subroutine generate_inputstr
c
      end module inputstrmod
