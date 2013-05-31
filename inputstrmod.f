      module inputstrmod
c     ------------------
      implicit none
************************************************************************
* Supernova atmospheric stratification
************************************************************************
      character(9),private :: fname='input.str'
      integer :: str_nabund=0
c
      real*8,allocatable :: str_vel(:) !(nr)
      real*8,allocatable :: str_mass(:) !(nr)
      real*8,allocatable :: str_massfr(:,:) !(nabund,nr)
      character(8),allocatable :: str_abundlabl(:) !(nabund)
c
      save
      public
c
      contains
c
      subroutine read_inputstr(nr)
c     ----------------------------
      use physconstmod
      use miscmod
      implicit none
      integer,intent(in) :: nr
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
      open(4,file=fname,status='old')
      if(ierr/=0) then
       call warn('read_inputstr','file missing: input.str')
       return
      endif
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy,nr_r,str_nabund
      if(ierr/=0) stop 'read_inputstr: input.str format err: dimensions'
c-- verify dimension
      if(nr_r/=nr) error stop 'read_inputstr: incompatible nr dimension'
c
c-- allocate arrays
      allocate(str_vel(str_nabund))
      allocate(str_mass(str_nabund))
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
      str_vel = raw(1,:)
      str_mass = raw(2,:)
      str_massfr = raw(3:,:)
c
c-- close file
      close(4)
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
      end module inputstrmod
