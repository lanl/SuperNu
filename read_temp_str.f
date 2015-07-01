      subroutine read_temp_str(igeom,ndim)
c     --------------------------------------------------------
      use inputstrmod
      implicit none
      integer,intent(in) :: igeom
      integer,intent(in) :: ndim(3)
************************************************************************
* Read the input gas temperature structure file if it exists.
* Like read_inputstr, but for 1 column (header info is sanity checked).
************************************************************************
      integer :: igeom_r,nx_r,ny_r,nz_r
      integer :: nx,ny,nz,l, ierr
      character(2) :: dmy
      real*8,allocatable :: raw(:)
c
c-- check for input.str_temp
      inquire(file='input.str_temp',exist=str_ltemp)
c-- quick exit
      if(.not.str_ltemp) return
c
c-- copy
      nx=ndim(1)
      ny=ndim(2)
      nz=ndim(3)
c
c-- open file
      open(4,file='input.str_temp',status='old')
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy, igeom_r,nx_r,ny_r,nz_r
c-- verify geometry and dimension
      if(igeom_r/=igeom) stop 'read_temp_str: incompatible igeom'
      if(nx_r/=nx) stop 'read_temp_str: incompatible nx dimension'
      if(ny_r/=ny) stop 'read_temp_str: incompatible ny dimension'
      if(nz_r/=nz) stop 'read_temp_str: incompatible nz dimension'
c
c-- allocate temperature structure
      allocate(str_temp(nx,ny,nz))
c
c-- read body
      read(4,*,iostat=ierr) str_temp
      if(ierr/=0) stop 'read_temp_str: input.str_temp fmt err: body'
      read(4,*,iostat=ierr) dmy
      if(ierr/=-1) stop 'read_temp_str: input.str_temp body too long'
c
c-- validity check
      if(any(temp/=temp)) stop 'read_temp_str: nan in input'
c
c-- close file
      close(4)
c
      end subroutine read_temp_str
