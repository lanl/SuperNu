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
      real*8,allocatable :: str_xleft(:) !(nx+1)
      real*8,allocatable :: str_yleft(:) !(ny+1)
      real*8,allocatable :: str_zleft(:) !(nz+1)
      real*8,allocatable :: str_mass(:,:,:) !(nx,ny,nz)
      real*8,allocatable :: str_massfr(:,:,:,:) !(nabund,nx,ny,nz)
c
c-- domain decomposition
      real*8,allocatable :: str_massdd(:) !(gas_ncell)
      real*8,allocatable :: str_massfrdd(:,:) !(nabund,gas_ncell)
c
      character(8),allocatable,private :: str_abundlabl(:) !(nabund)
c
      save
      public
c
      contains
c
c
      subroutine read_inputstr(igeom,ndim)
c     -------------------------------------!{{{
      implicit none
      integer,intent(in) :: igeom
      integer,intent(in) :: ndim(3)
************************************************************************
* wrapper around routines for different geometries
************************************************************************
      select case(igeom)
      case(1)
       call read_inputstr1(ndim)
      case(2)
       call read_inputstr2(ndim)
      case(3)
       call read_inputstr3(ndim)
      case default
       stop 'read_inputstr: invalid igeom'
      endselect
c
c-- allocate remaining arrays
      if(.not.allocated(str_yleft)) allocate(str_yleft(2))
      if(.not.allocated(str_zleft)) allocate(str_zleft(2))
c!}}}
      end subroutine read_inputstr
c
c
      subroutine generate_inputstr(igeom)
c     ---------------------------------------------!{{{
      implicit none
      integer,intent(in) :: igeom
************************************************************************
* wrapper around routines for different geometries
************************************************************************
      select case(igeom)
      case(1)
       call generate_inputstr1
      case(2)
       call generate_inputstr2
      case(3)
       call generate_inputstr3
      case default
       stop 'generate_inputstr: invalid igeom'
      endselect
c
c-- allocate remaining arrays
      if(.not.allocated(str_yleft)) allocate(str_yleft(2))
      if(.not.allocated(str_zleft)) allocate(str_zleft(2))
c!}}}
      end subroutine generate_inputstr
c
c
c
c
      subroutine read_inputstr1(ndim)
c     -----------------------------------!{{{
      use physconstmod
      use gasgridmod, only:gas_ini56,gas_ico56
      use miscmod
      implicit none
      integer,intent(in) :: ndim(3)
************************************************************************
* Read the input structure file
************************************************************************
      integer :: ierr,nx,nx_r,ini56
      character(2) :: dmy
      character(8) :: labl(2)
      real*8,allocatable :: raw(:,:)
      real*8 :: help
c
c-- copy
      nx = ndim(1)
c
c-- open file
      open(4,file=fname,status='old',iostat=ierr)
      if(ierr/=0) stop 'read_inputstr1: file missing: input.str'
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy,nx_r,str_nabund
      if(ierr/=0) stop 'read_inputstr1: input.str fmt err: dimensions'
c-- verify dimension
      if(nx_r/=nx) stop 'read_inputstr1: incompatible nx dimension'
c
c-- allocate arrays
      allocate(str_xleft(nx+1))
      allocate(str_mass(nx,1,1))
      allocate(str_massfr(str_nabund,nx,1,1))
      allocate(str_abundlabl(str_nabund))
      allocate(raw(str_nabund+2,nx))
c
c-- read labels
      read(4,*,iostat=ierr) dmy, labl, str_abundlabl
      if(ierr/=0) stop 'read_inputstr1: input.str fmt err: col labels'
c
c-- read body
      read(4,*,iostat=ierr) raw
      if(ierr/=0) stop 'read_inputstr1: input.str format err: body'
c
c-- transer data to final arrays
      str_xleft(1) = 0d0
      str_xleft(2:) = raw(1,:)
      str_mass(:,1,1) = raw(2,:)
      str_massfr(:,:,1,1) = raw(3:,:)
c
c-- close file
      close(4)
      deallocate(raw)
c
c-- convert abundlabl to element codes
      call elnam2elcode(ini56)
c
c-- output
      write(6,*)
      write(6,*) 'input structure 1D:'
      write(6,*) '===================='
      write(6,*) 'mass  :', sum(str_mass)/pc_msun, 'Msun'
c-- ni56 mass
      if(ini56>0) then
       help = sum(str_massfr(ini56,:,:,:)*str_mass)
      else
       help = 0d0
      endif
      write(6,*) 'm_ni56:', help/pc_msun, 'Msun'
c!}}}
      end subroutine read_inputstr1
c
c
      subroutine read_inputstr2(ndim)
c     -----------------------------------!{{{
      use physconstmod
      use gasgridmod, only:gas_ini56,gas_ico56
      use miscmod
      implicit none
      integer,intent(in) :: ndim(3)
************************************************************************
* Read the input structure file
************************************************************************
      integer :: ierr,nx,ny,nx_r,ny_r,nz_r,ini56,j
      character(2) :: dmy
      character(8) :: labl(5)
      real*8,allocatable :: raw(:,:)
      real*8 :: help
c
c-- copy
      nx = ndim(1)
      ny = ndim(2)
c
c-- open file
      open(4,file=fname,status='old',iostat=ierr)
      if(ierr/=0) stop 'read_inputstr2: file missing: input.str'
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy, nx_r,ny_r,nz_r,str_nabund
      if(ierr/=0) stop 'read_inputstr2: input.str fmt err: dimensions'
c-- verify dimension
      if(nx_r/=nx) stop 'read_inputstr2: incompatible nx dimension'
      if(ny_r/=ny) stop 'read_inputstr2: incompatible ny dimension'
c
c-- allocate arrays
      allocate(str_xleft(nx+1))
      allocate(str_yleft(ny+1))
      allocate(str_mass(nx,ny,1))
      allocate(str_massfr(str_nabund,nx,ny,1))
      allocate(str_abundlabl(str_nabund))
      allocate(raw(5+str_nabund,nx*ny))
c
c-- read labels
      read(4,*,iostat=ierr) dmy, labl, str_abundlabl
      if(ierr/=0) stop 'read_inputstr2: input.str fmt err: col labels'
c
c-- read body
      read(4,*,iostat=ierr) raw
      if(ierr/=0) stop 'read_inputstr2: input.str format err: body'
c
c-- transer data to final arrays
      str_xleft(1) = raw(1,1)
      str_xleft(2:) = raw(2,:nx)
      str_yleft(1) = raw(3,1)
      do j=1,ny
       str_yleft(j+1) = raw(4,nx*(j-1)+1)
      enddo
      str_mass(:,:,1) = reshape(raw(5,:),[nx,ny])
      str_massfr(:,:,:,1) = reshape(raw(6:,:),[str_nabund,nx,ny])
c
c-- close file
      close(4)
      deallocate(raw)
c
c-- convert abundlabl to element codes
      call elnam2elcode(ini56)
c
c-- output
      write(6,*)
      write(6,*) 'input structure 2D:'
      write(6,*) '===================='
      write(6,*) 'mass  :', sum(str_mass)/pc_msun, 'Msun'
c-- ni56 mass
      if(ini56>0) then
       help = sum(str_massfr(ini56,:,:,:)*str_mass)
      else
       help = 0d0
      endif
      write(6,*) 'm_ni56:', help/pc_msun, 'Msun'
c
      deallocate(str_abundlabl)
c!}}}
      end subroutine read_inputstr2
c
c
      subroutine read_inputstr3(ndim)
c     -----------------------------------!{{{
      use physconstmod
      use gasgridmod, only:gas_ini56,gas_ico56
      use miscmod
      implicit none
      integer,intent(in) :: ndim(3)
************************************************************************
* Read the input structure file
************************************************************************
      integer :: ini56
      real*8 :: help
c
c-- dummy
      stop 'read_inputstr3 is a stub'
      write(6,*) ndim !use ndim so compiler doesn't complain
c
c-- convert abundlabl to element codes
      call elnam2elcode(ini56)
c
c-- output
      write(6,*)
      write(6,*) 'input structure 3D:'
      write(6,*) '===================='
      write(6,*) 'mass  :', sum(str_mass)/pc_msun, 'Msun'
c-- ni56 mass
      if(ini56>0) then
       help = sum(str_massfr(ini56,:,:,:)*str_mass)
      else
       help = 0d0
      endif
      write(6,*) 'm_ni56:', help/pc_msun, 'Msun'
c!}}}
      end subroutine read_inputstr3
c
c
c
c
      subroutine generate_inputstr1
      use inputparmod!{{{
      implicit none
************************************************************************
* generate stratification from input.par variables
* if in_noreadstruct==.true.
************************************************************************
      real*8,allocatable :: rout(:) !(nx+1)
      integer :: i,nx
      real*8 :: help, dx
c
c-- 1D size
      nx = in_ndim(1)
c
c-- verifications (input.par)
      if(in_velout<=0d0.and.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_velout'
      if(in_lx<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_lx'
      if(in_ly<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_ly'
      if(in_lz<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_lz'
      if(in_totmass<0d0)
     &  stop 'generate_inputstr1: invalid in_totmass'
c
c-- allocate arrays
      allocate(rout(nx+1))
      allocate(str_xleft(nx+1))
      allocate(str_mass(nx,1,1))
c
c-- create unit sphere radii rout
      dx = 1d0/real(nx)
      forall(i=1:nx+1) rout(i) = (i-1)*dx
c
c-- outer shells
      if(in_isvelocity) then
       help = in_velout
      else
       help = in_lx
      endif
      str_xleft = help*rout
c
c-- mass
      if(in_dentype=='unif') then
       str_mass(:,1,1) = in_totmass*(rout(2:)**3 - rout(:nx)**3)
       str_mass(:,1,1) = str_mass(:,1,1)/(1d0 - rout(1)**3)
      elseif(in_dentype=='mass') then
       forall(i=1:nx)str_mass(i,1,1) = in_totmass/real(nx)
      else
       stop 'generate_inputstr1: invalid in_dentype'
      endif
      deallocate(rout)
c!}}}
      end subroutine generate_inputstr1
c
c
      subroutine generate_inputstr2
      use inputparmod!{{{
      implicit none
************************************************************************
* Read the input structure file
************************************************************************
      real*8,allocatable :: xout(:) !(nx+1)
      real*8,allocatable :: yout(:) !(ny+1)
      integer :: i,j,nx,ny
      real*8 :: helpx,helpy, dx,dy,rmax
c
c-- 2D size
      nx = in_ndim(1)
      ny = in_ndim(2)
c
c-- verifications (input.par)
      if(in_velout<=0d0.and.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_velout'
      if(in_lx<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_lx'
      if(in_ly<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_ly'
      if(in_lz<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr1: invalid in_lz'
      if(in_totmass<0d0)
     &  stop 'generate_inputstr1: invalid in_totmass'
c
c-- allocate arrays
      allocate(xout(nx+1))
      allocate(yout(ny+1))
      allocate(str_xleft(nx+1))
      allocate(str_yleft(ny+1))
      allocate(str_mass(nx,ny,1))
c
c-- create unit cylinder radii xout
      dx = 1d0/real(nx)
      forall(i=1:nx+1) xout(i) = (i-1)*dx
c
c-- create unit cylinder heights yout
      dy = 1d0/real(ny)
      forall(j=1:ny+1) yout(j) = -0.5d0+(j-1)*dy
c
c-- dimensional scaling
      if(in_isvelocity) then
       helpx = in_velout
       helpy = 2*in_velout
      else
       helpx = in_lx
       helpy = in_ly
      endif
      str_xleft = helpx*xout
      str_yleft = helpy*yout
c
c-- mass
      str_mass = 0d0
      if(in_dentype=='unif') then
c-- uniform density sphere
       rmax = min(helpy/2d0,helpx)
       do j = 1,ny
        do i = 1,nx
          helpx = 0.5*(str_xleft(i+1)+str_xleft(i))
          helpy = 0.5*(str_yleft(j+1)+str_yleft(j))
          if(helpy**2+helpx**2<rmax**2) then
           str_mass(i,j,1)=(str_yleft(j+1)-str_yleft(j)) *
     &       (str_xleft(i+1)**2-str_xleft(i)**2) *
     &       in_totmass/(pc_pi43*rmax**3)
          endif
        enddo
       enddo
      elseif(in_dentype=='mass') then
c-- spherical 1/r^2 mass shells
       rmax = min(helpy/2d0,helpx)
       do j = 1,ny
        do i = 1,nx
          helpx = 0.5*(str_xleft(i+1)+str_xleft(i))
          helpy = 0.5*(str_yleft(j+1)+str_yleft(j))
          if(helpy**2+helpx**2<rmax**2) then
           str_mass(i,j,1)=(str_yleft(j+1)-str_yleft(j)) *
     &       (str_xleft(i+1)**2-str_xleft(i)**2) *
     &       in_totmass/(pc_pi4*rmax*(helpy**2+helpx**2))
          endif
        enddo
       enddo
      elseif(in_dentype=='ufil') then
c-- uniform density for all cylinder
       forall(j=1:ny) str_mass(:,j,1) = in_totmass *
     &        (xout(2:)**2 - xout(:nx)**2)*dy
      elseif(in_dentype=='mfil') then
c-- equal mass per cell for all cylinder
       forall(j=1:ny)str_mass(:,j,1) = in_totmass*dx*dy
      else
       stop 'generate_inputstr2: invalid in_dentype'
      endif
c-- adjusting mass to correct total
      str_mass = str_mass*in_totmass/sum(str_mass)
c-- deallocating helper arrays
      deallocate(xout,yout)
c!}}}
      end subroutine generate_inputstr2
c
c
      subroutine generate_inputstr3
      use inputparmod!{{{
      implicit none
************************************************************************
* Read the input structure file
************************************************************************
      stop 'generate_inputstr3 is a stub'
c!}}}
      end subroutine generate_inputstr3
c
c
c
c
      subroutine elnam2elcode(ini56)
c     ------------------------------!{{{
      use miscmod, only:lcase
      use gasgridmod, only:gas_ini56,gas_ico56
      use elemdatamod
      implicit none
      integer,intent(out) :: ini56
************************************************************************
* convert the abundlabl labels to element codes (atomic z number), which
* also serve as indexing pointers to the mass0fr array.
************************************************************************
      integer :: l,j,iabund
      character(4) :: elname
c
c-- allocate element code array (pointer to mass0fr)
      allocate(str_iabund(str_nabund))
c
c-- default
      ini56 = 0
c
c-- determine atomic z number
      do l=1,str_nabund
       iabund = 0
       elname = lcase(trim(str_abundlabl(l)))
       if(elname=='ni56') then
c-- special case
        iabund = gas_ini56
        ini56 = l
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
c-- store element code (pointer to mass0fr)
       str_iabund(l) = iabund
      enddo!}}}
      end subroutine elnam2elcode
c
      end module inputstrmod
