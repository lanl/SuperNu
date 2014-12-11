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
      integer,private :: nx,ny,nz
      integer,private :: igeom
c
      save
      public
c
      contains
c
c
c
      subroutine inputstr_dealloc
c     ---------------------------!{{{
      implicit none
      str_nabund=0
      if(allocated(str_iabund)) deallocate(str_iabund)
      deallocate(str_xleft,str_yleft,str_zleft)
      deallocate(str_mass,str_massfr)
      deallocate(str_massdd,str_massfrdd)!}}}
      end subroutine inputstr_dealloc
c
c
      subroutine read_inputstr(igeomin,ndim,lvoidcorners)
c     ---------------------------------------------------!{{{
      use physconstmod
      use gasmod, only:gas_ini56,gas_ico56
      use miscmod
      implicit none
      integer,intent(in) :: igeomin
      integer,intent(in) :: ndim(3)
      logical,intent(in) :: lvoidcorners
************************************************************************
* Read the input structure file
************************************************************************
      integer :: i,j,k,ierr,nx_r,ny_r,nz_r,ini56,nvar,ncol,imass,nvoid
      character(2) :: dmy
      character(8),allocatable :: labl(:)
      real*8,allocatable :: raw(:,:)
      real*8 :: mni56
      real*8 :: r,rs
c-- statement functions
      real*8 :: x,y,z
      x(i) = .5d0*(str_xleft(i+1) + str_xleft(i))
      y(i) = .5d0*(str_yleft(i+1) + str_yleft(i))
      z(i) = .5d0*(str_zleft(i+1) + str_zleft(i))
c
c-- copy
      igeom = igeomin
      nx = ndim(1)
      ny = ndim(2)
      nz = ndim(3)
c
c-- open file
      open(4,file=fname,status='old',iostat=ierr)
      if(ierr/=0) stop 'read_inputstr: file missing: input.str'
c
c-- read dimensions
      read(4,*)
      read(4,*,iostat=ierr) dmy, nx_r,ny_r,nz_r,ncol,str_nabund
      if(ierr/=0) stop 'read_inputstr: input.str fmt err: dimensions'
c-- verify dimension
      if(nx_r/=nx) stop 'read_inputstr: incompatible nx dimension'
      if(ny_r/=ny) stop 'read_inputstr: incompatible ny dimension'
      if(nz_r/=nz) stop 'read_inputstr: incompatible nz dimension'
c
c-- allocate label arrays
      nvar = ncol - str_nabund
      allocate(str_abundlabl(str_nabund))
      allocate(labl(nvar))
c
c-- read labels
      read(4,*,iostat=ierr) dmy, labl, str_abundlabl
      if(ierr/=0) stop 'read_inputstr: input.str fmt err: col labels'
c
c-- allocate data arrays
      allocate(str_xleft(nx+1))
      allocate(str_yleft(ny+1))
      allocate(str_zleft(nz+1))
      allocate(str_mass(nx,ny,nz))
      allocate(str_massfr(str_nabund,nx,ny,nz))
      allocate(raw(ncol,nx*ny*nz))
c
c-- read body
      read(4,*,iostat=ierr) raw
      if(ierr/=0) stop 'read_inputstr: input.str format err: body'
      read(4,*,iostat=ierr) dmy
      if(ierr/=-1) stop 'read_inputstr: input.str body too long'
c
c-- validity check
      if(any(raw/=raw)) stop 'read_inputstr: nan in input'
c
c-- transer data to final arrays
c-- first dim
      if(igeom==1) then
       str_xleft(1) = 0d0
       str_xleft(2:) = raw(1,:nx)
      else
       str_xleft(1) = raw(1,1)
       str_xleft(2:) = raw(2,:nx)
      endif
c-- second dim
      if(igeom>1) then
       str_yleft(1) = raw(3,1)
       do i=1,ny
        str_yleft(i+1) = raw(4,nx*(i-1)+1)
       enddo
      else
       str_yleft = 0d0
      endif
c-- third dim
      if(igeom>2) then
       str_zleft(1) = raw(5,1)
       do i=1,nz
        str_zleft(i+1) = raw(6,nx*ny*(i-1)+1)
       enddo
      else
       str_zleft = 0d0
      endif
c-- var pointers
      imass = 0
      do i=1,nvar
       if(lcase(trim(labl(i)))=='mass') imass = i
      enddo
      if(imass==0) stop 'read_inputstr: mass label not found'
c-- vars
      str_mass(:,:,:) = reshape(raw(imass,:),[nx,ny,nz])
c-- abundances
      str_massfr(:,:,:,:) =reshape(raw(nvar+1:,:),[str_nabund,nx,ny,nz])
c
c-- close file
      close(4)
      deallocate(raw)
c
c-- convert abundlabl to element codes
      call elnam2elcode(ini56)
c
c
c-- Zero out the cell mass in the corners of the domain
c======================================================
c-- void cells
      nvoid = 0
      if(lvoidcorners .and. igeom>1) then
c-- sphere radius
       rs = min(str_xleft(nx+1),str_yleft(ny+1))
       if(igeom>2) rs = min(rs,str_zleft(nz+1))
c
       do k=1,nz
       do j=1,ny
       do i=1,nx
        r = sqrt(x(i)**2 + y(j)**2 + z(k)**2)
        if(r>rs) then
         nvoid = nvoid + 1
         str_mass(i,j,k) = 0d0
        endif
       enddo
       enddo
       enddo
       write(6,*) '# corner cells voided:', nvoid
      endif !lvoidcorners
c
c-- output
c
c-- ni56 mass
      if(ini56>0) then
       mni56 = sum(str_massfr(ini56,:,:,:)*str_mass)
      else
       mni56 = 0d0
      endif
!c-- kinetic energy
!      ekin = 
c
c-- output
      write(6,*)
      write(6,*) 'input structure:'
      write(6,*) '===================='
      write(6,*) 'igeom :',igeom
      write(6,*) 'ndim  :',nx,ny,nz, nx*nx*nz
      write(6,*) 'mass  :', sum(str_mass)/pc_msun, 'Msun'
      write(6,*) 'm_ni56:', mni56/pc_msun, 'Msun'
!     write(6,*) 'e_kin :', ekin, 'erg'
c!}}}
      end subroutine read_inputstr
c
c
c
      subroutine generate_inputstr(igeomin)
c     ---------------------------------------------!{{{
      implicit none
      integer,intent(in) :: igeomin
************************************************************************
* wrapper around routines for different geometries
************************************************************************
      igeom = igeomin
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
      subroutine generate_inputstr1
      use inputparmod!{{{
      implicit none
************************************************************************
* generate stratification from input.par variables
* if in_noreadstruct==.true.
************************************************************************
      real*8,allocatable :: rout(:) !(nx+1)
      integer :: i
      real*8 :: help, dx
c
c-- 1D size
      nx = in_ndim(1)
      ny = 1
      nz = 1
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
      use physconstmod
      implicit none
************************************************************************
* Read the input structure file
************************************************************************
      real*8,allocatable :: xout(:) !(nx+1)
      real*8,allocatable :: yout(:) !(ny+1)
      integer :: i,j
      real*8 :: helpx,helpy, dx,dy,rmax
c
c-- 2D size
      nx = in_ndim(1)
      ny = in_ndim(2)
      nz = 1
c
c-- verifications (input.par)
      if(in_velout<=0d0.and.in_isvelocity)
     &  stop 'generate_inputstr2: invalid in_velout'
      if(in_lx<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr2: invalid in_lx'
      if(in_ly<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr2: invalid in_ly'
      if(in_lz<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr2: invalid in_lz'
      if(in_totmass<0d0)
     &  stop 'generate_inputstr2: invalid in_totmass'
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
      use physconstmod
      implicit none
************************************************************************
* Read the input structure file
************************************************************************
      real*8,allocatable :: xout(:) !(nx+1)
      real*8,allocatable :: yout(:) !(ny+1)
      real*8,allocatable :: zout(:) !(nz+1)
      integer :: i,j,k
      real*8 :: helpx,helpy,helpz, dx,dy,dz,rmax
c
c-- 3D size
      nx = in_ndim(1)
      ny = in_ndim(2)
      nz = in_ndim(3)
c
c-- verifications (input.par)
      if(in_velout<=0d0.and.in_isvelocity)
     &  stop 'generate_inputstr3: invalid in_velout'
      if(in_lx<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr3: invalid in_lx'
      if(in_ly<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr3: invalid in_ly'
      if(in_lz<=0.and..not.in_isvelocity)
     &  stop 'generate_inputstr3: invalid in_lz'
      if(in_totmass<0d0)
     &  stop 'generate_inputstr3: invalid in_totmass'
c
c-- allocate arrays
      allocate(xout(nx+1))
      allocate(yout(ny+1))
      allocate(zout(nz+1))
      allocate(str_xleft(nx+1))
      allocate(str_yleft(ny+1))
      allocate(str_zleft(nz+1))
      allocate(str_mass(nx,ny,nz))
c
c-- create unit-length x array
      dx = 1d0/real(nx)
      forall(i=1:nx+1) xout(i) = -0.5d0+(i-1)*dx
c
c-- create unit-length y array
      dy = 1d0/real(ny)
      forall(j=1:ny+1) yout(j) = -0.5d0+(j-1)*dy
c
c-- create unit-length z array
      dz = 1d0/real(nz)
      forall(k=1:nz+1) zout(k) = -0.5d0+(k-1)*dz
c
c-- dimensional scaling
      if(in_isvelocity) then
       helpx = 2*in_velout
       helpy = 2*in_velout
       helpz = 2*in_velout
      else
       helpx = in_lx
       helpy = in_ly
       helpz = in_lz
      endif
      str_xleft = helpx*xout
      str_yleft = helpy*yout
      str_zleft = helpz*zout
c
c-- mass
      str_mass = 0d0
      if(in_dentype=='unif') then
c-- uniform density sphere
       rmax = min(helpx/2d0,helpy/2d0,helpz/2d0)
       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
           helpx = 0.5*(str_xleft(i+1)+str_xleft(i))
           helpy = 0.5*(str_yleft(j+1)+str_yleft(j))
           helpz = 0.5*(str_zleft(k+1)+str_zleft(k))
           if(helpx**2+helpy**2+helpz**2<rmax**2) then
            str_mass(i,j,k)=(str_xleft(i+1)-str_xleft(i)) *
     &        (str_yleft(j+1)-str_yleft(j)) *
     &        (str_zleft(k+1)-str_zleft(k)) *
     &        in_totmass/(pc_pi43*rmax**3)
           endif
         enddo
        enddo
       enddo
      elseif(in_dentype=='mass') then
c-- spherical 1/r^2 mass shells
       rmax = min(helpz/2d0,helpy/2d0,helpz/2d0)
       do k = 1,nz
        do j = 1,ny
         do i = 1,nx
           helpx = 0.5*(str_xleft(i+1)+str_xleft(i))
           helpy = 0.5*(str_yleft(j+1)+str_yleft(j))
           helpz = 0.5*(str_zleft(k+1)+str_zleft(k))
           if(helpx**2+helpy**2+helpz**2<rmax**2) then
            str_mass(i,j,k)=(str_xleft(i+1)-str_xleft(i)) *
     &        (str_yleft(j+1)-str_yleft(j)) *
     &        (str_zleft(k+1)-str_zleft(k)) *
     &        in_totmass/(pc_pi4*rmax*(helpy**2+helpx**2 +
     &        helpz**2))
           endif
         enddo
        enddo
       enddo
      elseif(in_dentype=='ufil'.or.in_dentype=='mfil') then
c-- uniform density for all box
       str_mass = in_totmass/(nx*ny*nz)
      else
       stop 'generate_inputstr3: invalid in_dentype'
      endif
c-- adjusting mass to correct total
      str_mass = str_mass*in_totmass/sum(str_mass)
c-- deallocating helper arrays
      deallocate(xout,yout,zout)
c!}}}
      end subroutine generate_inputstr3
c
c
c
c
      subroutine elnam2elcode(ini56)
c     ------------------------------!{{{
      use miscmod, only:lcase
      use gasmod, only:gas_ini56,gas_ico56
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
