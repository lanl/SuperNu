      module groupmod
      implicit none
c
c-- wavelength grid (gridmod has a copy as well)
      integer,target :: grp_ng=0
      real*8,allocatable :: grp_wl(:) !(grp_ng) wavelength grid
      real*8,allocatable :: grp_wlinv(:) !(grp_ng) wavelength grid
c
      type grp_t_cache
       integer :: ic=0
       integer :: nlump !number of groups in the lump
       real*8 :: tempinv,capgreyinv
       real*8 :: speclump,emitlump,caplump
       real*8,pointer :: specarr(:) !(grp_ng)
       integer,pointer :: glumps(:) !(grp_ng)
       logical,pointer :: llumps(:) !(grp_ng) this group belongs to the lump
       integer :: istat
      end type grp_t_cache
c
      save
c
      contains
c
c
      subroutine groupmod_init(ng,wldex,wlmin,wlmax)
c     -------------------------------------------!{{{
      implicit none
      integer,intent(in) :: ng,wldex
      real*8,intent(in) :: wlmin,wlmax
************************************************************************
* setup wavelength grid
************************************************************************
      integer :: ig
c
c-- read wavelength grid from file
      if(ng==0) then
       call read_wlgrid(grp_ng)
      else
       grp_ng = ng
       allocate(grp_wl(grp_ng+1))
       forall(ig=1:grp_ng+1) grp_wl(ig) =
     &   wlmin*(wlmax/dble(wlmin))**((ig-1d0)/grp_ng)
      endif
c
      allocate(grp_wlinv(grp_ng+1))
      grp_wlinv = 1d0/grp_wl
c
      contains
c
      subroutine read_wlgrid(ng)
c     --------------------------!{{{
      use gasmod
      use inputparmod
      implicit none
      integer,intent(out) :: ng
************************************************************************
* read wavelength grid from file
*
* EXAMPLE FILE LAYOUT:
* --------------------
* #wavelength bin boundaries. units: [cm]
* #ncell ngroupmax
* #icell ngroup wlbound
* 10 5
*  1 1 1e-5 32e-5
*  2 5 1e-5 2e-5 4e-5 8e-5 16e-5 32e-5
* 
************************************************************************
      real*8 :: help
      real*8,allocatable :: wlstore(:)
      integer :: ngm,nrm,irr,l,ig
c
      open(4,file='input.wlgrid',status='old')
c
c-- strip header
      do l=1,3
       read(4,*)
      enddo
c-- read dimensions
      read(4,*) nrm,ngm !not used
c
      do l=1,wldex-1
       read(4,*)
      enddo
      read(4,*) irr, ng
c
      allocate(grp_wl(ng+1))
      allocate(wlstore(ng+3))
      rewind(4)
      do l=1,wldex+3
       read(4,*)
      enddo
      read(4,*) wlstore(:)
      close(4)
      grp_wl = wlstore(3:)
      deallocate(wlstore)
c
c-- check monotonicity
      help = 0d0
      do ig=1,ng+1
       if(grp_wl(ig)<=help) stop 'read_wlgrid: wlgrid not increasing'
       help = grp_wl(ig)
      enddo
c
      write(6,*)
      write(6,*) 'read_wlgrid: ng=',ng
!     write(6,*) 'wavelength grid in [cm]'
!     write(6,*) grp_wl
c!}}}
      end subroutine read_wlgrid
c!}}}
      end subroutine groupmod_init
c
c
      elemental function specint0(tempinv,ig)
c     ---------------------------------------!{{{
      use physconstmod
      implicit none
      real*8 :: specint0
      real*8,intent(in) :: tempinv
      integer,intent(in) :: ig
************************************************************************
* Calculate x**3/(exp(x) - 1), where x = h*c/(wl*k*T)
************************************************************************
      real*8,parameter :: ftpi4=15d0/pc_pi**4
      real*8,parameter :: hck=pc_h*pc_c/pc_kb
      real*8 :: x,dx
c
      x = hck*tempinv
      dx = x*abs(grp_wlinv(ig+1) - grp_wlinv(ig))
      x = x*.5d0*(grp_wlinv(ig+1) + grp_wlinv(ig))
c
      specint0 = ftpi4 * dx * x**3/(exp(x) - 1d0)
c!}}}
      end function specint0
c
c
      pure function specintv(tempinv,mode) result(ss)
c     -----------------------------------------------!{{{
      use physconstmod
      implicit none
      real*8 :: ss(grp_ng)
      real*8,intent(in) :: tempinv
      integer,intent(in),optional :: mode
************************************************************************
* Integrate normalized Planck spectrum using Newton-Cotes formulae of
* different degrees.
************************************************************************
      real*8,parameter :: ftpi4=15d0/pc_pi**4
      real*8,parameter :: hck=pc_h*pc_c/pc_kb
      real*8,parameter :: one6th=ftpi4/6d0
      real*8,parameter :: one90th=ftpi4/90d0
      real*8,parameter :: onehalf=ftpi4/2d0
      real*8,parameter :: one=ftpi4
      real*8 :: xarr(grp_ng+1)
      real*8 :: farr(grp_ng+1)
      real*8 :: f2arr(grp_ng)
      integer :: imode
c
c-- default mode is linear
      imode = 1
      if(present(mode)) imode = mode
c
c-- x
      xarr = hck*tempinv*grp_wlinv
c
c-- edge values are always used
c
      select case(imode)
c-- constant
      case(0)
       ss = one*abs(xarr(2:)-xarr(:grp_ng)) *
     &   f(.5d0*(xarr(2:)+xarr(:grp_ng)))
c
c-- linear
      case(1)
       farr = f(xarr) !edge values
       ss = onehalf*abs(xarr(2:)-xarr(:grp_ng)) *
     &   (farr(2:)+farr(:grp_ng))
c
c-- quadratic
      case(2)
       farr = f(xarr) !edge values
       f2arr = f(.5d0*(xarr(2:) + xarr(:grp_ng))) !midpoints
c-- simpson's rule
       ss = one6th*abs(xarr(2:)-xarr(:grp_ng)) *
     &  (farr(2:) + farr(:grp_ng) + 4d0*f2arr)
c
c-- cubic
      case(11)
       farr = f(xarr) !edge values
c-- quarter points
       f2arr = 12d0*f(.5d0*(xarr(2:) + xarr(:grp_ng))) + 32d0*(
     &   f(.25d0*xarr(2:) + .75d0*xarr(:grp_ng)) +
     &   f(.75d0*xarr(2:) + .25d0*xarr(:grp_ng)))
c-- Boole's rule
       ss = one90th*abs(xarr(2:)-xarr(:grp_ng)) *
     &   (7d0*(farr(2:) + farr(:grp_ng)) + f2arr)
c
c-- invalid
      case default
       ss = 0d0
      endselect
c
      contains
c
      elemental real*8 function f(x)
      real*8,intent(in) :: x
      f = x**3/(exp(x) - 1d0)
      end function
c!}}}
      end function specintv
c
c
      pure subroutine specintw(tempinv,ss,mode,mask)
c     ----------------------------------------------------!{{{
      use physconstmod
      implicit none
      real*8,intent(inout) :: ss(grp_ng)
      real*8,intent(in) :: tempinv
      integer,intent(in),optional :: mode
      logical,intent(in),optional :: mask(grp_ng)
************************************************************************
* Integrate normalized Planck spectrum using Newton-Cotes formulae of
* different degrees.
************************************************************************
      real*8,parameter :: ftpi4=15d0/pc_pi**4
      real*8,parameter :: hck=pc_h*pc_c/pc_kb
      real*8,parameter :: one6th=ftpi4/6d0
      real*8,parameter :: one90th=ftpi4/90d0
      real*8,parameter :: onehalf=ftpi4/2d0
      real*8,parameter :: one=ftpi4
      real*8 :: xarr(grp_ng+1)
      real*8 :: farr(grp_ng+1)
      real*8 :: f2arr(grp_ng)
      integer :: i,imode
c
c-- default mode is linear
      imode = 1
      if(present(mode)) imode = mode
c
c-- x
      xarr = hck*tempinv*grp_wlinv
c
c-- edge values are always used
c
      select case(imode)
c-- constant
      case(0)
       if(present(mask)) then
        where(mask) ss = one*abs(xarr(2:)-xarr(:grp_ng)) *
     &    f(.5d0*(xarr(2:)+xarr(:grp_ng)))
       else
        ss = one*abs(xarr(2:)-xarr(:grp_ng)) *
     &    f(.5d0*(xarr(2:)+xarr(:grp_ng)))
       endif
c
c-- linear
      case(1)
       farr = f(xarr) !edge values
       ss = onehalf*abs(xarr(2:)-xarr(:grp_ng)) *
     &   (farr(2:)+farr(:grp_ng))
c
c-- quadratic
      case(2)
       farr = f(xarr) !edge values
       f2arr = f(.5d0*(xarr(2:) + xarr(:grp_ng))) !midpoints
c-- simpson's rule
       ss = one6th*abs(xarr(2:)-xarr(:grp_ng)) *
     &  (farr(2:) + farr(:grp_ng) + 4d0*f2arr)
c
c-- cubic
      case(11)
       farr = f(xarr) !edge values
c-- quarter points
       f2arr = 12d0*f(.5d0*(xarr(2:) + xarr(:grp_ng))) + 32d0*(
     &   f(.25d0*xarr(2:) + .75d0*xarr(:grp_ng)) +
     &   f(.75d0*xarr(2:) + .25d0*xarr(:grp_ng)))
c-- Boole's rule
       ss = one90th*abs(xarr(2:)-xarr(:grp_ng)) *
     &   (7d0*(farr(2:) + farr(:grp_ng)) + f2arr)
c
c-- invalid
      case default
       ss = 0d0
      endselect
c
      contains
c
      elemental real*8 function f(x)
      real*8,intent(in) :: x
      f = x**3/(exp(x) - 1d0)
      end function
c!}}}
      end subroutine specintw
c
c
      pure function specintvp(tempinv,i1,i2,mode) result(ss)
c     -----------------------------------------------!{{{
      use physconstmod
      implicit none
      real*8 :: ss(i2-i1+1)
      real*8,intent(in) :: tempinv
      integer,intent(in) :: i1,i2
      integer,intent(in),optional :: mode
************************************************************************
* Integrate normalized Planck spectrum using Newton-Cotes formulae of
* different degrees.
************************************************************************
      real*8,parameter :: ftpi4=15d0/pc_pi**4
      real*8,parameter :: hck=pc_h*pc_c/pc_kb
      real*8,parameter :: one6th=ftpi4/6d0
      real*8,parameter :: one90th=ftpi4/90d0
      real*8,parameter :: onehalf=ftpi4/2d0
      real*8,parameter :: one=ftpi4
      real*8 :: xarr(i2-i1+2)
      real*8 :: farr(i2-i1+2)
      real*8 :: f2arr(i2-i1+1)
      integer :: imode,n
c
c-- default mode is linear
      imode = 1
      if(present(mode)) imode = mode
c
c-- x
      n = i2 - i1 + 1
      xarr = hck*tempinv*grp_wlinv(i1:i2+1)
c
c-- edge values are always used
c
      select case(imode)
c-- constant
      case(0)
       ss = one*abs(xarr(2:)-xarr(:n)) *
     &   f(.5d0*(xarr(2:)+xarr(:n)))
c
c-- linear
      case(1)
       farr = f(xarr) !edge values
       ss = onehalf*abs(xarr(2:)-xarr(:n)) *
     &   (farr(2:)+farr(:n))
c
c-- quadratic
      case(2)
       farr = f(xarr) !edge values
       f2arr = f(.5d0*(xarr(2:) + xarr(:n))) !midpoints
c-- simpson's rule
       ss = one6th*abs(xarr(2:)-xarr(:n)) *
     &  (farr(2:) + farr(:n) + 4d0*f2arr)
c
c-- cubic
      case(11)
       farr = f(xarr) !edge values
c-- quarter points
       f2arr = 12d0*f(.5d0*(xarr(2:) + xarr(:n))) + 32d0*(
     &   f(.25d0*xarr(2:) + .75d0*xarr(:n)) +
     &   f(.75d0*xarr(2:) + .25d0*xarr(:n)))
c-- Boole's rule
       ss = one90th*abs(xarr(2:)-xarr(:n)) *
     &   (7d0*(farr(2:) + farr(:n)) + f2arr)
c
c-- invalid
      case default
       ss = 0d0
      endselect
c
      contains
c
      elemental real*8 function f(x)
      real*8,intent(in) :: x
      f = x**3/(exp(x) - 1d0)
      end function
c!}}}
      end function specintvp
c
c
      end module groupmod
