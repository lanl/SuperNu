*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module groupmod
      implicit none
c
c-- wavelength grid (gridmod has a copy as well)
      integer,target :: grp_ng,grp_ngs
      real*8 :: grp_wlmin,grp_wlmax
      real*8,allocatable :: grp_wl(:) !(grp_ng) wavelength grid
      real*8,allocatable :: grp_wlinv(:) !(grp_ng) wavelength grid
c
      type grp_t_cache
       integer :: ic=0
       integer :: nlump !number of groups in the lump
       real*8 :: capgreyinv
       real*8 :: speclump,emitlump,caplump,doplump
       real*8,pointer :: specarr(:) !(grp_ng)
       integer :: istat
       integer*2,pointer :: glumps(:) !(grp_ng)
       logical*2,pointer :: llumps(:) !(grp_ng) this group belongs to the lump
      end type grp_t_cache
c
      save
c
      contains
c
c
      subroutine groupmod_init(wldex)
c     -------------------------------!{{{
      implicit none
      integer,intent(in) :: wldex
************************************************************************
* setup wavelength grid
************************************************************************
      integer :: ig
c
c-- verify data range
      if(grp_ng>2**15-1) stop 'groupmod_init: grp_ng > 2^15. '//
     &  'Increase glumps int*2 data range'
c
c-- read wavelength grid from file
      if(grp_ng==0) then
       call read_wlgrid(wldex,grp_ng)
      else
       allocate(grp_wl(grp_ng+1))
       forall(ig=1:grp_ng+1) grp_wl(ig) =
     &   grp_wlmin*(grp_wlmax/dble(grp_wlmin))**((ig-1d0)/grp_ng)
      endif
c
      allocate(grp_wlinv(grp_ng+1))
      grp_wlinv = 1d0/grp_wl
c
      contains
c
      subroutine read_wlgrid(wldex,ng)
c     --------------------------------!{{{
      implicit none
      integer,intent(in) :: wldex
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
      elemental function dopspeccalc(tempinv,ig)
c     ---------------------------------------!{{{
      use physconstmod
      implicit none
      real*8 :: dopspeccalc
      real*8,intent(in) :: tempinv
      integer,intent(in) :: ig
************************************************************************
* Calculate x**3/(exp(x) - 1), where x = h*c/(wl*k*T)
************************************************************************
      real*8,parameter :: ftpi4=15d0/pc_pi**4
      real*8,parameter :: hck=pc_h*pc_c/pc_kb
      real*8 :: x
c
      x = hck*tempinv
      x = x*grp_wlinv(ig+1)
c
      dopspeccalc = ftpi4 * x**4/(exp(x) - 1d0)
c!}}}
      end function dopspeccalc
c
c
      pure subroutine specintv(tempinv,n,ss,offset,mode,mask,maskval)
c     -----------------------------------------------!{{{
      use physconstmod
      implicit none
      real*8,intent(in) :: tempinv
      integer,intent(in) :: n
      real*8,intent(out) :: ss(n)
      integer,intent(in),optional :: offset
      integer,intent(in),optional :: mode
      logical*2,intent(in),optional :: mask(n)
      logical,intent(in),optional :: maskval
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
      real*8 :: xarr(n+1)
      real*8 :: farr(n+1)
      real*8 :: f2arr(n)
      integer :: imode,ioff
c
c-- defaults
      ioff = 1
      if(present(offset)) ioff = offset
      imode = 0
      if(present(mode)) imode = mode
c
c-- sanity check
      if(ioff+n>grp_ng+1 .or. present(mask).neqv.present(maskval)) then
       ss = -1d0
       return
      endif
c
c-- x
      xarr = hck*tempinv*grp_wlinv(ioff:ioff+n)
c
      select case(imode)
c-- constant
      case(0)
       if(present(mask)) then
        where(mask.eqv.maskval)
         ss = one*abs(xarr(2:)-xarr(:n)) *
     &    f(.5d0*(xarr(2:)+xarr(:n)))
        endwhere
       else
         ss = one*abs(xarr(2:)-xarr(:n)) *
     &     f(.5d0*(xarr(2:)+xarr(:n)))
       endif
c
c-- linear
      case(1)
       if(present(mask)) then
        where(mask.eqv.maskval)
         farr = f(xarr) !edge values
         ss = onehalf*abs(xarr(2:)-xarr(:n)) *
     &     (farr(2:)+farr(:n))
        endwhere
       else
         farr = f(xarr) !edge values
         ss = onehalf*abs(xarr(2:)-xarr(:n)) *
     &     (farr(2:)+farr(:n))
       endif
c
c-- quadratic
      case(2)
       if(present(mask)) then
        where(mask.eqv.maskval)
         farr = f(xarr) !edge values
         f2arr = f(.5d0*(xarr(2:) + xarr(:n))) !midpoints
c-- simpson's rule
         ss = one6th*abs(xarr(2:)-xarr(:n)) *
     &    (farr(2:) + farr(:n) + 4d0*f2arr)
        endwhere
       else
         farr = f(xarr) !edge values
         f2arr = f(.5d0*(xarr(2:) + xarr(:n))) !midpoints
c-- simpson's rule
         ss = one6th*abs(xarr(2:)-xarr(:n)) *
     &    (farr(2:) + farr(:n) + 4d0*f2arr)
       endif
c
c-- cubic
      case(11)
       if(present(mask)) then
        where(mask.eqv.maskval)
         farr = f(xarr) !edge values
c-- quarter points
         f2arr = 12d0*f(.5d0*(xarr(2:) + xarr(:n))) + 32d0*(
     &     f(.25d0*xarr(2:) + .75d0*xarr(:n)) +
     &     f(.75d0*xarr(2:) + .25d0*xarr(:n)))
c-- Boole's rule
         ss = one90th*abs(xarr(2:)-xarr(:n)) *
     &     (7d0*(farr(2:) + farr(:n)) + f2arr)
        endwhere
       else
         farr = f(xarr) !edge values
c-- quarter points
         f2arr = 12d0*f(.5d0*(xarr(2:) + xarr(:n))) + 32d0*(
     &     f(.25d0*xarr(2:) + .75d0*xarr(:n)) +
     &     f(.75d0*xarr(2:) + .25d0*xarr(:n)))
c-- Boole's rule
         ss = one90th*abs(xarr(2:)-xarr(:n)) *
     &     (7d0*(farr(2:) + farr(:n)) + f2arr)
       endif
c
c-- invalid
      case default
       ss = -1d0
      endselect
c
      contains
c
      elemental real*8 function f(x)
      real*8,intent(in) :: x
      f = x**3/(exp(x) - 1d0)
      end function
c!}}}
      end subroutine specintv
c
c
      end module groupmod
c vim: fdm=marker
