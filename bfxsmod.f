      module bfxsmod
c     --------------
      implicit none
c
      integer,private :: l
      integer,private :: ll(7),ninn(30),ntot(30)
      real,private :: ph1(6,7,30,30),ph2(7,30,30)
c
      data (ll(l),l=1,7) /0,0,1,0,1,2,0/
c
      data (ninn(l),l=1,30) /0,0,1,1,1,1,1,1,1,1,3,3,
     & 3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5/
c
      data (ntot(l),l=1,30) /1,1,2,2,3,3,3,3,3,3,4,4,
     & 5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/
c
c-- block data
!     include 'bf_verner.blk'
c
      save
c
      contains
c
c
c
      subroutine bfxs_read_data
c------------------------------
      implicit none
************************************************************************
* Read Verner 1995 (ph1) and 1996 (ph2) data.
************************************************************************
      character(14),parameter :: fname='data.bf_verner'
      integer :: ne,nz,ns,istat
      real :: ph1d(9,1699),ph2d(9,465)
c-- read
      open(4,file=fname,action='read',status='old',
     &  iostat=istat)
      if(istat/=0) stop 'bfxs_read_data: cannot read bf_verner.dat'
      read(4,*) !comments
      read(4,*) !comments
      read(4,*,err=99) ph1d
      read(4,*) !comments
      read(4,*) !comments
      read(4,*,err=99) ph2d
      close(4)
c-- parse
      ph1 = 0.
      do l=1,1699
       ns = nint(ph1d(1,l))
       ne = nint(ph1d(2,l))
       nz = nint(ph1d(3,l))
       ph1(:,ns,nz,ne) = ph1d(4:,l)
      enddo
      ph2 = 0.
      do l=1,465
       ne = nint(ph2d(1,l))
       nz = nint(ph2d(2,l))
       ph2(:,nz,ne) = ph2d(3:,l)
      enddo
      return
99    stop 'read error in file: bf_verner.dat'
      end subroutine bfxs_read_data
c
c
c
      pure function bfxs(nz,ne,e) result(xs)
c     ---------------------------
      implicit none
      integer,intent(in) :: nz,ne
      real*8,intent(in) :: e
      real*8 :: xs
****************************************************************************
* Version 2. March 25, 1996.
* Written by D. A. Verner, verner@pa.uky.edu
* Inner-shell ionization energies of some low-ionized species are slightly
* improved to fit smoothly the experimental inner-shell ionization energies
* of neutral atoms.
****************************************************************************
* This subroutine calculates partial photoionization cross sections
* for all ionization stages of all atoms from H to Zn (Z=30) by use of
* the following fit parameters:
* Outer shells of the Opacity Project (OP) elements:
*    Verner, Ferland, Korista, Yakovlev, 1996, ApJ, in press.
* Inner shells of all elements, and outer shells of the non-OP elements:
*    Verner and Yakovlev, 1995, A&AS, 109, 125
* Input parameters:  nz - atomic number from 1 to 30 (integer)
*                    ne - number of electrons from 1 to iz (integer)
*                    !is - shell number (integer)
*                    e - photon energy, eV
* Output parameter:  xs - photoionization cross section, Mb
* Shell numbers:
* 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s.
* If a species in the ground state has no electrons on the given shell,
* the subroutine returns xs=0.
****************************************************************************
      integer :: nout,nint,is
      real*8 :: p1,x,y,z,q,a,b,einn
c
      xs = 0d0
c
      if(nz<1 .or. nz>30) return
      if(ne<1 .or. ne>nz) return
c
c-- outermost shell state
      nout = ntot(ne)
c-- special cases
      if(nz>18) then
       if(nz==ne) nout = 7 !neutral K,Ca,...
       if(nz==ne+1) then
        if(nz==20 .or. nz==21 .or. nz==22 .or. !singly-ionized Ca,Sc,Ti,Mn,Fe
     &   nz==25 .or. nz==26) nout = 7
       endif
      endif
c
c     if(is>nout) return
c     if(e<ph1(1,nz,ne,is)) return
      if(e<ph1(1,nout,nz,ne)) return
c
c-- core state
      nint = ninn(ne)
      if(nz==15 .or. nz==17 .or. nz==19 .or. (nz>20 .and. nz/=26)) then
       einn = 0.0
      else
       if(ne<3) then
        einn = 1.0e30
       else
        einn = ph1(1,nint,nz,ne)
       endif
      endif
c
      xs = 0d0
      do is=1,nout
       if(e<ph1(1,is,nz,ne)) cycle
       if(is<nout .and. is>nint .and. e<einn) cycle
c--
       if(is<=nint .or. e>=einn) then
        p1 = -ph1(5,is,nz,ne)
        y = e/ph1(2,is,nz,ne)
        q = -0.5*p1 - ll(is) - 5.5
        a = ph1(3,is,nz,ne)*((y - 1.0)**2 + ph1(6,is,nz,ne)**2)
        b = sqrt(y/ph1(4,is,nz,ne)) + 1.0
        xs = xs + a*y**q*b**p1
       else
        p1 = -ph2(4,nz,ne)
        q = -0.5*p1 - 5.5
        x = e/ph2(1,nz,ne) - ph2(6,nz,ne)
        z = sqrt(x*x + ph2(7,nz,ne)**2)
        a = ph2(2,nz,ne)*((x - 1.0)**2 + ph2(5,nz,ne)**2)
        b = 1.0 + sqrt(z/ph2(3,nz,ne))
        xs = xs + a*z**q*b**p1
       endif
      enddo !is
c
      end function bfxs
c
c
c
!      subroutine write_verner
!c     -----------------------
!      implicit none
!************************************************************************
!* Write the Verner 1995 and Verner 1996 data to file.
!* The include inflates the executable.
!************************************************************************
!      real :: ph1(6,7,30,30),ph2(7,30,30) !local
!      integer :: l,ns,ne,nz,nout
!c
!      include 'bf_verner.blk'
!c
!c-- open file handles
!      open(41,file='bfxs1',action='write')
!      open(42,file='bfxs2',action='write')
!c
!      do nz=1,30
!      do ne=1,nz
!       nout = ntot(ne) !general
!c-- special cases
!       if(nz>18) then
!        if(nz==ne) nout = 7 !neutral K,Ca,...
!        if(nz==ne+1) then
!         if(nz==20 .or. nz==21 .or. nz==22 .or. !singly-ionized Ca,Sc,Ti,Mn,Fe
!     &    nz==25 .or. nz==26) nout = 7
!        endif
!       endif
!c
!       do ns=1,nout
!        write(41,'(3i3,1p,6e12.4)') ns,ne,nz,ph1(:,ns,nz,ne)
!       enddo
!       write(42,'(2i3,1p,7e12.4)') ne,nz,ph2(:,nz,ne)
!      enddo !ne
!      enddo !nz
!c
!c-- close file handles
!      close(41)
!      close(42)
!      end subroutine write_verner
c
      end module bfxsmod
