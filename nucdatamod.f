*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module nucdatamod
c     -----------------
      use physconstmod, only:pc_kev
      private pc_kev
c
c-- nuclear data: http://ie.lbl.gov/
c-- nuclear decay half-life times
      real*8,parameter :: nuc_thl_ni56 = 7.575d5 !sec, 6.077 days, converted to 1/exp life
      real*8,parameter :: nuc_thl_co56 = 9.632d6 !sec, 77.27 days, converted to 1/exp life
      real*8,parameter :: nuc_thl_fe52 = 4.298d4 !sec, 8.275(8) h, converted to 1/exp life
      real*8,parameter :: nuc_thl_mn52 = 1.828d3 !sec, 21.1(2) min, converted to 1/exp life
      real*8,parameter :: nuc_thl_cr48 = 1.119d5 !sec, 21.56 h, converted to 1/exp life
      real*8,parameter :: nuc_thl_v48 = 1.991d6 !sec, 15.973 days, converted to 1/exp life
c
c-- nuclear decay gamma emission, 1MeV=1.602176d-6 erg
      real*8,parameter :: nuc_qhl_ni56 = 1720d0*pc_kev                  !total gamma ray production
      real*8,parameter :: nuc_qhl_co56 = (3440d0 + .19d0*2*511d0)*pc_kev!direct gamma's plus positron annihilation
c-- average kinetic energy of co56->fe56 emergent positron
      real*8,parameter :: nuc_q_poskin = .19d0*116*pc_kev
c
      contains
c
c
c
      subroutine nucdecay3(ncell,t,thl3,thl2,x)
c     -------------------------------------------
      implicit none
      integer,intent(in) :: ncell
      real*8,intent(in) :: t,thl3,thl2    !half-life of ni56 and co56 (or other 3-level chains)
      real*8,intent(inout) :: x(ncell,3)  !(fe56,co56,ni56)
************************************************************************
* Calculate x 3-level decay chain, such as ni56->co56->fe56
* x input: ni56/co56 at t=0
* x output: ni56/co56/fe56 at t=t
************************************************************************
      real*8 :: help,exp3,exp2
c
      exp3 = exp(-t/thl3) !ni56
      exp2 = exp(-t/thl2) !co56
c
c-- calculate fe
      help = 1d0 + (thl2*exp2 - thl3*exp3)/(thl3 - thl2)
      if(help.lt.0) stop 'nucdecay3: ni->fe < 0'
      x(:,1) = x(:,3)*help + x(:,2)*(1d0-exp2) !initial ni56 + initial co56
c
c-- update co56
      help = thl2*(exp3 - exp2)/(thl3 - thl2)
      if(help.lt.0) stop 'nucdecay3: ni->co < 0'
      x(:,2) = x(:,3)*help + x(:,2)*exp2 !initial ni56 + initial co56
c
c-- update ni56
      x(:,3) = x(:,3)*exp3  !initial ni56
      end subroutine nucdecay3
c
      end module nucdatamod
c vim: fdm=marker
