!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine initialnumbers(nmpi)

  use gridmod
  use sourcemod
  use totalsmod
  use inputparmod
  implicit none
  integer,intent(in) :: nmpi

!##################################################
  !This subroutine computes the distribution of initial particle energy
  !before the first time step.  A fraction of the total initial particle
  !number is given to each cell based on the amount of inital radiative
  !energy profile.
!##################################################

  integer :: l,iimpi
  real*8 :: base
  integer :: n,ndone
  real*8,parameter :: basefrac=.1d0
  integer*8 :: nstot,nsavail,nsbase
  real*8 :: etot,einv,pwr,edone,en,invn
  integer :: nemit,nvol,nvolex,ncell
  
!-- shortcut
  pwr = in_srcepwr

!-- total particle number
  nstot = nmpi*int(src_ninit,8)
!
  call analytic_initial

  etot = sum(grd_evolinit**pwr)
  if(etot/=etot) stop 'initialnumbers: etot nan'
!
  tot_eext = etot
  
!-- base (flat,constant) particle number per cell over ALL RANKS
  n = count(grd_evolinit>0d0)  !number of cells with nonzero energy
  base = dble(nstot)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase
  
!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/(nsavail + nsbase)
  einv = 1d0/etot
  do l=1,grd_ncell
     en = grd_evolinit(l)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     grd_nvolinit(l) = n
     edone = edone + en
     ndone = ndone + n
  enddo

!-- from total nvol (over ALL RANKS) to nvol PER RANK
!-- also convert emit to energy PER PARTICLE
  src_ninitnew = 0
  iimpi = 0
  ncell=grd_ncell
  if(grd_ivoid>0) ncell=ncell-1
  do l=1,ncell
     call sourcenumbers_roundrobin(iimpi,grd_evolinit(l)**pwr, &
        0.0d0,grd_nvolinit(l),nemit,nvol,nvolex)
!-- particle counts
     src_ninitnew = src_ninitnew + nvol
  enddo

end subroutine initialnumbers
! vim: fdm=marker
