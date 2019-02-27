!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine sourcenumbers(keephigh)
!{{{
  use sourcemod
  use totalsmod
  use gridmod
  use inputparmod
  use mpimod
  implicit none
  logical,intent(in) :: keephigh

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: l,iimpi
  integer*8 :: n,ndone
  integer :: nextra,nsmean
  integer*8 :: nvacant(nmpi)
  real*8,parameter :: basefrac=.1d0
  integer*8 :: nstot,nstotd,nsavail
  integer*8 :: nvacantall,nnewvacant
  integer*8 :: ncactive
  real*8 :: etot,einv,pwr,edone,en,invn
  integer :: nemit,nvol,nvolex
! tot_esurf for any new prt_particles from a surface source
! src_nsurf = number of surface prt_particles
! src_nnew = total number of new prt_particles~=src_ns

!-- initialize volume numbers
  grd_nvol = 0

!-- shortcut
  pwr = in_srcepwr

!-- total number of vacancies on all ranks
  nvacantall = sum(src_nvacantall)
!-- number of new vacancies
  nnewvacant = nvacantall - src_nvacantmin
!-- keep record of maximum number of active particles
  src_nvacantmin = min(src_nvacantmin,nvacantall)

!-- total particle number
  nstot = nmpi*src_ns
!-- keep at least at level of minimum vacancies
  if(keephigh) nstot = max(nmpi*src_ns,nnewvacant,src_nflux)
!-- limit total particle number
  nstotd = nstot - src_nflux
  nstotd = min(nstotd,(nvacantall-src_nflux)/2 + 1)
  nstot = src_nflux + nstotd
!-- limit total particle number
  nstot = min(nstot,nvacantall)

!-- etot
  etot = sum(grd_emit**pwr) + sum(grd_emitex**pwr) + tot_esurf**pwr
  if(etot/=etot) stop 'sourcenumbers: etot nan'

!-- number of boundary particles (if any)
  src_nsurftot = nint(tot_esurf**pwr*nstot/etot)

!-- number of cells with nonzero energy
  ncactive = count(grd_emit>0d0 .or. grd_emitex>0d0)
  if(src_nsurftot>0) ncactive = ncactive+1
  if(nvacantall<ncactive) stop &
     'sourcenumbers: insufficient vacancies in particle array'

!-- number of particles available for proportional distribution
  nsavail = nstot - src_nsurftot

!-- particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/nsavail
  einv = 1d0/etot
  if(etot==0d0) einv = 0d0
  do l=1,grd_ncell
     en = grd_emit(l)**pwr + grd_emitex(l)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = max(1,int(en*nsavail*einv))  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     grd_nvol(l) = int(n)
     edone = edone + en
     ndone = ndone + n
  enddo
  if(ndone==0) stop 'sourcenumbers: no volume source particles distributed'

!-- too many particles
  nextra = int(nstot - sum(grd_nvol) - src_nsurftot)
  if(nextra>grd_ncell+src_nsurftot) stop &
       'sourcenumbers: nextra>grd_ncell+src_nsurftot'
!-- correct to exact target number
  nsmean = int(nstot/ncactive)
  do while(nextra/=0)
     if(src_nsurftot>=nsmean) then
        src_nsurftot = src_nsurftot + sign(1,nextra)
        nextra = nextra - sign(1,nextra)
     endif
     do l=1,grd_ncell
        if(nextra==0) exit
        if(grd_nvol(l)<nsmean) cycle
        grd_nvol(l) = grd_nvol(l) + sign(1,nextra)
        nextra = nextra - sign(1,nextra)
     enddo
  enddo

!-- from total nvol (over ALL RANKS) to nvol PER RANK
!-- also convert emit to energy PER PARTICLE
  src_nsurf = 0
  iimpi = -1
  if(src_nsurftot>0) then
!-- from total nsurf to nsurf per rank
     call sourcenumbers_roundrobin_limit(iimpi,src_nvacantall, &
          tot_esurf**pwr,0d0,src_nsurftot,nemit,nvol,nvolex)
     src_nsurf = nvol
  endif

  src_nnew = src_nsurf
  src_nnonth = 0
  nvacant = src_nvacantall
  iimpi = -1
  do l=1,grd_ncell
     call sourcenumbers_roundrobin_limit(iimpi,nvacant,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,nvol,nvolex)
!-- particle counts
     src_nnew = src_nnew + nvol + nvolex
     src_nnonth = src_nnonth + nvolex
  enddo
!}}}
end subroutine sourcenumbers



subroutine sourcenumbers_roundrobin_limit(iimpi,nvacant,evol,evolex,ntot,mvol,nvol,nvolex)
  use mpimod!{{{
  implicit none
  integer,intent(inout) :: iimpi
  integer*8,intent(inout) :: nvacant(0:nmpi-1)
  real*8,intent(in) :: evol,evolex
  integer,intent(in) :: ntot
  integer,intent(out) :: mvol  !particle number on all ranks
  integer,intent(out) :: nvol,nvolex !particle numbers on this rank
!-----------------------------------------------------------------------
! Distribute the source particle numbers over the mpi ranks in a
! round-robin fashion using iimpi as tracker.
!-----------------------------------------------------------------------
  real*8 :: help
  integer :: i,n,l
  integer :: nmpiavail,neach,nhere
!
!-- quick exit
  if(evol+evolex==0d0) then
     mvol = 0
     nvol = 0
     nvolex = 0
     return
  endif
!
!-- 
  help = evol/(evol + evolex)
  mvol = nint(ntot*help)
  nvolex = ntot - mvol
  n = mod(mvol,nmpi) !remainder
!-- constant part
  nvol = mvol/nmpi !on each rank
!-- add to remainder whatever does not fit in vacancies
  do l=0,nmpi-1
     n = n + max(0,nvol-int(nvacant(l)))
     if(l==impi) nvol = min(nvol,int(nvacant(l))) !current rank
     nvacant(l) = max(0,nvacant(l)-nvol)
  enddo
!-- sanity check
  if(n>sum(nvacant)) stop 'srcnr_roundrobin_limit: not enough vacancies'
!-- remaining particles are distributed round robin over ranks
  nmpiavail = count(nvacant>0)
  do while(n>0)
     neach = max(1,n/nmpiavail) !number of particles to add to this rank
     do i=1,nmpi
        iimpi = iimpi + 1
        if(iimpi==nmpi) iimpi = 0
!-- no space on current rank
        if(nvacant(iimpi)==0) cycle
!-- number of particles that can be added to this rank
        nhere = min(neach,int(nvacant(iimpi)))
        if(iimpi==impi) nvol = nvol + nhere  !current rank
!-- remaining
        nvacant(iimpi) = nvacant(iimpi) - nhere
        n = n - nhere
        if(n==0) exit
        if(nvacant(iimpi)==0) nmpiavail = nmpiavail - 1
     enddo
  enddo

  n = mod(nvolex,nmpi) !remainder
!-- constant part
  nvolex = nvolex/nmpi !on each rank
!-- add to remainder whatever does not fit in vacancies
  do l=0,nmpi-1
     n = n + max(0,nvolex-int(nvacant(l)))
     if(l==impi) nvolex = min(nvolex,int(nvacant(l))) !current rank
     nvacant(l) = max(0,nvacant(l)-nvolex)
  enddo
!-- sanity check
  if(n>sum(nvacant)) stop 'srcnr_roundrobin_limit: not enough vacancies'
!-- remaining particles are distributed round robin over ranks
  nmpiavail = count(nvacant>0)
  do while(n>0)
     neach = max(1,n/nmpiavail) !number of particles to add to this rank
     do i=1,nmpi
        iimpi = iimpi + 1
        if(iimpi==nmpi) iimpi = 0
!-- no space on current rank
        if(nvacant(iimpi)==0) cycle
!-- number of particles that can be added to this rank
        nhere = min(neach,int(nvacant(iimpi)))
        if(iimpi==impi) nvolex = nvolex + nhere  !current rank
!-- remaining
        nvacant(iimpi) = nvacant(iimpi) - nhere
        n = n - nhere
        if(n==0) exit
        if(nvacant(iimpi)==0) nmpiavail = nmpiavail - 1
     enddo
  enddo
!}}}
end subroutine sourcenumbers_roundrobin_limit



subroutine sourcenumbers_roundrobin(iimpi,evol,evolex,ntot,mvol,nvol,nvolex)
  use mpimod!{{{
  implicit none
  integer,intent(inout) :: iimpi
  real*8,intent(in) :: evol,evolex
  integer,intent(in) :: ntot
  integer,intent(out) :: mvol  !particle number on all ranks
  integer,intent(out) :: nvol,nvolex !particle numbers on this rank
!-----------------------------------------------------------------------
! Distribute the source particle numbers over the mpi ranks in a
! round-robin fashion using iimpi as tracker.
!-----------------------------------------------------------------------
  real*8 :: help
  integer :: n,l
!
!-- quick exit
  if(evol+evolex==0d0) then
     mvol = 0
     nvol = 0
     nvolex = 0
     return
  endif
!
!-- 
  help = evol/(evol + evolex)
  mvol = nint(ntot*help)
  nvolex = ntot - mvol
  n = mod(mvol,nmpi) !remainder
!-- constant part
  nvol = mvol/nmpi !on each rank
!-- remaining particles are distributed round robin over ranks
  do l=1,n
     if(iimpi==impi) nvol = nvol + 1
     iimpi = iimpi + 1
     if(iimpi==nmpi) iimpi = 0
  enddo
  n = mod(nvolex,nmpi) !remainder
!-- constant part
  nvolex = nvolex/nmpi !on each rank
!-- remaining particles are distributed round robin over ranks
  do l=1,n
     if(iimpi==impi) nvolex = nvolex + 1
     iimpi = iimpi + 1
     if(iimpi==nmpi) iimpi = 0
  enddo
!}}}
end subroutine sourcenumbers_roundrobin
! vim: fdm=marker
