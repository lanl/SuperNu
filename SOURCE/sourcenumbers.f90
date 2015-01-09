subroutine sourcenumbers(nmpi)
!{{{
  use sourcemod
  use totalsmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  integer,intent(in) :: nmpi

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number src_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: l,iimpi
  real*8 :: base
  integer :: n,ndone
  real*8,parameter :: basefrac=.1d0
  integer*8 :: nstot,nsavail,nsbase
  real*8 :: etot,einv,pwr,edone,en,invn
  integer :: nemit,nvol,nvolex
! tot_esurf for any new prt_particles from a surface source
! src_nsurf = number of surface prt_particles
! src_nnew = total number of new prt_particles~=src_ns

!-- initialize volume numbers
  grd_nvol = 0

!-- shortcut
  pwr = in_srcepwr

!-- total particle number
  nstot = nmpi*int(src_ns,8)

!-- etot
  etot = sum(grd_emit**pwr) + sum(grd_emitex**pwr) + tot_esurf**pwr
  if(etot/=etot) stop 'sourcenumber: etot nan'

!-- calculating number of boundary particles (if any)
  src_nsurf = nint(tot_esurf**pwr*nstot/etot)

!-- base (flat,constant) particle number per cell over ALL RANKS
  n = count(grd_emit>0d0 .or. grd_emitex>0d0)  !number of cells with nonzero energy
  base = dble(nstot - src_nsurf)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase - src_nsurf


!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/(nsavail + nsbase)
  einv = 1d0/etot
  do l=1,grd_ncell
     en = grd_emit(l)**pwr + grd_emitex(l)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     grd_nvol(l) = n
     edone = edone + en
     ndone = ndone + n
  enddo


!-- from total nvol (over ALL RANKS) to nvol PER RANK
!-- also convert emit to energy PER PARTICLE
  src_nnew = src_nsurf
  src_nnonth = 0
  iimpi = 0
  do l=1,grd_ncell
     call sourcenumbers_roundrobin(iimpi,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,nvol,nvolex)
!-- particle counts
     src_nnew = src_nnew + nvol + nvolex
     src_nnonth = src_nnonth + nvolex
  enddo
!}}}
end subroutine sourcenumbers



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
