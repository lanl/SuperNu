subroutine sourcenumbers
!{{{
  use mpimod
  use totalsmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: i,j,k,iimpi
  real*8 :: base
  integer :: n,ndone
  real*8,parameter :: basefrac=.1d0
  integer*8 :: nstot,nsavail,nsbase
  real*8 :: etot,einv,pwr,edone,en,invn
  integer :: nemit,nvol,nvolex
! tot_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

!-- initialize volume numbers
  grd_nvol=0

!-- shortcut
  pwr = in_srcepwr

!-- total particle number
  nstot = nmpi*int(prt_ns,8)

!-- etot
  etot = sum(grd_emit**pwr) + sum(grd_emitex**pwr) + tot_esurf**pwr

!-- calculating number of boundary particles (if any)
  prt_nsurf = nint(tot_esurf**pwr*nstot/etot)

!-- base (flat,constant) particle number per cell over ALL RANKS
  n = count(grd_emit>0d0 .or. grd_emitex>0d0)  !number of cells with nonzero energy
  base = dble(nstot - prt_nsurf)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase - prt_nsurf


!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/(nsavail + nsbase)
  einv = 1d0/etot
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     en = grd_emit(i,j,k)**pwr + grd_emitex(i,j,k)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     grd_nvol(i,j,k) = n
     edone = edone + en
     ndone = ndone + n
  enddo
  enddo
  enddo


!-- from total nvol (over ALL RANKS) to nvol PER RANK
!-- also convert emit to energy PER PARTICLE
  prt_nnew = prt_nsurf
  prt_nexsrc = 0
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     call sourcenumbers_roundrobin(iimpi,grd_emit(i,j,k)**pwr, &
        grd_emitex(i,j,k)**pwr,grd_nvol(i,j,k),nemit,nvol,nvolex)
!-- particle counts
     prt_nnew = prt_nnew + nvol + nvolex
     prt_nexsrc = prt_nexsrc + nvolex
  enddo
  enddo
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
