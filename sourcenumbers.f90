subroutine sourcenumbers(impi,nmpi)

  use totalsmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  integer,intent(in) :: impi,nmpi

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: i,j,k,l,iimpi
  integer :: ibase,nsbase,n
  integer*8 :: nsavail
  real*8 :: esqrt,pwr,help
  integer :: nvol,nvolex
! tot_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

!-- shortcut
  pwr = in_srcepwr

!-- esqrt
  esqrt = sum(grd_emit**pwr) + sum(grd_emitex**pwr) + tot_esurf**pwr
  
!-- calculating number of boundary particles (if any)
  prt_nsurf = nint(tot_esurf**pwr*prt_ns/esqrt)

!-- base (flat,constant) particle number per cell over ALL RANKS
  ibase = (nmpi*int(prt_ns,8))/(grd_nx*grd_ny*grd_nz)  !uniform distribution
  ibase = max(10,ibase/10)

!-- number of particles available for proportional distribution
  n = count(grd_emit>0d0 .or. grd_emitex>0d0)  !number of cells with nonzero energy
  nsbase = n*ibase  !total number of base particles
  nsavail = nmpi*prt_ns - nsbase

!-- total particle number per cell
  grd_nvol = nint((grd_emit**pwr + grd_emitex**pwr)*nsavail/esqrt) + ibase


!-- from total nvol (over ALL RANKS) to nvol PER RANK
!-- also convert emit to energy PER PARTICLE
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     help = grd_emit(i,j,k)**pwr/(grd_emit(i,j,k)**pwr + grd_emitex(i,j,k)**pwr)
     nvol = nint(grd_nvol(i,j,k)*help)
     nvolex = grd_nvol(i,j,k) - nvol
!
!-- energy per particle
     grd_emit(i,j,k) = grd_emit(i,j,k)/nvol
     grd_emitex(i,j,k) = grd_emitex(i,j,k)/nvolex
     if(nvol==0) grd_emit(i,j,k) = 0
     if(nvolex==0) grd_emitex(i,j,k) = 0
!
!-- overwrite grd_nvol below
!-- insert the constant part
     grd_nvol(i,j,k) = nvol/nmpi
     grd_nvolex(i,j,k) = nvolex/nmpi
!-- remaining particles are distributed round robin over ranks
     n = mod(nvol,nmpi)
     do l=1,n
        iimpi = iimpi + 1
        if(iimpi>nmpi) iimpi = 1
        if(iimpi==impi) grd_nvol(i,j,k) = grd_nvol(i,j,k) + 1
     enddo
!-- remaining particles are distributed round robin over ranks
     n = mod(nvolex,nmpi)
     do l=1,n
        iimpi = iimpi + 1
        if(iimpi>nmpi) iimpi = 1
        if(iimpi==impi) grd_nvolex(i,j,k) = grd_nvolex(i,j,k) + 1
     enddo
  enddo
  enddo
  enddo

  prt_nnew = prt_nnew + sum(grd_nvol)
  prt_nexsrc = sum(grd_nvolex)

end subroutine sourcenumbers
