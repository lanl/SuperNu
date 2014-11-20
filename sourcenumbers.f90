subroutine sourcenumbers(nmpi)

  use totalsmod
  use gridmod
  use particlemod
  use inputparmod
  implicit none
  integer,intent(in) :: nmpi

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: i,j,k
  integer :: ns,nsfloor
  real*8 :: etot
! tot_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

!-- sanity check

  grd_nvol = 0
  grd_nvolex = 0

!-- etot
  etot = sum(grd_emit) + sum(grd_emitex) + tot_esurf
  
!-- calculating number of boundary particles (if any)
  prt_nsurf = nint(tot_esurf*prt_ns/etot)
  prt_nnew = prt_nsurf

!-- floor under number of source particles
  nsfloor = prt_ns/(grd_nx*grd_ny*grd_nz)  !uniform distribution
  nsfloor = max(1,nsfloor/10)

  ! Calculating number of particles per cell (dd_nvol): loop
  prt_nexsrc=0
  do k = 1, grd_nz
  do j = 1, grd_ny
  do i = 1, grd_nx

     !thermal volume source numbers
     if(grd_emit(i,j,k)<=0d0) then
        grd_nvol(i,j,k)=0
     else
        ns = nint(grd_emit(i,j,k)*prt_ns/etot) 
        grd_nvol(i,j,k) = max(ns,nsfloor)
     endif
     prt_nnew = prt_nnew + grd_nvol(i,j,k)

     !external volume source numbers
     if(grd_emitex(i,j,k)<=0d0) then
        grd_nvolex(i,j,k)=0
     else
        ns = nint(grd_emitex(i,j,k)*prt_ns/etot)
        grd_nvolex(i,j,k) = max(ns,nsfloor)
     endif
     prt_nexsrc = prt_nexsrc + grd_nvolex(i,j,k)
  enddo
  enddo
  enddo

!-- add to total
  prt_nnew = prt_nnew + prt_nexsrc

end subroutine sourcenumbers
