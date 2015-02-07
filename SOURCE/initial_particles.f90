subroutine initial_particles

  use randommod
  use sourcemod
  use gridmod
  use groupmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine instantiates the initial particles before
  !the first time step.
!##################################################
  type(packet) :: ptcl
!
  logical :: lhelp
  integer :: ig, i,j,k,l, iig, ipart
  integer :: ix,iy,iz
  integer :: iused
  real*8 :: wl0, mu0, om0, ep0
  real*8 :: denom2
  real*8 :: r1
  real*8 :: wl1,wl2,wl3,wl4

!-- helper quantities
  wl1=grp_wlinv(grp_ng+1)
  wl2=grp_wlinv(1)

!-- instantiating initial particles
  i = 1
  j = 1
  k = 1
  iused = 0
  l = grd_icell(i,j,k)
  do ipart=1,src_ninitnew

!-----------------------------------------------------------------------
!TODO: make this an outer loop, like in interior_source
!-- incrementing to next vacant cell
     loopk: do k=k,grd_nz
        do j=j,grd_ny
           do i=i,grd_nx
              l = grd_icell(i,j,k)
              lhelp = iused<grd_nvolinit(l)
              if(lhelp) then
                 iused = iused + 1
                 exit loopk
              endif
              iused = 0
           enddo
           iused = 0
           i = 1
        enddo
        iused = 0
        j = 1
     enddo loopk
!-----------------------------------------------------------------------
!
!-- sanity check
     if(.not.lhelp) stop 'initial_particles: invalid particle'
!
!-- setting 1st cell index
     ix = i

!-- setting particle index to not vacant
     prt_isvacant(ipart) = .false.
!
!-- calculating particle time
     ptcl%t = tsp_t

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     do ig = 1, grp_ng
        wl3 = grp_wlinv(ig+1)
        wl4 = grp_wlinv(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2 + (wl4-wl3)/(wl2-wl1)) exit
        denom2 = denom2 + (wl4-wl3)/(wl2-wl1)
     enddo
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_evolinit(l)/dble(grd_nvolinit(l))
     ptcl%e = ep0
     ptcl%e0 = ep0

!-- calculating direction cosine (comoving)
     call rnd_r(r1,rnd_state)
     mu0 = 1d0-2d0*r1
     ptcl%mu = mu0

!-- sampling azimuthal angle of direction
     call rnd_r(r1,rnd_state)
     om0 = pc_pi2*r1
     ptcl%om = om0

!
!-- selecting geometry
     select case(in_igeom)

!-- 3D spherical
     case(1)
!-- calculating position
        call rnd_r(r1,rnd_state)
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)
        call rnd_r(r1,rnd_state)
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        call rnd_r(r1,rnd_state)
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))

!-- 2D
     case(2)
!-- setting 2nd cell index
        iy = j
!-- calculating position
        call rnd_r(r1,rnd_state)
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
        call rnd_r(r1,rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + &
             (1d0-r1)*grd_yarr(j)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))

!-- 3D
     case(3)
!-- setting 2nd, 3rd cell index
        iy = j
        iz = k
!-- calculating position
        call rnd_r(r1,rnd_state)
        ptcl%x = r1*grd_xarr(i+1) + &
             (1d0-r1)*grd_xarr(i)
        call rnd_r(r1,rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + &
             (1d0-r1)*grd_yarr(j)
        call rnd_r(r1,rnd_state)
        ptcl%z = r1*grd_zarr(k+1) + &
             (1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))

!-- 1D
     case(11)
!-- calculating position
        call rnd_r(r1,rnd_state)
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
     endselect

!-- save particle result
!-----------------------
     prt_particles(ipart) = ptcl

  enddo !ipart


end subroutine initial_particles
