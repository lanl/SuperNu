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
  integer :: ig, i,j,k,l, iig, ipart, ii
  integer :: nhere,ndmy,iimpi,nemit
  integer :: ix,iy,iz
  real*8 :: wl0, mu0, om0, ep0
  real*8 :: denom2
  real*8 :: pwr
  real*8 :: r1
  real*8 :: wl1,wl2,wl3,wl4
  real*8 :: specarr(grp_ng)

!-- helper quantities
  wl1=grp_wlinv(grp_ng+1)
  wl2=grp_wlinv(1)
  if(in_tempradinit>0d0) then
     specarr = specintv(1d0/in_tempradinit,1)
     specarr = specarr/sum(specarr)
  else
     specarr = 1d0
  endif
!-- shortcut
  pwr = in_srcepwr
  
!-- instantiating initial particles
  ipart = 0
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     call sourcenumbers_roundrobin(iimpi,grd_evolinit(l)**pwr, &
        0.0d0,grd_nvolinit(l),nemit,nhere,ndmy)
  do ii=1,nhere
     ipart = ipart + 1!{{{
!
!-- setting 1st cell index
     ix = i

!-- setting particle index to not vacant
     prt_isvacant(ipart) = .false.

!-- particle origin
     ptcl%icorig = l
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
        if(r1>=denom2.and.r1<denom2 + specarr(ig)) exit
        denom2 = denom2 + specarr(ig)
     enddo
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_evolinit(l)/dble(nemit)
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
        ptcl%z = grd_zarr(1)

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
        ptcl%y = grd_yarr(1)
        ptcl%z = grd_zarr(1)
     endselect

!-- save particle result
!-----------------------
     prt_particles(ipart) = ptcl

  enddo!}}} !ipart
!
  enddo !i
  enddo !j
  enddo !k
  
end subroutine initial_particles
