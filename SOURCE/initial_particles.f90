!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
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
  real*8 :: wl0, mu0, om0, ep0
  real*8 :: denom2
  real*8 :: pwr
  real*8 :: r1
  real*8 :: tradinv
  real*8 :: specarr(grp_ng)

!-- shortcut
  pwr = in_srcepwr
  
!-- instantiating initial particles
  ipart = 0
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     if(l==grd_ivoid) cycle
     call sourcenumbers_roundrobin(iimpi,grd_evolinit(l)**pwr, &
        0.0d0,grd_nvolinit(l),nemit,nhere,ndmy)
     tradinv = (pc_acoef*grd_vol(l)/grd_evolinit(l))**(.25d0)
     call specintv(tradinv,grp_ng,specarr,1)
     specarr = specarr/sum(specarr)
  do ii=1,nhere
     ipart = ipart + 1!{{{
!
!-- setting particle index to not vacant
     prt_isvacant(ipart) = .false.
!
!-- calculating particle time
     ptcl%t = tsp_t

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     do ig = 1, grp_ng
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
!-- x position
     select case(grd_igeom)
!-- spherical
     case(1,11)
        call rnd_r(r1,rnd_state)
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- cylindrical
     case(2)
        call rnd_r(r1,rnd_state)
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + (1d0-r1)*grd_xarr(i)**2)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- cartesian
     case(3)
        call rnd_r(r1,rnd_state)
        ptcl%x = r1*grd_xarr(i+1) + (1d0-r1)*grd_xarr(i)
     endselect
!
!-- x,y position
     call rnd_r(r1,rnd_state)
     ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
     call rnd_r(r1,rnd_state)
     ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)

!-- save particle result
!-----------------------
     prt_particles(ipart) = ptcl

  enddo!}}} !ipart
!
  enddo !i
  enddo !j
  enddo !k
  
end subroutine initial_particles
! vim: fdm=marker
