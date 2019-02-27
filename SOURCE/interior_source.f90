!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine interior_source

  use randommod
  use sourcemod
  use miscmod
  use groupmod
  use gridmod
  use totalsmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  use manufacmod
  use countersmod

  implicit none

!##################################################
  !This subroutine instantiates new volume (cell) particle properties.
  !Composed of external source particle loop (1st) and thermal source
  !particle loop (2nd).
!##################################################
  integer :: i,j,k, ipart,ivac,ig,ii
  integer :: nhere,ndmy,iimpi,nemit
  integer*8,allocatable :: nvacantall(:)
  real*8 :: pwr
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, ep0, wl0
  real*8 :: denom2,x1,x2,x3,x4, thelp
!-- neighbor emit values (for source tilting)
  integer :: icnb(6)
!
  real*8 :: emitprob(grp_ng)
  type(packet),pointer :: ptcl
!-- statement functions
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- shortcut
  pwr = in_srcepwr

  x1=grp_wlinv(grp_ng+1)
  x2=grp_wlinv(1)

  allocate(nvacantall(size(src_nvacantall)))

!Volume particle instantiation: loop
!Loop run over the number of new particles that aren't surface source
!particles.
  ipart = src_nsurf
  iimpi = -1
  nvacantall = src_nvacantall
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     call sourcenumbers_roundrobin_limit(iimpi,nvacantall,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,ndmy,nhere)
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = src_ivacant(ipart)
     ptcl => prt_particles(ivac)

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.

!
!-- calculating particle time
     call rnd_r(r1,rnd_state)
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     do ig = 1, grp_ng-1
        x3=grp_wlinv(ig+1)
        x4=grp_wlinv(ig)
        if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
        denom2 = denom2+(x4-x3)/(x2-x1)
     enddo
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_emitex(l)/dble(grd_nvol(l)-nemit)
     ptcl%e = ep0
     ptcl%e0 = ep0

!-- calculating direction cosine (comoving)
     call rnd_r(r1,rnd_state)
     mu0 = 1d0-2d0*r1
     ptcl%mu = mu0 !overwrite when isvelocity
!-- sampling azimuthal angle of direction
     call rnd_r(r1,rnd_state)
     om0 = pc_pi2*r1
     ptcl%om = om0 !overwrite when isvelocity

!
!-- x position
     select case(grd_igeom)
     case(1,11) !-- spherical
        call rnd_r(r1,rnd_state)
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
     case(2) !-- cylindrical
        call rnd_r(r1,rnd_state)
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + (1d0-r1)*grd_xarr(i)**2)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
     case(3) !-- cartesian
        call rnd_r(r1,rnd_state)
        ptcl%x = r1*grd_xarr(i+1) + (1d0-r1)*grd_xarr(i)
     endselect
!
!-- y,z position
     call rnd_r(r1,rnd_state)
     ptcl%y = r1*grd_yarr(j+1) + (1d0-r1) * grd_yarr(j)
     call rnd_r(r1,rnd_state)
     ptcl%z = r1*grd_zarr(k+1) + (1d0-r1) * grd_zarr(k)
!}}}
  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=src_nsurf+src_nnonth) stop 'interior_source: n/=nexecsrc'
  

!-- Thermal volume particle instantiation: loop
  iimpi = -1
  nvacantall = src_nvacantall
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     call sourcenumbers_roundrobin_limit(iimpi,nvacantall,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,nhere,ndmy)
     if(nhere<1) cycle
!-- integrate planck function over each group
     call specintv(grd_tempinv(l),grp_ng,emitprob)
     emitprob = emitprob*grd_cap(:,l)/grd_capgrey(l)
!
!-- neighbors
     icnb(1) = grd_icell(max(i-1,1),j,k)      !left neighbor
     icnb(2) = grd_icell(min(i+1,grd_nx),j,k) !right neighbor
     icnb(3) = grd_icell(i,max(j-1,1),k)      !left neighbor
     icnb(4) = grd_icell(i,min(j+1,grd_ny),k) !right neighbor
     icnb(5) = grd_icell(i,j,max(k-1,1))      !left neighbor
     icnb(6) = grd_icell(i,j,min(k+1,grd_nz)) !right neighbor
!
!
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = src_ivacant(ipart)
     ptcl => prt_particles(ivac)

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     call rnd_r(r1,rnd_state)
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     do ig = 1, grp_ng-1
        if (r1>=denom2.and.r1<denom2+emitprob(ig)) exit
        denom2 = denom2+emitprob(ig)
     enddo
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_emit(l)/dble(nemit)
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
!-- x position:
     select case(grd_igeom)
!!{{{
!-- spherical
     case(1,11)
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%x = x0
!
!-- cylindrical
     case(2)
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           x0 = sqrt(r1*grd_xarr(i+1)**2+(1.0-r1)*grd_xarr(i)**2)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%x = x0
!
!-- cartesian
     case(3)
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%x = r1*grd_xarr(i+1)+(1d0-r1)*grd_xarr(i)!}}}
     endselect
!-- must be inside cell
     ptcl%x = min(ptcl%x,grd_xarr(i+1))
     ptcl%x = max(ptcl%x,grd_xarr(i))
!
!-- y,z position
     if(grd_igeom==11) then
        ptcl%y = grd_yarr(1)
        ptcl%z = grd_zarr(1)
     else
!-- source tilting in y!{{{
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(3)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(4)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
!-- must be inside cell
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
!
!-- source tilting in z
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(5)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(6)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_r(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_r(r2,rnd_state)
        enddo
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))!}}}
     endif !grd_igeom
!}}}
  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=src_nnew) stop 'interior_source: n/=nnew'

  call counterreg(ct_npcreate, src_nnew)

end subroutine interior_source
! vim: fdm=marker
