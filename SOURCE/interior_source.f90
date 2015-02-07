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

  implicit none

!##################################################
  !This subroutine instantiates new volume (cell) particle properties.
  !Composed of external source particle loop (1st) and thermal source
  !particle loop (2nd).
!##################################################
  integer :: i,j,k, ipart,ivac,ig,ii
  integer :: nhere,ndmy,iimpi,nemit
  real*8 :: pwr
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, y0, ep0, wl0
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

!Volume particle instantiation: loop
!Loop run over the number of new particles that aren't surface source
!particles.
  ipart = src_nsurf
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     call sourcenumbers_roundrobin(iimpi,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,ndmy,nhere)
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = src_ivacant(ipart)
     ptcl => prt_particles(ivac)

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.

!
!-- calculating particle time
     call rnd_rp(r1,rnd_state)
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     call rnd_rp(r1,rnd_state)
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grp_ng-1
        x3=grp_wlinv(ig+1)
        x4=grp_wlinv(ig)
        if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
        denom2 = denom2+(x4-x3)/(x2-x1)
     enddo
     call rnd_rp(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_emitex(l)/dble(grd_nvol(l)-nemit)
     ptcl%e = ep0
     ptcl%e0 = ep0

!-- calculating direction cosine (comoving)
     call rnd_rp(r1,rnd_state)
     mu0 = 1d0-2d0*r1
     ptcl%mu = mu0 !overwrite when isvelocity

!-- sampling azimuthal angle of direction
     call rnd_rp(r1,rnd_state)
     om0 = pc_pi2*r1
     ptcl%om = om0 !overwrite when isvelocity

!
!-- selecting geometry
     select case(in_igeom)

!-- 3D spherical
     case(1)
!-- calculating position!{{{
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))!}}}
!-- 2D
     case(2)
!-- calculating position!{{{
        call rnd_rp(r1,rnd_state)
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
        call rnd_rp(r1,rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))!}}}
!-- 3D
     case(3)
!-- calculating position!{{{
        call rnd_rp(r1,rnd_state)
        ptcl%x = r1*grd_xarr(i+1) + (1d0-r1) * &
             grd_xarr(i)
        call rnd_rp(r1,rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
        call rnd_rp(r1,rnd_state)
        ptcl%z = r1*grd_zarr(k+1) + (1d0-r1) * &
             grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k)) !}}}
!-- 1D
     case(11)
!-- calculating position!{{{
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i)) !}}}

     endselect
!}}}
  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=src_nsurf+src_nnonth) stop 'interior_source: n/=nexecsrc'
  

!-- Thermal volume particle instantiation: loop
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
     call sourcenumbers_roundrobin(iimpi,grd_emit(l)**pwr, &
        grd_emitex(l)**pwr,grd_nvol(l),nemit,nhere,ndmy)
     if(nhere<1) cycle
!-- integrate planck function over each group
     emitprob = specintv(1d0/grd_temp(l),0)
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
     call rnd_rp(r1,rnd_state)
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     call rnd_rp(r1,rnd_state)
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grp_ng-1
        if (r1>=denom2.and.r1<denom2+emitprob(ig)) exit
        denom2 = denom2+emitprob(ig)
     enddo
     call rnd_rp(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
     ptcl%wl = wl0

!-- calculating particle energy
     ep0 = grd_emit(l)/dble(nemit)
     ptcl%e = ep0
     ptcl%e0 = ep0

!-- calculating direction cosine (comoving)
     call rnd_rp(r1,rnd_state)
     mu0 = 1d0-2d0*r1
     ptcl%mu = mu0

!-- sampling azimuthal angle of direction
     call rnd_rp(r1,rnd_state)
     om0 = pc_pi2*r1
     ptcl%om = om0

!
!-- selecting geometry
     select case(in_igeom)

!-- 3D spherical
     case(1)
!-- calculating position:!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           prt_tlyrand = prt_tlyrand+1
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_rp(r2,rnd_state)
           prt_tlyrand = prt_tlyrand+1
        enddo
        ptcl%x = x0
!-- uniform in angles
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        call rnd_rp(r1,rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))
!}}}
!-- 2D
     case(2)
!-- calculating position:!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           x0 = sqrt(r1*grd_xarr(i+1)**2+(1.0-r1)*grd_xarr(i)**2)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_rp(r2,rnd_state)
        enddo
        ptcl%x = x0
!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(3)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(4)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_rp(r2,rnd_state)
        enddo
        y0 = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        ptcl%y = y0
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
!}}}
!-- 3D
     case(3)
!-- source tilting in x!{{{
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_rp(r2,rnd_state)
        enddo
        ptcl%x = r1*grd_xarr(i+1)+(1d0-r1)*grd_xarr(i)

!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(3)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(4)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_rp(r2,rnd_state)
        enddo
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)

!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(5)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(6)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           r3 = r1*uur+(1d0-r1)*uul
           call rnd_rp(r2,rnd_state)
        enddo
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))
!}}}
!-- 1D
     case(11)
!-- calculating position!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        uul = .5d0*(grd_emit(icnb(1)) + grd_emit(l))
        uur = .5d0*(grd_emit(icnb(2)) + grd_emit(l))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           call rnd_rp(r1,rnd_state)
           prt_tlyrand = prt_tlyrand+1
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           call rnd_rp(r2,rnd_state)
           prt_tlyrand = prt_tlyrand+1
        enddo
        ptcl%x = x0
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!}}}
     endselect
!}}}
  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=src_nnew) stop 'interior_source: n/=nnew'


end subroutine interior_source
