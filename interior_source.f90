subroutine interior_source

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
  logical :: lhelp
  integer :: i,j,k,il,ir,ipart,ivac,ig,ii
  integer :: nhere,ndmy,iimpi,nemit
  real*8 :: pwr
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, y0, z0, ep0, wl0
  real*8 :: denom2,x1,x2,x3,x4, thelp
  real*8 :: cmffact,mu1,mu2,gm
  real*8 :: emitprob(grp_ng)
  type(packet),target :: ptcl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)

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
  ipart = prt_nsurf
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     call sourcenumbers_roundrobin(iimpi,grd_emit(i,j,k)**pwr, &
        grd_emitex(i,j,k)**pwr,grd_nvol(i,j,k),nemit,ndmy,nhere)
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = prt_vacantarr(ipart)

!-- setting cell index
     ptcl%ix = i
     ptcl%iy = j
     ptcl%iz = k

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.

!-- default, recalculated for isvelocity and itype==1
     cmffact = 1d0
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grp_ng-1
        x3=grp_wlinv(ig+1)
        x4=grp_wlinv(ig)
        if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
        denom2 = denom2+(x4-x3)/(x2-x1)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = grd_emitex(i,j,k)/dble(grd_nvol(i,j,k)-nemit)

!
!-- selecting geometry
     select case(in_igeom)

!-- 3D spherical
     case(1)
!-- calculating position!{{{
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- sampling azimuthal angle of direction
        r1 = rand()
        ptcl%om = pc_pi2*r1
!-- setting IMC logical
        lhelp = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k)) * &
             min(dx(i),xm(i)*dyac(j),xm(i)*ym(j)*dz(k)) * &
             thelp < prt_tauddmc).or.(in_puretran)

!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
           ptcl%itype = 1 !IMC
        else
           ptcl%itype = 2 !DDMC
           ptcl%mu = mu0
        endif
!}}}
!-- 2D
     case(2)
!-- calculating position!{{{
        r1 = rand()
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
        r1 = rand()
        ptcl%y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- setting IMC logical
        lhelp = ((grd_sig(i,j,1)+grd_cap(ig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           gm = 1d0/sqrt(1d0-(x0**2+y0**2)/pc_c**2)
!-- om
           ptcl%om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+(gm*x0/pc_c) * &
                (1d0+gm*(cmffact-1d0)/(gm+1d0)))
           if(ptcl%om<0d0) ptcl%om=ptcl%om+pc_pi2
!-- mu
           ptcl%mu = (mu0+(gm*y0/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                (gm*cmffact)
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%om = om0
           ptcl%itype = 2 !DDMC
        endif
!}}}
!-- 3D
     case(3)
!-- calculating position!{{{
        r1 = rand()
        ptcl%x = r1*grd_xarr(i+1) + (1d0-r1) * &
             grd_xarr(i)
        r1 = rand()
        ptcl%y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
        r1 = rand()
        ptcl%z = r1*grd_zarr(k+1) + (1d0-r1) * &
             grd_zarr(k)
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- setting IMC logical
        lhelp = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k)) * &
             min(dx(i),dy(j),dz(k))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
           z0 = ptcl%z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           ptcl%mu = (mu0+z0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           ptcl%om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(ptcl%om<0d0) ptcl%om = ptcl%om+pc_pi2
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%om = om0
           ptcl%itype = 2 !DDMC
        endif
!}}}
!-- 1D
     case(4)
!-- calculating position!{{{
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- setting IMC logical
        lhelp = ((grd_sig(i,1,1)+grd_cap(ig,i,1,1))*dx(i)* &
             thelp < prt_tauddmc).or.(in_puretran)

!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
           ptcl%itype = 1 !IMC
        else
           ptcl%itype = 2 !DDMC
           ptcl%mu = mu0
        endif
!}}}
     endselect

     ptcl%e = ep0*cmffact
     ptcl%e0 = ep0*cmffact
     ptcl%wl = wl0/cmffact
!-- velocity effects accounting
     tot_evelo=tot_evelo+ep0*(1d0-cmffact)

!-- save particle result
!-----------------------
     prt_particles(ivac) = ptcl
!}}}
  enddo !ipart
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=prt_nsurf+prt_nexsrc) stop 'interior_source: n/=nexecsrc'
  

!-- Thermal volume particle instantiation: loop
  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     call sourcenumbers_roundrobin(iimpi,grd_emit(i,j,k)**pwr, &
        grd_emitex(i,j,k)**pwr,grd_nvol(i,j,k),nemit,nhere,ndmy)
     if(nhere<1) cycle
!-- integrate planck function over each group
     emitprob = specintv(1d0/grd_temp(i,j,k),0)
     emitprob = emitprob*grd_cap(:,i,j,k)/grd_capgrey(i,j,k)
!
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = prt_vacantarr(ipart)
!
!-- setting cell index
     ptcl%ix = i
     ptcl%iy = j
     ptcl%iz = k

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.

!-- default, recalculated for isvelocity and itype==1
     cmffact = 1d0

!-- default IMC, reset if DDMC
     ptcl%itype = 1
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1     
     do ig = 1, grp_ng-1
        if (r1>=denom2.and.r1<denom2+emitprob(ig)) exit
        denom2 = denom2+emitprob(ig)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = grd_emit(i,j,k)/dble(nemit)

!
!-- selecting geometry
     select case(in_igeom)

!-- 1D
     case(1)
!-- calculating position:!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,grd_nx)  !-- right neighbor
        uul = .5d0*(grd_emit(il,j,k) + grd_emit(i,j,k))
        uur = .5d0*(grd_emit(ir,j,k) + grd_emit(i,j,k))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
        enddo
        ptcl%x = x0
!-- uniform in angles
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- setting IMC logical
        lhelp = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k)) * &
             min(dx(i),xm(i)*dyac(j),xm(i)*ym(j)*dz(k)) * &
             thelp < prt_tauddmc).or.(in_puretran)
!write(0,*) i,grd_sig(i,1,1),grd_cap(ig,i,1,1),dx(i),thelp,prt_tauddmc

!-- if velocity-dependent, transforming direction
        if (lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%itype = 2 !DDMC
        endif
!}}}
!-- 2D
     case(2)
!-- calculating position:!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,grd_nx)  !-- right neighbor
        uul = .5d0*(grd_emit(il,j,1) + grd_emit(i,j,1))
        uur = .5d0*(grd_emit(ir,j,1) + grd_emit(i,j,1))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           x0 = sqrt(r1*grd_xarr(i+1)**2+(1.0-r1)*grd_xarr(i)**2)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           r2 = rand()
        enddo
        ptcl%x = x0
!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        il = max(j-1,1)  !-- lower neighbor
        ir = min(j+1,grd_ny)  !-- upper neighbor
        uul = .5d0*(grd_emit(i,il,1) + grd_emit(i,j,1))
        uur = .5d0*(grd_emit(i,ir,1) + grd_emit(i,j,1))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           r3 = r1*uur+(1d0-r1)*uul
           r2 = rand()
        enddo
        y0 = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        ptcl%y = y0
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1

!-- setting IMC logical
        lhelp = ((grd_sig(i,j,1)+grd_cap(ig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           gm = 1d0/sqrt(1d0-(x0**2+y0**2)/pc_c**2)
!-- om
           ptcl%om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+(gm*x0/pc_c) * &
                (1d0+gm*(cmffact-1d0)/(gm+1d0)))
           if(ptcl%om<0d0) ptcl%om=ptcl%om+pc_pi2
!-- mu
           ptcl%mu = (mu0+(gm*y0/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                (gm*cmffact)
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%om = om0
           ptcl%itype = 2 !DDMC
        endif
!}}}
!-- 3D
     case(3)
!-- source tilting in x!{{{
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,grd_nx)  !-- right neighbor
        uul = .5d0*(grd_emit(il,j,k) + grd_emit(i,j,k))
        uur = .5d0*(grd_emit(ir,j,k) + grd_emit(i,j,k))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           r3 = r1*uur+(1d0-r1)*uul
           r2 = rand()
        enddo
        ptcl%x = r1*grd_xarr(i+1)+(1d0-r1)*grd_xarr(i)

!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        il = max(j-1,1)  !-- lower neighbor
        ir = min(j+1,grd_ny)  !-- upper neighbor
        uul = .5d0*(grd_emit(i,il,k) + grd_emit(i,j,k))
        uur = .5d0*(grd_emit(i,ir,k) + grd_emit(i,j,k))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           r3 = r1*uur+(1d0-r1)*uul
           r2 = rand()
        enddo
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)

!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        il = max(k-1,1)  !-- lower neighbor
        ir = min(k+1,grd_nz)  !-- upper neighbor
        uul = .5d0*(grd_emit(i,j,il) + grd_emit(i,j,k))
        uur = .5d0*(grd_emit(i,j,ir) + grd_emit(i,j,k))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           r3 = r1*uur+(1d0-r1)*uul
           r2 = rand()
        enddo
        ptcl%z = r1*grd_zarr(k+1) + (1d0-r1) * &
             grd_zarr(k)

!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- setting IMC logical
        lhelp = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k)) * &
             min(dx(i),dy(j),dz(k))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
           z0 = ptcl%z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           ptcl%mu = (mu0+z0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           ptcl%om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(ptcl%om<0d0) ptcl%om = ptcl%om+pc_pi2
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%om = om0
           ptcl%itype = 2 !DDMC
        endif!}}}
!-- 1D spherical
     case(4)
!-- calculating position:!{{{
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,grd_nx)  !-- right neighbor
        uul = .5d0*(grd_emit(il,1,1) + grd_emit(i,1,1))
        uur = .5d0*(grd_emit(ir,1,1) + grd_emit(i,1,1))
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           x0 = (r1*grd_xarr(i+1)**3+(1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-grd_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
        enddo
        ptcl%x = x0
!-- setting IMC logical
        lhelp = ((grd_sig(i,1,1)+grd_cap(ig,i,1,1))*dx(i)* &
             thelp < prt_tauddmc).or.(in_puretran)
!write(0,*) i,grd_sig(i,1,1),grd_cap(ig,i,1,1),dx(i),thelp,prt_tauddmc

!-- if velocity-dependent, transforming direction
        if (lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
           ptcl%itype = 1 !IMC
        else
           ptcl%mu = mu0
           ptcl%itype = 2 !DDMC
        endif
!}}}
     endselect

!-- transformation into lab frame (in static grids cmffact==1d0)
     ptcl%e = ep0*cmffact
     ptcl%e0 = ep0*cmffact
     ptcl%wl = wl0/cmffact
!-- velocity effects accounting
     tot_evelo=tot_evelo+ep0*(1d0-cmffact)

!-- save particle result
!-----------------------
     prt_particles(ivac) = ptcl

!}}}
  enddo !ipart
!
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=prt_nnew) stop 'interior_source: n/=nnew'


end subroutine interior_source
