subroutine interior_source

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
  integer :: i,j,k,il,ir,ipart,ivac,ig,iig,ii
  integer :: nhere,ndmy,iimpi,nemit
  real*8 :: pwr
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, y0, z0, ep0, wl0
  real*8 :: denom2,x1,x2,x3,x4, thelp
  real*8 :: cmffact,azitrfm,mu1,mu2
  type(packet),pointer :: ptcl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- shortcut
  pwr = in_srcepwr

  x1=1d0/grd_wl(grd_ng+1)
  x2=1d0/grd_wl(1)

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
     ptcl => prt_particles(ivac)

!-- setting 1st cell index
     ptcl%ix = i

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grd_ng
        x3=1d0/grd_wl(ig+1)
        x4=1d0/grd_wl(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
        denom2 = denom2+(x4-x3)/(x2-x1)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/grd_wl(iig)+r1/grd_wl(iig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = grd_emitex(i,j,k)/dble(grd_nvol(i,j,k)-nemit)
     tot_eext=tot_eext+ep0

!
!-- selecting geometry
     select case(in_igeom)

!-- 1D
     case(1)
!-- calculating position
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1.0-r1)*grd_xarr(i)**3)**(1.0/3.0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- setting IMC logical
        lhelp = ((grd_sig(i,1,1)+grd_cap(iig,i,1,1))*dx(i)* &
             thelp < prt_tauddmc).or.(in_puretran)

!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
        else
           ptcl%mu = mu0
        endif

!-- 2D
     case(2)
!-- setting 2nd cell index
        ptcl%iy = j
!-- calculating position
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
        lhelp = ((grd_sig(i,j,1)+grd_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%mu = (mu0+y0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           if(azitrfm >= 0d0) then
              ptcl%om = azitrfm
           else
              ptcl%om = azitrfm+pc_pi2
           endif
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif

!-- 3D
     case(3)
!-- setting 2nd,3rd cell index
        ptcl%iy = j
        ptcl%iz = k
!-- calculating position
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
        lhelp = ((grd_sig(i,j,k)+grd_cap(iig,i,j,k)) * &
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
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif

     endselect

     if (lhelp) then
!-- IMC
        if(grd_isvelocity) then
           ptcl%e = ep0*cmffact
           ptcl%e0 = ep0*cmffact
           ptcl%wl = wl0/cmffact
!-- velocity effects accounting
           tot_evelo=tot_evelo+ep0*(1d0-cmffact)
        else
           ptcl%e = ep0
           ptcl%e0 = ep0
           ptcl%wl = wl0
        endif
        ptcl%itype = 1
     else
!-- DDMC
        ptcl%e = ep0
        ptcl%e0 = ep0
        ptcl%wl = wl0
        ptcl%itype = 2
     endif
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
  do ii=1,nhere
     ipart = ipart + 1!{{{
     ivac = prt_vacantarr(ipart)
     ptcl => prt_particles(ivac)
!
!-- setting 1st cell index
     ptcl%ix = i

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1     
     do ig = 1, grd_ng
        iig = ig
        if (r1>=denom2.and.r1<denom2+grd_emitprob(ig,i,j,k)) exit
        denom2 = denom2+grd_emitprob(ig,i,j,k)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/grd_wl(iig)+r1/grd_wl(iig+1))

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
        lhelp = ((grd_sig(i,1,1)+grd_cap(iig,i,1,1))*dx(i)* &
             thelp < prt_tauddmc).or.(in_puretran)
!write(0,*) i,grd_sig(i,1,1),grd_cap(iig,i,1,1),dx(i),thelp,prt_tauddmc

!-- if velocity-dependent, transforming direction
        if (lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
        else
           ptcl%mu = mu0
        endif
!}}}
!-- 2D
     case(2)
!-- setting 2nd cell index!{{{
        ptcl%iy = j
!-- calculating position:
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
        lhelp = ((grd_sig(i,j,1)+grd_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%mu = (mu0+y0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           if(azitrfm >= 0d0) then
              ptcl%om = azitrfm
           else
              ptcl%om = azitrfm+pc_pi2
           endif
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif
!}}}
!-- 3D
     case(3)
!-- setting 2nd,3rd cell index!{{{
        ptcl%iy = j
        ptcl%iz = k
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
        lhelp = ((grd_sig(i,j,k)+grd_cap(iig,i,j,k)) * &
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
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif!}}}
     endselect


     if (lhelp) then
!-- IMC
        if(grd_isvelocity) then
           ptcl%e = ep0*cmffact
           ptcl%e0 = ep0*cmffact
           ptcl%wl = wl0/cmffact
!-- velocity effects accounting
           tot_evelo=tot_evelo+ep0*(1d0-cmffact)
        else
           ptcl%e = ep0
           ptcl%e0 = ep0
           ptcl%wl = wl0
        endif
        ptcl%itype = 1
     else
!-- DDMC
        ptcl%e = ep0
        ptcl%e0 = ep0
        ptcl%wl = wl0
        ptcl%itype = 2
     endif

!}}}
  enddo !ipart
  enddo !i
  enddo !j
  enddo !k
  if(ipart/=prt_nnew) stop 'interior_source: n/=nnew'


end subroutine interior_source
