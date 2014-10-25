subroutine interior_source

  use gasgridmod
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
  integer :: i,j,k,il,ir,ipart,ivac,ig,iig
  integer, dimension(gas_nx,gas_ny,gas_nz) :: ijkused
  real*8 :: r1, r2, r3, uul, uur, uumax
  real*8 :: om0, mu0, x0, y0, z0, ep0, wl0
  real*8 :: denom2,x1,x2,x3,x4, help
  real*8 :: cmffact,azitrfm
  type(packet),pointer :: ptcl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)

  if(gas_isvelocity) then
     help = tsp_t
  else
     help = 1d0
  endif

  i = 1
  j = 1
  k = 1
  ijkused = 0
  !Volume particle instantiation: loop
  !Loop run over the number of new particles that aren't surface source
  !particles.
  
  x1=1d0/gas_wl(gas_ng+1)
  x2=1d0/gas_wl(1)
  do ipart = prt_nsurf+1, prt_nsurf+prt_nexsrc
     ivac = prt_vacantarr(ipart)!{{{
     ptcl => prt_particles(ivac)
!
!-- check for available particle space to populate in cell
     do k=k,gas_nz
        do j=j,gas_ny
           do i=i,gas_nx
              lhelp = ijkused(i,j,k)<gas_nvolex(i,j,k)
              if (lhelp) exit
           enddo
           if (lhelp) exit
        enddo
        if (lhelp) exit
     enddo
!-- increasing cell occupancy
     ijkused(i,j,k) = ijkused(i,j,k)+1

!-- setting 1st cell index
     ptcl%zsrc = i

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%tsrc = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, gas_ng
        x3=1d0/gas_wl(ig+1)
        x4=1d0/gas_wl(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
        denom2 = denom2+(x4-x3)/(x2-x1)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = gas_emitex(i,j,k)/real(gas_nvolex(i,j,k))
     gas_eext=gas_eext+ep0

!
!-- selecting geometry
     select case(in_igeom)

!-- 1D
     case(1)
!-- calculating position
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        ptcl%rsrc = (r1*gas_xarr(i+1)**3 + &
             (1.0-r1)*gas_xarr(i)**3)**(1.0/3.0)
!-- setting IMC logical
        lhelp = ((gas_sig(i,1,1)+gas_cap(iig,i,1,1))*dx(i)* &
             help < prt_tauddmc).or.(in_puretran)

!-- if velocity-dependent, transforming direction
        if(lhelp.and.gas_isvelocity) then
           x0 = ptcl%rsrc
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%musrc = (mu0+x0/pc_c)/cmffact
        else
           ptcl%musrc = mu0
        endif

!-- 2D
     case(2)
!-- setting 2nd cell index
        ptcl%iy = j
!-- calculating position
        r1 = rand()
        ptcl%rsrc = sqrt(r1*gas_xarr(i+1)**2 + &
             (1d0-r1)*gas_xarr(i)**2)
        r1 = rand()
        ptcl%y = r1*gas_yarr(j+1) + (1d0-r1) * &
             gas_yarr(j)
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- setting IMC logical
        lhelp = ((gas_sig(i,j,1)+gas_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*help < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.gas_isvelocity) then
           x0 = ptcl%rsrc
           y0 = ptcl%y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%musrc = (mu0+y0/pc_c)/cmffact
!-- om
           if(azitrfm >= 0d0) then
              ptcl%om = azitrfm
           else
              ptcl%om = azitrfm+pc_pi2
           endif
        else
           ptcl%musrc = mu0
           ptcl%om = om0
        endif

!-- 3D
     case(3)
        stop 'interior_source: no 3D transport'
     endselect

     if (lhelp) then
!-- IMC
        if(gas_isvelocity) then
           ptcl%esrc = ep0*cmffact
           ptcl%ebirth = ep0*cmffact
           ptcl%wlsrc = wl0/cmffact
!-- velocity effects accounting
           gas_evelo=gas_evelo-ep0*x0*mu0/pc_c
        else
           ptcl%esrc = ep0
           ptcl%ebirth = ep0
           ptcl%wlsrc = wl0
        endif
        ptcl%rtsrc = 1
     else
!-- DDMC
        ptcl%esrc = ep0
        ptcl%ebirth = ep0
        ptcl%wlsrc = wl0
        ptcl%rtsrc = 2
     endif
!
  enddo
  

!-- Thermal volume particle instantiation: loop
  i = 1
  j = 1
  k = 1
  ijkused = 0

  do ipart = prt_nsurf+prt_nexsrc+1, prt_nnew
     ivac = prt_vacantarr(ipart)!{{{
     ptcl => prt_particles(ivac)

!
!-- check for available particle space to populate in cell
     do k=k,gas_nz
        do j=j,gas_ny
           do i=i,gas_nx
              lhelp = ijkused(i,j,k)<gas_nvol(i,j,k)
              if (lhelp) exit
           enddo
           if (lhelp) exit
        enddo
        if (lhelp) exit
     enddo
!-- increasing cell occupancy
     ijkused(i,j,k) = ijkused(i,j,k)+1
!
!-- setting 1st cell index
     ptcl%zsrc = i

!-- setting particle index to not vacant
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%tsrc = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1     
     do ig = 1, gas_ng
        iig = ig
        if (r1>=denom2.and.r1<denom2+gas_emitprob(ig,i,j,k)) exit
        denom2 = denom2+gas_emitprob(ig,i,j,k)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = gas_emit(i,j,k)/real(gas_nvol(i,j,k))

!
!-- selecting geometry
     select case(in_igeom)

!-- 1D
     case(1)
!-- calculating position:
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,gas_nx)  !-- right neighbor
        uul = .5d0*(gas_temp(il,1,1)**4 + gas_temp(i,1,1)**4)
        uur = .5d0*(gas_temp(ir,1,1)**4 + gas_temp(i,1,1)**4)
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           x0 = (r1*gas_xarr(i+1)**3+(1.0-r1)*gas_xarr(i)**3)**(1.0/3.0)
           r3 = (x0-gas_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
        enddo
        ptcl%rsrc = x0
!-- setting IMC logical
        lhelp = ((gas_sig(i,1,1)+gas_cap(iig,i,1,1))*dx(i)* &
             help < prt_tauddmc).or.(in_puretran)

!-- if velocity-dependent, transforming direction
        if (lhelp.and.gas_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+mu0*x0/pc_c
!-- mu
           ptcl%musrc = (mu0+x0/pc_c)/cmffact
        else
           ptcl%musrc = mu0
        endif

!-- 2D
     case(2)
!-- setting 2nd cell index
        ptcl%iy = j
!-- calculating position:
!-- source tilting in x
        r3 = 0d0
        r2 = 1d0
        il = max(i-1,1)  !-- left neighbor
        ir = min(i+1,gas_nx)  !-- right neighbor
        uul = .5d0*(gas_temp(il,j,1)**4 + gas_temp(i,j,1)**4)
        uur = .5d0*(gas_temp(ir,j,1)**4 + gas_temp(i,j,1)**4)
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           x0 = sqrt(r1*gas_xarr(i+1)**2+(1.0-r1)*gas_xarr(i)**2)
           r3 = (x0-gas_xarr(i))/dx(i)
           r3 = r3*uur+(1.0-r3)*uul
           r2 = rand()
        enddo
        ptcl%rsrc = x0
!- source tilting in y
        r3 = 0d0
        r2 = 1d0
        il = max(j-1,1)  !-- lower neighbor
        ir = min(j+1,gas_ny)  !-- upper neighbor
        uul = .5d0*(gas_temp(i,il,1)**4 + gas_temp(i,j,1)**4)
        uur = .5d0*(gas_temp(i,ir,1)**4 + gas_temp(i,j,1)**4)
        uumax = max(uul,uur)
        uul = uul/uumax
        uur = uur/uumax
        do while (r2 > r3)
           r1 = rand()
           r3 = r1*uur+(1d0-r1)*uul
           r2 = rand()
        enddo
        y0 = r1*gas_yarr(j+1)+(1d0-r1)*gas_yarr(j)
        ptcl%y = y0
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1

!-- setting IMC logical
        lhelp = ((gas_sig(i,j,1)+gas_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*help < prt_tauddmc) &
             .or.in_puretran
!-- if velocity-dependent, transforming direction
        if(lhelp.and.gas_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%musrc = (mu0+y0/pc_c)/cmffact
!-- om
           if(azitrfm >= 0d0) then
              ptcl%om = azitrfm
           else
              ptcl%om = azitrfm+pc_pi2
           endif
        else
           ptcl%musrc = mu0
           ptcl%om = om0
        endif

!-- 3D
     case(3)
        stop 'interior_source: no 3D transport'
     endselect

     if (lhelp) then
!-- IMC
        if(gas_isvelocity) then
           ptcl%esrc = ep0*cmffact
           ptcl%ebirth = ep0*cmffact
           ptcl%wlsrc = wl0/cmffact
!-- velocity effects accounting
           gas_evelo=gas_evelo-ep0*x0*mu0/pc_c
        else
           ptcl%esrc = ep0
           ptcl%ebirth = ep0
           ptcl%wlsrc = wl0
        endif
        ptcl%rtsrc = 1
     else
!-- DDMC
        ptcl%esrc = ep0
        ptcl%ebirth = ep0
        ptcl%wlsrc = wl0
        ptcl%rtsrc = 2
     endif

!}}}
  enddo


end subroutine interior_source
