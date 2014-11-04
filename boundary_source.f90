subroutine boundary_source

  use particlemod
  use timestepmod
  use physconstmod
  use gasgridmod
  use inputparmod
  use miscmod, only:specint
  implicit none

  logical :: lhelp
  integer :: ipart, ivac, ig, iig, i,j,k
  real*8 :: r1, r2, P, mu0, x0, y0, esurfpart, wl0, om0
  real*8 :: denom2, wl1, wl2, thelp, mfphelp
  real*8 :: srftemp = 1d4
  real*8 :: cmffact,azitrfm
  integer, external :: binsrch
  real*8, dimension(gas_ng) :: emitsurfprobg  !surface emission probabilities 

!
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
  dz(l) = gas_zarr(l+1) - gas_zarr(l)

!
  gas_eleft = 0d0
  gas_eright = 0d0
!

  esurfpart = gas_esurf/dble(prt_nsurf)

  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- calculating grouped thermal emission probabilities
  if(gas_opacanaltype=='pick') then
     emitsurfprobg(1) = gas_ppick(1)
     emitsurfprobg(2) = gas_ppick(2)
     do ig = 3, gas_ng
        emitsurfprobg(3) = 0d0
     enddo
  else
     select case(in_igeom)
!-- 1D: sphere surface
     case(1)
        i = gas_nx
        j = 1
        k = 1
!-- 2D: cylinder base
     case(2)
        j = 1
        k = 1
!-- 3D: box base
     case(3)
        k = 1
     endselect
     do ig = 1, gas_ng
        wl1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*srftemp)
        wl2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*srftemp)
        emitsurfprobg(ig) = 15d0*specint(wl1,wl2,3)/pc_pi**4 
     enddo
  endif
  

!-- instantiating surface particles:
  do ipart = 1, prt_nsurf

!-- filling vacant spot in vacancy array
     ivac = prt_vacantarr(ipart)
     prt_isvacant(ivac) = .false.

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, gas_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     if(gas_isvelocity.and.gas_srctype=='manu') then
        iig = 2
     endif
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = (1d0-r1)*gas_wl(iig)+r1*gas_wl(iig+1)

!-- sampling surface projection
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     r2 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = max(r1,r2)

!
!-- selecting geometry and surface
     select case(in_igeom)

!-- 1D (outer surface)
     case(1)
!-- calculating position
        prt_particles(ivac)%rsrc = gas_xarr(i+1)
        x0 = prt_particles(ivac)%rsrc
!-- setting cell index
        prt_particles(ivac)%zsrc = i
!-- calculating albedo
        mfphelp = (gas_cap(iig,i,1,1)+gas_sig(i,1,1))*dx(i)*thelp
        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((gas_sig(i,1,1)+gas_cap(iig,i,1,1)) * &
             dx(i)*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.gas_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0-mu0*x0/pc_c
!-- mu
           prt_particles(ivac)%musrc = (-mu0+x0/pc_c)/cmffact
        else
           prt_particles(ipart)%musrc = -mu0
        endif

!-- 2D (cylinder base)
     case(2)
!-- calculating position
        r1 = rand()
        prt_particles(ivac)%rsrc = sqrt(r1)*gas_xarr(gas_nx+1)
        x0 = prt_particles(ivac)%rsrc
        prt_particles(ivac)%y = gas_yarr(j)
        y0 = prt_particles(ivac)%y
!-- setting cell index
        prt_particles(ivac)%zsrc = binsrch(x0,gas_xarr,gas_nx+1,0)
        i = prt_particles(ivac)%zsrc
        prt_particles(ivac)%iy = j
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- calculating albedo
        mfphelp = (gas_cap(iig,i,j,1)+gas_sig(i,j,1))*dy(j)*thelp
        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((gas_sig(i,j,1)+gas_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.gas_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           prt_particles(ivac)%musrc = (mu0+y0/pc_c)/cmffact
           if(prt_particles(ivac)%musrc>1d0) then
              prt_particles(ivac)%musrc = 1d0
           elseif(prt_particles(ivac)%musrc<-1d0) then
              prt_particles(ivac)%musrc = -1d0
           endif
!-- om
           if(azitrfm >= 0d0) then
              prt_particles(ivac)%om = azitrfm
           else
              prt_particles(ivac)%om = azitrfm+pc_pi2
           endif
        else
           prt_particles(ivac)%musrc = mu0
           prt_particles(ivac)%om = om0
        endif

!-- 3D (box base)
     case(3)
        stop 'boundary_source: no 3D transport'
     endselect

     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     prt_particles(ivac)%tsrc = tsp_t+r1*tsp_dt


     if(lhelp) then
!-- IMC
        if(gas_isvelocity) then
           gas_eext = gas_eext+esurfpart
           prt_particles(ivac)%esrc = esurfpart*cmffact
           prt_particles(ivac)%ebirth = esurfpart*cmffact
           prt_particles(ivac)%wlsrc = wl0/cmffact
!-- velocity effects accounting
           gas_evelo=gas_evelo-esurfpart*(cmffact-1d0)
        else
           prt_particles(ivac)%esrc = esurfpart
           prt_particles(ivac)%ebirth = esurfpart
           prt_particles(ivac)%wlsrc = wl0
        endif
        prt_particles(ivac)%rtsrc = 1
     else
!-- DDMC
        prt_particles(ivac)%esrc = P*esurfpart
        prt_particles(ivac)%ebirth = P*esurfpart
        gas_eext = gas_eext+prt_particles(ivac)%esrc
        prt_particles(ivac)%wlsrc = wl0
        prt_particles(ivac)%rtsrc = 2
     endif


  enddo


end subroutine boundary_source
