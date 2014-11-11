subroutine boundary_source

  use particlemod
  use timestepmod
  use physconstmod
  use gridmod
  use totalsmod
  use inputparmod
  use miscmod, only:specint
  implicit none

  logical :: lhelp
  integer :: ipart, ivac, ig, iig, i,j,k
  real*8 :: r1, r2, P, mu0, x0,y0,z0, esurfpart, wl0, om0
  real*8 :: denom2, wl1, wl2, thelp, mfphelp, mu1, mu2
  real*8 :: srftemp = 1d4
  real*8 :: cmffact,azitrfm
  type(packet),pointer :: ptcl
  integer, external :: binsrch
  real*8, dimension(grd_ng) :: emitsurfprobg  !surface emission probabilities 

!
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)

!
  tot_eleft = 0d0
  tot_eright = 0d0
!

  esurfpart = tot_esurf/dble(prt_nsurf)

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- calculating grouped thermal emission probabilities
  if(in_opacanaltype=='pick') then
     emitsurfprobg(1) = in_suolpick1
     emitsurfprobg(2) = 1d0 - in_suolpick1
     do ig = 3, grd_ng
        emitsurfprobg(3) = 0d0
     enddo
  else
     select case(in_igeom)
!-- 1D: sphere surface
     case(1)
        i = grd_nx
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
     do ig = 1, grd_ng
        wl1 = pc_h*pc_c/(grd_wl(ig+1)*pc_kb*srftemp)
        wl2 = pc_h*pc_c/(grd_wl(ig)*pc_kb*srftemp)
        emitsurfprobg(ig) = 15d0*specint(wl1,wl2,3)/pc_pi**4 
     enddo
  endif
  

!-- instantiating surface particles:
  do ipart = 1, prt_nsurf

!-- filling vacant spot in vacancy array
     ivac = prt_vacantarr(ipart)
     prt_isvacant(ivac) = .false.
     ptcl => prt_particles(ivac)
!
!-- calculating particle time
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     ptcl%tsrc = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grd_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     if(grd_isvelocity.and.in_srctype=='manu') then
        iig = 2
     endif
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/grd_wl(iig)+r1/grd_wl(iig+1))

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
        ptcl%rsrc = grd_xarr(i+1)
        x0 = ptcl%rsrc
!-- setting cell index
        ptcl%zsrc = i
!-- calculating albedo
        mfphelp = (grd_cap(iig,i,1,1)+grd_sig(i,1,1))*dx(i)*thelp
        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(i,1,1)+grd_cap(iig,i,1,1)) * &
             dx(i)*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0-mu0*x0/pc_c
!-- mu
           ptcl%musrc = (-mu0+x0/pc_c)/cmffact
        else
           ptcl%musrc = -mu0
        endif

!-- 2D (cylinder base)
     case(2)
!-- calculating position
        r1 = rand()
        ptcl%rsrc = sqrt(r1)*grd_xarr(grd_nx+1)
        x0 = ptcl%rsrc
        ptcl%y = grd_yarr(j)
        y0 = ptcl%y
!-- setting cell index
        ptcl%zsrc = binsrch(x0,grd_xarr,grd_nx+1,0)
        i = ptcl%zsrc
        ptcl%iy = j
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- calculating albedo
        mfphelp = (grd_cap(iig,i,j,1)+grd_sig(i,j,1))*dy(j)*thelp
        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(i,j,1)+grd_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%musrc = (mu0+y0/pc_c)/cmffact
           if(ptcl%musrc>1d0) then
              ptcl%musrc = 1d0
           elseif(ptcl%musrc<-1d0) then
              ptcl%musrc = -1d0
           endif
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

!-- 3D (box base)
     case(3)
!-- calculating position
        r1 = rand()
        ptcl%rsrc = r1*grd_xarr(grd_nx+1)
        x0 = ptcl%rsrc
        r1 = rand()
        ptcl%y = r1*grd_yarr(grd_ny+1)
        y0 = ptcl%y
        ptcl%z = grd_zarr(k)
        z0 = ptcl%z
!-- setting cell index
        ptcl%zsrc = binsrch(x0,grd_xarr,grd_nx+1,0)
        i = ptcl%zsrc
        ptcl%iy = binsrch(y0,grd_yarr,grd_ny+1,0)
        j = ptcl%iy
        ptcl%iz = k
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- calculating albedo
        mfphelp = (grd_cap(iig,i,j,k)+grd_sig(i,j,k))*dz(k)*thelp
        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(i,j,k)+grd_cap(iig,i,j,k)) * &
             min(dx(i),dy(j),dz(k))*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           ptcl%musrc = (mu0+z0/pc_c)/cmffact
           if(ptcl%musrc>1d0) then
              ptcl%musrc = 1d0
           elseif(ptcl%musrc<-1d0) then
              ptcl%musrc = -1d0
           endif
!-- om
           ptcl%om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(ptcl%om<0d0) ptcl%om = ptcl%om+pc_pi2
        else
           ptcl%musrc = mu0
           ptcl%om = om0
        endif
     endselect

     if(lhelp) then
!-- IMC
        if(grd_isvelocity) then
           tot_eext = tot_eext+esurfpart
           ptcl%esrc = esurfpart*cmffact
           ptcl%ebirth = esurfpart*cmffact
           ptcl%wlsrc = wl0/cmffact
!-- velocity effects accounting
           tot_evelo=tot_evelo-esurfpart*(cmffact-1d0)
        else
           ptcl%esrc = esurfpart
           ptcl%ebirth = esurfpart
           ptcl%wlsrc = wl0
        endif
        ptcl%rtsrc = 1
     else
!-- DDMC
        ptcl%esrc = P*esurfpart
        ptcl%ebirth = P*esurfpart
        tot_eext = tot_eext+ptcl%esrc
        ptcl%wlsrc = wl0
        ptcl%rtsrc = 2
     endif


  enddo


end subroutine boundary_source
