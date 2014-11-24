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
  real*8 :: denom2, wl1, wl2, thelp, mfphelp, help, mu1, mu2
  real*8 :: srftemp = 1d4
  real*8 :: cmffact, alb,beta,eps
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
!-- 1D
     case(1)
        i = grd_nx
        j = 1
        k = 1
!-- 2D
     case(2)
        if(in_surfsrcloc=='down') then
           j = 1
        elseif(in_surfsrcloc=='up') then
           j = grd_ny
        else
           i = grd_nx
        endif
        k = 1
!-- 3D
     case(3)
        if(in_surfsrcloc=='in') then
           i = 1
        elseif(in_surfsrcloc=='out') then
           i = grd_nx
        elseif(in_surfsrcloc=='down') then
           j = 1
        elseif(in_surfsrcloc=='up') then
           j = grd_ny
        elseif(in_surfsrcloc=='botm') then
           k = 1
        else
           k = grd_nz
        endif
     endselect
     if(in_srcmax>0d0.and.in_srctype=='surf') srftemp = in_srcmax
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
     ptcl%t = tsp_t+r1*tsp_dt

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
     if(in_surfsrcmu=='beam') then
        mu0 = 1d0
     elseif(in_surfsrcmu=='isot') then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        r2 = rand()
        prt_tlyrand = prt_tlyrand+1
        mu0 = max(r1,r2)
     endif
!
!-- selecting geometry and surface
     select case(in_igeom)

!-- 1D (outer surface)
     case(1)
!-- calculating position
        ptcl%x = grd_xarr(i+1)
        x0 = ptcl%x
!-- setting cell index
        ptcl%ix = i
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
           ptcl%mu = (-mu0+x0/pc_c)/cmffact
        else
           ptcl%mu = -mu0
        endif

!-- 2D
     case(2)
!-- calculating position
        if(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- flat surface
           r1 = rand()
           ptcl%x = sqrt(r1)*grd_xarr(grd_nx+1)
           x0 = ptcl%x
           if(j==1) then
              ptcl%y = grd_yarr(j)
           else
              ptcl%y = grd_yarr(j+1)
              mu0 = -mu0
           endif
           y0 = ptcl%y
!-- setting cell index
           ptcl%ix = binsrch(x0,grd_xarr,grd_nx+1,0)
           i = ptcl%ix
           ptcl%iy = j
!-- sampling direction helpers
           r1 = rand()
           om0 = pc_pi2*r1
!-- setting albedo helpers
           help = dy(j)
           mu2 = mu0
        else
!-- curved surface
           ptcl%x = grd_xarr(grd_nx+1)
           x0 = ptcl%x
           r1 = rand()
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
!-- setting cell index
           ptcl%ix = i
           ptcl%iy = binsrch(y0,grd_yarr,grd_ny+1,0)
           j = ptcl%iy
!-- sampling direction helpers
           r1 = rand()
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),-mu0)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dx(i)
           mu2 = mu0
           mu0 = mu1
        endif

!-- calculating albedo
        mfphelp = (grd_cap(iig,i,j,1)+grd_sig(i,j,1))*help*thelp
        P = 4d0*(1.0+1.5*abs(mu2))/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(i,j,1)+grd_cap(iig,i,j,1)) * &
             min(dx(i),dy(j))*thelp < prt_tauddmc) &
             .or.in_puretran.or.P>1d0.or.P<0d0
!-- if velocity-dependent, transforming direction
        if(lhelp.and.grd_isvelocity) then
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           ptcl%om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           ptcl%mu = (mu0+y0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           if(ptcl%om < 0d0) ptcl%om=ptcl%om+pc_pi2
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif

!-- 3D
     case(3)
!-- calculating position
        if(in_surfsrcloc=='in'.or.in_surfsrcloc=='out') then
!-- x surface
           if(i==1) then
              ptcl%x = grd_xarr(i)
           else
              ptcl%x = grd_xarr(i+1)
              mu0 = -mu0
           endif
           x0 = ptcl%x
           r1 = rand()
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
           r1 = rand()
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           ptcl%ix = i
           ptcl%iy = binsrch(y0,grd_yarr,grd_ny+1,0)
           j = ptcl%iy
           ptcl%iz = binsrch(z0,grd_zarr,grd_nz+1,0)
           k = ptcl%iz
!-- sampling direction helpers
           r1 = rand()
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu1,mu0)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dx(i)
           mu2 = mu0
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        elseif(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- y surface
           r1 = rand()
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           if(j==1) then
              ptcl%y = grd_yarr(j)
           else
              ptcl%y = grd_yarr(j+1)
              mu0 = -mu0
           endif
           y0 = ptcl%y
           r1 = rand()
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           ptcl%ix = binsrch(x0,grd_xarr,grd_nx+1,0)
           i = ptcl%ix
           ptcl%iy = j
           ptcl%iz = binsrch(z0,grd_zarr,grd_nz+1,0)
           k = ptcl%iz
!-- sampling direction helpers
           r1 = rand()
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu0,mu1)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dy(j)
           mu2 = mu0
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        else
!-- z surface
           r1 = rand()
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           r1 = rand()
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
           if(k==1) then
              ptcl%z = grd_zarr(k)
           else
              ptcl%z = grd_zarr(k+1)
              mu0 = -mu0
           endif
           z0 = ptcl%z
!-- setting cell index
           ptcl%ix = binsrch(x0,grd_xarr,grd_nx+1,0)
           i = ptcl%ix
           ptcl%iy = binsrch(y0,grd_yarr,grd_ny+1,0)
           j = ptcl%iy
           ptcl%iz = k
!-- sampling azimuthal angle of direction
           r1 = rand()
           om0 = pc_pi2*r1
!-- setting albedo helpers
           help = dz(k)
           mu2 = mu0
        endif

!-- calculating albedo
        mfphelp = (grd_cap(iig,i,j,k)+grd_sig(i,j,k))*help*thelp
        alb = grd_fcoef(i,j,k)*grd_cap(iig,i,j,k)/ &
             (grd_cap(iig,i,j,k)+grd_sig(i,j,k))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        P = 0.5d0*eps*beta*(1d0+1.5d0*abs(mu2))/(beta-0.75*eps*mfphelp)
!        P = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
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

     if(lhelp) then
!-- IMC
        if(grd_isvelocity) then
           tot_eext = tot_eext+esurfpart
           ptcl%e = esurfpart*cmffact
           ptcl%e0 = esurfpart*cmffact
           ptcl%wl = wl0/cmffact
!-- velocity effects accounting
           tot_evelo=tot_evelo-esurfpart*(cmffact-1d0)
        else
           ptcl%e = esurfpart
           ptcl%e0 = esurfpart
           ptcl%wl = wl0
        endif
        ptcl%itype = 1
     else
!-- DDMC
        ptcl%e = P*esurfpart
        ptcl%e0 = P*esurfpart
        tot_eext = tot_eext+ptcl%e
        ptcl%wl = wl0
        ptcl%itype = 2
     endif


  enddo


end subroutine boundary_source
