!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine boundary_source

  use randommod
  use particlemod
  use transportmod
  use sourcemod
  use timestepmod
  use physconstmod
  use groupmod
  use gridmod
  use totalsmod
  use inputparmod
  use miscmod
  implicit none

  logical :: lhelp
  integer :: ipart, ivac, ig, iig, i,j,k
  real*8 :: r1, r2, p, mu0, x0,y0,z0, esurfpart, wl0, om0
  real*8 :: denom2, wl1, wl2, thelp, mfphelp, help, mu1, mu2
  real*8 :: srftemp = 1d4
  real*8 :: alb,beta,eps
  type(packet) :: ptcl
  real*8, dimension(grp_ng) :: emitsurfprobg  !surface emission probabilities 

!
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz,xm,dyac,ym
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)

  esurfpart = tot_esurf/dble(src_nsurftot)

  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- calculating grouped thermal emission probabilities
  if(in_opacanaltype=='pick') then
     emitsurfprobg(1) = in_suolpick1
     emitsurfprobg(2) = 1d0 - in_suolpick1
     do ig = 3, grp_ng
        emitsurfprobg(3) = 0d0
     enddo
  else
     select case(grd_igeom)
!-- 3D spherical
     case(1)
        i = grd_nx
!-- 2D
     case(2)
        if(in_surfsrcloc=='down') then
           j = 1
        elseif(in_surfsrcloc=='up') then
           j = grd_ny
        else
           i = grd_nx
        endif
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
!-- 1D spherical
     case(11)
        i = grd_nx
        j = 1
        k = 1
        l = grd_icell(i,j,k)
     endselect
     if(in_srcmax>0d0.and.in_srctype=='surf') srftemp = in_srcmax
     do ig = 1, grp_ng
        wl1 = pc_h*pc_c*grp_wlinv(ig+1)/(pc_kb*srftemp)
        wl2 = pc_h*pc_c*grp_wlinv(ig)/(pc_kb*srftemp)
        emitsurfprobg(ig) = 15d0*specint(wl1,wl2,3)/pc_pi**4 
     enddo
  endif
  

!-- instantiating surface particles:
  do ipart=1,src_nsurf

!-- filling vacant spot in vacancy array
     ivac = src_ivacant(ipart)
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     call rnd_r(r1,rnd_state)
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     call rnd_r(r1,rnd_state)
     do ig = 1, grp_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     if(grd_isvelocity.and.in_srctype=='manu') then
        iig = 2
     endif
     call rnd_r(r1,rnd_state)
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))

!-- sampling surface projection
     if(in_surfsrcmu=='beam') then
        mu0 = 1d0
     elseif(in_surfsrcmu=='isot') then
        call rnd_r(r1,rnd_state)
        call rnd_r(r2,rnd_state)
        mu0 = max(r1,r2)
     endif
!
!-- selecting geometry and surface
     select case(grd_igeom)

!-- 3D (outer surface)
     case(1)
!-- calculating position!{{{
        ptcl%x = grd_xarr(i+1)
        x0 = ptcl%x
        call rnd_r(r1,rnd_state)
        ptcl%y = 1d0-2d0*r1
        call rnd_r(r1,rnd_state)
        ptcl%z = pc_pi2*r1
!-- cell index
        j = binsrch(ptcl%y,grd_yarr,grd_ny+1,.false.)
        k = binsrch(ptcl%z,grd_zarr,grd_nz+1,.false.)
!-- sampling azimuthal direction angle
        ptcl%om = r1*pc_pi2
!-- calculating albedo
        l = grd_icell(i,j,k)
        mfphelp = (grd_cap(iig,l)+grd_sig(l))*dx(i)*thelp
        p = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(l)+grd_cap(iig,l)) * &
             min(dx(i),xm(i)*dyac(j),xm(i)*ym(j)*dz(k)) * &
             thelp < trn_tauddmc) &
             .or.in_puretran.or.p>1d0.or.p<0d0
        ptcl%mu = -mu0!}}}

!-- 2D
     case(2)
!-- calculating position!{{{
        if(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- flat surface
           call rnd_r(r1,rnd_state)
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
           i = binsrch(x0,grd_xarr,grd_nx+1,.false.)
!-- sampling direction helpers
           call rnd_r(r1,rnd_state)
           om0 = pc_pi2*r1
!-- setting albedo helpers
           help = dy(j)
           mu2 = mu0
        else
!-- curved surface
           ptcl%x = grd_xarr(grd_nx+1)
           x0 = ptcl%x
           call rnd_r(r1,rnd_state)
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
!-- setting cell index
           j = binsrch(y0,grd_yarr,grd_ny+1,.false.)
!-- sampling direction helpers
           call rnd_r(r1,rnd_state)
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),-mu0)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dx(i)
           mu2 = mu0
           mu0 = mu1
        endif
        call rnd_r(r1,rnd_state)
        ptcl%z = pc_pi2*r1
        k = binsrch(ptcl%z,grd_zarr,grd_nz+1,.false.)

!-- calculating albedo
        l = grd_icell(i,j,k)
        mfphelp = (grd_cap(iig,l)+grd_sig(l))*help*thelp
        p = 4d0*(1.0+1.5*abs(mu2))/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(l)+grd_cap(iig,l)) * &
             min(dx(i),dy(j),xm(i)*dz(k))*thelp < trn_tauddmc) &
             .or.in_puretran.or.p>1d0.or.p<0d0

        ptcl%mu = mu0
        ptcl%om = om0!}}}

!-- 3D
     case(3)
!-- calculating position!{{{
        if(in_surfsrcloc=='in'.or.in_surfsrcloc=='out') then
!-- x surface
           if(i==1) then
              ptcl%x = grd_xarr(i)
           else
              ptcl%x = grd_xarr(i+1)
              mu0 = -mu0
           endif
           x0 = ptcl%x
           call rnd_r(r1,rnd_state)
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
           call rnd_r(r1,rnd_state)
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           j = binsrch(y0,grd_yarr,grd_ny+1,.false.)
           k = binsrch(z0,grd_zarr,grd_nz+1,.false.)
!-- sampling direction helpers
           call rnd_r(r1,rnd_state)
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu1,mu0)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dx(i)
           mu2 = mu0
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        elseif(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- y surface
           call rnd_r(r1,rnd_state)
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           if(j==1) then
              ptcl%y = grd_yarr(j)
           else
              ptcl%y = grd_yarr(j+1)
              mu0 = -mu0
           endif
           y0 = ptcl%y
           call rnd_r(r1,rnd_state)
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           i = binsrch(x0,grd_xarr,grd_nx+1,.false.)
           k = binsrch(z0,grd_zarr,grd_nz+1,.false.)
!-- sampling direction helpers
           call rnd_r(r1,rnd_state)
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu0,mu1)
           if(om0<0d0) om0=om0+pc_pi2
!-- setting albedo helpers
           help = dy(j)
           mu2 = mu0
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        else
!-- z surface
           call rnd_r(r1,rnd_state)
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           call rnd_r(r1,rnd_state)
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
           i = binsrch(x0,grd_xarr,grd_nx+1,.false.)
           j = binsrch(y0,grd_yarr,grd_ny+1,.false.)
!-- sampling azimuthal angle of direction
           call rnd_r(r1,rnd_state)
           om0 = pc_pi2*r1
!-- setting albedo helpers
           help = dz(k)
           mu2 = mu0
        endif

!-- calculating albedo
        l = grd_icell(i,j,k)
        mfphelp = (grd_cap(iig,l)+grd_sig(l))*help*thelp
        alb = grd_fcoef(l)*grd_cap(iig,l)/ &
             (grd_cap(iig,l)+grd_sig(l))
        eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
        beta = 1.5d0*alb*mfphelp**2+sqrt(3d0*alb*mfphelp**2 + &
             2.25d0*alb**2*mfphelp**4)
        p = 0.5d0*eps*beta*(1d0+1.5d0*abs(mu2))/(beta-0.75*eps*mfphelp)
!        p = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(l)+grd_cap(iig,l)) * &
             min(dx(i),dy(j),dz(k))*thelp < trn_tauddmc) &
             .or.in_puretran.or.p>1d0.or.p<0d0

        ptcl%mu = mu0
        ptcl%om = om0!}}}

!-- 1D (outer surface)
     case(11)
!-- calculating position!{{{
        ptcl%x = grd_xarr(i+1)
        x0 = ptcl%x
        ptcl%y = grd_yarr(1)
        ptcl%z = grd_zarr(1)
!-- calculating albedo
        mfphelp = (grd_cap(iig,l)+grd_sig(l))*dx(i)*thelp
        p = 4d0*(1.0+1.5*mu0)/(3d0*mfphelp+6d0*pc_dext)
        lhelp = ((grd_sig(l)+grd_cap(iig,l)) * &
             dx(i)*thelp < trn_tauddmc) &
             .or.in_puretran.or.p>1d0.or.p<0d0

        ptcl%mu = -mu0!}}}
     endselect

!-- particle properties are saved in comoving frame
     ptcl%e = esurfpart
     ptcl%e0 = esurfpart
     ptcl%wl = wl0

!-- DDMC
     if(.not.lhelp) then
        ptcl%e = ptcl%e*p
        ptcl%e0 = ptcl%e0*p
     endif

!-- save particle result
!-----------------------
     prt_particles(ivac) = ptcl


  enddo


end subroutine boundary_source
! vim: fdm=marker
