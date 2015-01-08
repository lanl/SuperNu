subroutine boundary_source

  use randommod
  use particlemod
  use timestepmod
  use physconstmod
  use groupmod
  use gridmod
  use totalsmod
  use inputparmod
  use miscmod
  implicit none

  logical :: lhelp
  integer :: ipart, ivac, ig, iig, i,j,k, ix,iy,iz, l
  real*8 :: r1, r2, P, mu0, x0,y0,z0, esurfpart, wl0, om0
  real*8 :: denom2, wl1, wl2, thelp, mfphelp, help, mu1, mu2
  real*8 :: srftemp = 1d4
  real*8 :: cmffact, alb,beta,eps,gm
  type(packet) :: ptcl
  real*8, dimension(grp_ng) :: emitsurfprobg  !surface emission probabilities 


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
     do ig = 3, grp_ng
        emitsurfprobg(3) = 0d0
     enddo
  else
     select case(in_igeom)
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
  do ipart=1,prt_nsurf

!-- filling vacant spot in vacancy array
     ivac = prt_vacantarr(ipart)
     prt_isvacant(ivac) = .false.
!
!-- calculating particle time
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     ptcl%t = tsp_t+r1*tsp_dt

!-- calculating wavelength
     denom2 = 0d0
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grp_ng
        iig = ig
        if(r1>=denom2.and.r1<denom2+emitsurfprobg(ig)) exit
        denom2 = denom2+emitsurfprobg(ig)
     enddo
     if(grd_isvelocity.and.in_srctype=='manu') then
        iig = 2
     endif
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))

!-- sampling surface projection
     if(in_surfsrcmu=='beam') then
        mu0 = 1d0
     elseif(in_surfsrcmu=='isot') then
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        r2 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        mu0 = max(r1,r2)
     endif
!
!-- selecting geometry and surface
     select case(in_igeom)

!-- 3D (outer surface)
     case(1)
!-- calculating position!{{{
        ptcl%x = grd_xarr(i+1)
        x0 = ptcl%x
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = 1d0-2d0*r1
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = pc_pi2*r1
        j = binsrch(ptcl%y,grd_yarr,grd_ny+1,.false.)
        k = binsrch(ptcl%z,grd_zarr,grd_nz+1,.false.)
!-- setting cell index
        ix = i
        iy = j
        iz = k
!-- sampling azimuthal direction angle
        ptcl%om = r1*pc_pi2
        ptcl%mu = -mu0!}}}

!-- 2D
     case(2)
!-- calculating position!{{{
        if(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- flat surface
           r1 = rnd_r(rnd_state)
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
           ix = binsrch(x0,grd_xarr,grd_nx+1,.false.)
           i = ix
           iy = j
           iz = k
!-- sampling direction helpers
           r1 = rnd_r(rnd_state)
           om0 = pc_pi2*r1
        else
!-- curved surface
           ptcl%x = grd_xarr(grd_nx+1)
           x0 = ptcl%x
           r1 = rnd_r(rnd_state)
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
!-- setting cell index
           ix = i
           iy = binsrch(y0,grd_yarr,grd_ny+1,.false.)
           j = iy
!-- sampling direction helpers
           r1 = rnd_r(rnd_state)
           mu0 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),-mu0)
           if(om0<0d0) om0=om0+pc_pi2
        endif
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
           r1 = rnd_r(rnd_state)
           ptcl%y = r1*grd_yarr(grd_ny+1)+(1d0-r1)*grd_yarr(1)
           y0 = ptcl%y
           r1 = rnd_r(rnd_state)
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           ix = i
           iy = binsrch(y0,grd_yarr,grd_ny+1,.false.)
           j = iy
           iz = binsrch(z0,grd_zarr,grd_nz+1,.false.)
           k = iz
!-- sampling direction helpers
           r1 = rnd_r(rnd_state)
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu1,mu0)
           if(om0<0d0) om0=om0+pc_pi2
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        elseif(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- y surface
           r1 = rnd_r(rnd_state)
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           if(j==1) then
              ptcl%y = grd_yarr(j)
           else
              ptcl%y = grd_yarr(j+1)
              mu0 = -mu0
           endif
           y0 = ptcl%y
           r1 = rnd_r(rnd_state)
           ptcl%z = r1*grd_zarr(grd_nz+1)+(1d0-r1)*grd_zarr(1)
           z0 = ptcl%z
!-- setting cell index
           ix = binsrch(x0,grd_xarr,grd_nx+1,.false.)
           i = ix
           iy = j
           iz = binsrch(z0,grd_zarr,grd_nz+1,.false.)
           k = iz
!-- sampling direction helpers
           r1 = rnd_r(rnd_state)
           mu1 = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om0 = atan2(mu0,mu1)
           if(om0<0d0) om0=om0+pc_pi2
           mu0 = sqrt(1d0-mu0**2)*sin(pc_pi2*r1)

        else
!-- z surface
           r1 = rnd_r(rnd_state)
           ptcl%x = r1*grd_xarr(grd_nx+1)+(1d0-r1)*grd_xarr(1)
           x0 = ptcl%x
           r1 = rnd_r(rnd_state)
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
           ix = binsrch(x0,grd_xarr,grd_nx+1,.false.)
           i = ix
           iy = binsrch(y0,grd_yarr,grd_ny+1,.false.)
           j = iy
           iz = k
!-- sampling azimuthal angle of direction
           r1 = rnd_r(rnd_state)
           om0 = pc_pi2*r1
        endif

        ptcl%mu = mu0
        ptcl%om = om0!}}}

!-- 1D (outer surface)
     case(11)
!-- calculating position!{{{
        ptcl%x = grd_xarr(i+1)
        x0 = ptcl%x
!-- setting cell index
        ix = i
        iy = j
        iz = k
        ptcl%mu = -mu0!}}}
     endselect

!-- particle properties are saved in comoving frame
     ptcl%e = esurfpart
     ptcl%e0 = esurfpart
     ptcl%wl = wl0

!-- DDMC
     if(.not.lhelp) then
        ptcl%e = ptcl%e*P
        ptcl%e0 = ptcl%e0*P
     endif

!-- save particle result
!-----------------------
     prt_particles(ivac) = ptcl


  enddo


end subroutine boundary_source
