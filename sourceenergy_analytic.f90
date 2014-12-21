subroutine sourceenergy_analytic

  use gridmod
  use mpimod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  use totalsmod
  implicit none

  integer :: i,j,k,nhelp
  real*8 :: srcren
  real*8 :: thelp, help, xcent, ycent, zcent

  tot_esurf = 0d0
  grd_emitex = 0d0

!-- setting source helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

  if(in_srctype=='none') then
    return
  elseif(in_srctype=='surf') then
!-- time-integrated surface flux [erg/cm^2]
     help = 0.25d0*pc_acoef*pc_c*tsp_dt*in_srcmax**4 * &
          thelp**2
     select case(in_igeom)
!-- [123]D spherical
     case(1,11)
        tot_esurf = help*pc_pi4*grd_xarr(grd_nx+1)**2
!-- 2D
     case(2)
        if(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- flat surface
           tot_esurf = help*pc_pi*grd_xarr(grd_nx+1)**2
        else
!-- curved surface
           tot_esurf = help*2d0*pc_pi*grd_xarr(grd_nx+1) * &
                (grd_yarr(grd_ny+1)-grd_yarr(1))
        endif
!-- 3D
     case(3)
        if(in_surfsrcloc=='in'.or.in_surfsrcloc=='out') then
!-- x surface
           tot_esurf = help*(grd_yarr(grd_ny+1)-grd_yarr(1)) * &
                (grd_zarr(grd_nz+1)-grd_zarr(1))
        elseif(in_surfsrcloc=='down'.or.in_surfsrcloc=='up') then
!-- y surface
           tot_esurf = help*(grd_xarr(grd_nx+1)-grd_xarr(1)) * &
                (grd_zarr(grd_nz+1)-grd_zarr(1))
        elseif(in_surfsrcloc=='botm'.or.in_surfsrcloc=='top') then
!-- z surface
           tot_esurf = help*(grd_xarr(grd_nx+1)-grd_xarr(1)) * &
                (grd_yarr(grd_ny+1)-grd_yarr(1))
        endif
     endselect
  elseif(in_srctype=='heav') then
     !Heaviside source (uniform source sphere)!{{{
     if (tsp_t<=(in_tfirst+in_theav)*pc_day) then
        select case(in_igeom)
!-- [123]D spherical
        case(1,11)
           do i=1,min(in_nheav,grd_nx)
              grd_emitex(i,:,:) = in_srcmax * &
                   grd_vol(i,:,:)*tsp_dt/thelp**3
              !write(0,*) impi,grd_emitex(i,1,1),grd_vol(i,1,1)
           enddo

!-- 2D
        case(2)
!
!-- using min distance to cylinder bound
           help = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1))
           if(help == grd_xarr(grd_nx+1)) then
              nhelp = grd_nx
           else
              nhelp = grd_ny
           endif
!-- Heaviside radius <= distance to cylinder bound
           help = dble(min(in_nheav,nhelp))*help / &
                dble(nhelp)
!-- non-zero source within Heaviside sphere
           do j=1,grd_ny
           ycent = 0.5d0*(grd_yarr(j+1)+grd_yarr(j))
           do i=1,grd_nx
              xcent = 0.5d0*(grd_xarr(i+1)+grd_xarr(i))
              if(xcent**2+ycent**2<help**2) then
                 grd_emitex(i,j,1) = in_srcmax * &
                      grd_vol(i,j,1)*tsp_dt/thelp**3
              endif
           enddo
           enddo
!-- adjusting bulk source energy
!           grd_emitex=grd_emitex*in_srcmax*tsp_dt*pc_pi43*help**3 / &
!                sum(grd_emitex)

!-- 3D
        case(3)
!
!-- using min distance to bound
           help = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1) , &
                grd_zarr(grd_nz+1))
           if(help == grd_xarr(grd_nx+1)) then
              nhelp = grd_nx/2
           elseif(help == grd_yarr(grd_ny+1)) then
              nhelp = grd_ny/2
           else
              nhelp = grd_nz/2
           endif
!-- Heaviside radius <= distance to bound
           help = dble(min(in_nheav,nhelp))*help / &
                dble(nhelp)
!-- non-zero source within Heaviside sphere
           do k=1,grd_nz
           zcent = 0.5d0*(grd_zarr(k+1)+grd_zarr(k))
           do j=1,grd_ny
           ycent = 0.5d0*(grd_yarr(j+1)+grd_yarr(j))
           do i=1,grd_nx
              xcent = 0.5d0*(grd_xarr(i+1)+grd_xarr(i))
              if(xcent**2+ycent**2+zcent**2<help**2) then
                 grd_emitex(i,j,k) = in_srcmax * &
                      grd_vol(i,j,k)*tsp_dt/thelp**3
              endif
           enddo
           enddo
           enddo
!-- adjusting bulk source energy
!           grd_emitex=grd_emitex*in_srcmax*tsp_dt*pc_pi43*help**3 / &
!                sum(grd_emitex)
        endselect
     endif
!-- no temp source for heav (matsrc=0.0)
!--
     !!}}}
  elseif(in_srctype=='strt') then
     !Linear source profile!{{{
     if(grd_ny>1) stop 'sourceenergy_analytic: strt: no 2D'
     do i=1,grd_nx
        srcren = in_srcmax*(grd_xarr(grd_nx+1)- &
             0.5d0*(grd_xarr(i)+grd_xarr(i+1)))/ & 
             (grd_xarr(grd_nx+1)-grd_xarr(1))
        grd_emitex(i,1,1) = srcren * grd_vol(i,1,1)*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(in_srctype=='manu') then
     !!{{{
     if(grd_ny>1) stop 'sourceenergy_analytic: manu: no 2D'
!
!-- radiation source
     call generate_manuradsrc(in_totmass,in_gas_capcoef,tsp_t,tsp_dt)
!
!--  NOW DONE IN SOURCEENERGY
!!-- temperature source
!     call generate_manutempsrc(in_totmass,in_gas_capcoef,tsp_t,tsp_dt)
!     
     !}}}
  else
     stop 'sourceenergy_analytic: in_srctype invalid'
  endif

!
!-- add analytic radiation source tot total
  if(impi==impi0) tot_eext = tot_eext + sum(grd_emitex)

end subroutine sourceenergy_analytic
