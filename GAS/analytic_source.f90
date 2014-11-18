subroutine analytic_source

  use gridmod
  use mpimod
  use gasmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  use totalsmod
  implicit none

  integer :: i,j,k,l,ll,nhelp
  integer :: l1,l2
  real*8 :: srcren
  real*8 :: thelp, help, xcent, ycent, zcent

  tot_esurf = 0d0
  gas_emitex = 0d0

!-- setting source helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

  l1 = impi*gas_ncell + 1
  l2 = (impi+1)*gas_ncell

  if(in_srctype=='none') then
    return
  elseif(in_srctype=='surf') then
!-- time-integrated surface flux [erg/cm^2]
     help = 0.25d0*pc_acoef*pc_c*tsp_dt*in_srcmax**4 * &
          thelp**2
     select case(in_igeom)
!-- 1D
     case(1)
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
!-- 1D
        case(1)
           ll = 0
           do i = 1, min(in_nheav,grd_nx)
              l = i
              if(l<l1 .or. l>l2) cycle
              ll = ll + 1
              gas_emitex(ll) = in_srcmax * &
                   gas_vol(ll)*tsp_dt/thelp**3
              !write(0,*) impi,ll,gas_emitex(ll),gas_vol(ll)
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
           l = 0
           ll = 0
           do j = 1,grd_ny
           do i = 1,grd_nx
              l = l + 1
              if(l<l1 .or. l>l2) cycle
              ll = ll + 1
              xcent = 0.5d0*(grd_xarr(i+1)+grd_xarr(i))
              ycent = 0.5d0*(grd_yarr(j+1)+grd_yarr(j))
              if(xcent**2+ycent**2<help**2) then
                 gas_emitex(ll) = in_srcmax * &
                      gas_vol(ll)*tsp_dt/thelp**3
              endif
           enddo
           enddo

!-- 3D
        case(3)
!
!-- using min distance to bound
           help = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1) , &
                grd_zarr(grd_nz+1))
           if(help == grd_xarr(grd_nx+1)) then
              nhelp = grd_nx
           elseif(help == grd_yarr(grd_ny+1)) then
              nhelp = grd_ny
           else
              nhelp = grd_nz
           endif
!-- Heaviside radius <= distance to bound
           help = dble(min(in_nheav,nhelp))*help / &
                dble(nhelp)
!-- non-zero source within Heaviside sphere
           l = 0
           ll = 0
           do k = 1,grd_nz
           do j = 1,grd_ny
           do i = 1,grd_nx
              l = l + 1
              if(l<l1 .or. l>l2) cycle
              ll = ll + 1
              xcent = 0.5d0*(grd_xarr(i+1)+grd_xarr(i))
              ycent = 0.5d0*(grd_yarr(j+1)+grd_yarr(j))
              zcent = 0.5d0*(grd_zarr(k+1)+grd_zarr(k))
              if(xcent**2+ycent**2+zcent**2<help**2) then
                 gas_emitex(ll) = in_srcmax * &
                      gas_vol(ll)*tsp_dt/thelp**3
              endif
           enddo
           enddo
           enddo
        endselect
     endif
!-- no temp source for heav (matsrc=0.0)
!--
     !!}}}
  elseif(in_srctype=='strt') then
     !Linear source profile!{{{
     if(grd_ny>1) stop 'analytic_source: strt: no 2D'
     ll = 0
     do i=1,grd_nx
        l = i
        if(l<l1 .or. l>l2) cycle
        ll = ll + 1
        srcren = in_srcmax*(grd_xarr(grd_nx+1)- &
             0.5d0*(grd_xarr(i)+grd_xarr(i+1)))/ & 
             (grd_xarr(grd_nx+1)-grd_xarr(1))
        gas_emitex(ll) = srcren * gas_vol(ll)*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(in_srctype=='manu') then
     !!{{{
     if(grd_ny>1) stop 'analytic_source: manu: no 2D'
!
!-- radiation source
     call generate_manuradsrc(in_totmass,in_sigcoef,tsp_t,tsp_dt)
!
!-- temperature source
     call generate_manutempsrc(in_totmass,in_sigcoef,tsp_t,tsp_dt)
!     
     !}}}
  else
     stop 'analytic_source: in_srctype invalid'
  endif

  !write(*,*) grd_siggrey(grd_nx), grd_cap(1,grd_nx)

end subroutine analytic_source
