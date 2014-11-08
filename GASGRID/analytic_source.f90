subroutine analytic_source

  use gridmod
  use mpimod
  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: i,j,l,ll,nhelp
  integer :: l1,l2
  real*8 :: srcren
  real*8 :: thelp, help, xcent, ycent

  dd_emitex = 0d0

!-- setting source helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

  l1 = impi*dd_ncell + 1
  l2 = (impi+1)*dd_ncell

  if(gas_srctype=='none') then
    return
  elseif(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)!{{{
     if (tsp_t<=(in_tfirst+gas_theav)*pc_day) then
        select case(in_igeom)
!-- 1D
        case(1)
           ll = 0
           do i = 1, min(gas_nheav,grd_nx)
              l = i
              if(l<l1 .or. l>l2) cycle
              ll = ll + 1
              dd_emitex(ll) = gas_srcmax * &
                   dd_vol(ll)*tsp_dt/thelp**3
              !write(0,*) impi,ll,dd_emitex(ll),dd_vol(ll)
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
           help = dble(min(gas_nheav,nhelp))*help / &
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
                 dd_emitex(ll) = gas_srcmax * &
                      dd_vol(ll)*tsp_dt/thelp**3
              endif
           enddo
           enddo
        case(3)
           stop 'analytic_source: no 3D transport'
        endselect
     endif
!-- no temp source for heav (matsrc=0.0)
!--
     !!}}}
  elseif(gas_srctype=='strt') then
     !Linear source profile!{{{
     if(grd_ny>1) stop 'analytic_source: strt: no 2D'
     ll = 0
     do i=1,grd_nx
        l = i
        if(l<l1 .or. l>l2) cycle
        ll = ll + 1
        srcren = gas_srcmax*(grd_xarr(grd_nx+1)- &
             0.5d0*(grd_xarr(i)+grd_xarr(i+1)))/ & 
             (grd_xarr(grd_nx+1)-grd_xarr(1))
        dd_emitex(ll) = srcren * dd_vol(ll)*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(gas_srctype=='manu') then
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
     stop 'analytic_source: gas_srctype invalid'
  endif

  !write(*,*) grd_siggrey(grd_nx), grd_cap(1,grd_nx)

end subroutine analytic_source
