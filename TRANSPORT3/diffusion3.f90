subroutine diffusion3(ptcl,isvacant)

  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  use totalsmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch
!
!-- lumped quantities
  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleakllump, opacleakrlump !x-leakage
  real*8 :: opacleakdlump, opacleakulump !y-leakage
  real*8 :: opacleakblump, opacleaktlump !z-leakage
  real*8 :: mfphelp, pp
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer :: glumps(grd_ng)
  real*8 :: glumpinv,dtinv,capinv(grd_ng)
  real*8 :: help
!
  integer,pointer :: ix,iy,iz
  real*8,pointer :: x,y,z,xi,om,ep,ep0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dx(l) = grd_zarr(l+1) - grd_zarr(l)

  ix => ptcl%zsrc
  iy => ptcl%iy
  iz => ptcl%iz
  x => ptcl%rsrc
  y => ptcl%y
  z => ptcl%z
  xi => ptcl%musrc
  om => ptcl%om
  ep => ptcl%esrc
  ep0 => ptcl%ebirth
  wl => ptcl%wlsrc
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  capinv = 1d0/grd_cap(:,ix,iy,iz)

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- looking up initial group
  g = binsrch(wl,grd_wl,grd_ng+1,in_ng)
!-- checking group bounds
  if(g>grd_ng.or.g<1) then
     if(g==grd_ng+1) then
        g = grd_ng
     elseif(g==0) then
        g = 1
     else
        stop 'diffusion3: particle group invalid'
     endif
  endif

end subroutine diffusion3
