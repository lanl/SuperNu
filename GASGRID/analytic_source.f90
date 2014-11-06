subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: i,j, nhelp
  real*8 :: srcren
  real*8 :: thelp, help, xcent, ycent

  gas_emitex = 0d0

!-- setting source helper
  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

  if(gas_srctype=='none') then
    return
  elseif(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)!{{{
     if (tsp_t<=(in_tfirst+gas_theav)*pc_day) then
        select case(in_igeom)
!-- 1D
        case(1)
           do i = 1, min(gas_nheav,gas_nx)
              gas_emitex(i,1,1) = gas_srcmax * &
                   gas_vol(i,1,1)*tsp_dt/thelp**3
           enddo

!-- 2D
        case(2)
!
!-- using min distance to cylinder bound
           help = min(gas_xarr(gas_nx+1),gas_yarr(gas_ny+1))
           if(help == gas_xarr(gas_nx+1)) then
              nhelp = gas_nx
           else
              nhelp = gas_ny
           endif
!-- Heaviside radius <= distance to cylinder bound
           help = dble(min(gas_nheav,nhelp))*help / &
                dble(nhelp)
!-- non-zero source within Heaviside sphere
           do j = 1,gas_ny
              do i = 1,gas_nx
                 xcent = 0.5d0*(gas_xarr(i+1)+gas_xarr(i))
                 ycent = 0.5d0*(gas_yarr(j+1)+gas_yarr(j))
                 if(xcent**2+ycent**2<help**2) then
                    gas_emitex(i,j,1) = gas_srcmax * &
                         gas_vol(i,j,1)*tsp_dt/thelp**3
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
     if(gas_ny>1) stop 'analytic_source: strt: no 2D'
     do i=1,gas_nx
        srcren = gas_srcmax*(gas_xarr(gas_nx+1)- &
             0.5d0*(gas_xarr(i)+gas_xarr(i+1)))/ & 
             (gas_xarr(gas_nx+1)-gas_xarr(1))
        gas_emitex(i,1,1) = srcren * gas_vol(i,1,1)*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(gas_srctype=='manu') then
     !!{{{
     if(gas_ny>1) stop 'analytic_source: manu: no 2D'
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

  !write(*,*) gas_siggrey(gas_nx), gas_cap(1,gas_nx)

end subroutine analytic_source
