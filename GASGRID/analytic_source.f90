subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: i, ig
  real*8 :: x1, x2, x3, x4, srcren
  real*8 :: thelp

  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)

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
        do i = 1, min(gas_nheav,gas_nx)
           gas_emitex(i,1,1) = gas_srcmax * &
                gas_vals2(i,1,1)%vol*tsp_dt/thelp**3
        enddo
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
        gas_emitex(i,1,1) = srcren * gas_vals2(i,1,1)%vol*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(gas_srctype=='manu') then
     !!{{{
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
