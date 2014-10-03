subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: ir, ig
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
        do ir = 1, min(gas_nheav,gas_nr)
           gas_emitex(ir) = gas_srcmax * &
                gas_vals2(ir)%vol*tsp_dt/thelp**3
        enddo
     endif
!-- no temp source for heav (matsrc=0.0)
!--
     !!}}}
  elseif(gas_srctype=='strt') then
     !Linear source profile!{{{
     do ir = 1, gas_nr
        srcren = gas_srcmax*(gas_rarr(gas_nr+1)- &
             0.5d0*(gas_rarr(ir)+gas_rarr(ir+1)))/ & 
             (gas_rarr(gas_nr+1)-gas_rarr(1))
        gas_emitex(ir) = srcren * gas_vals2(ir)%vol*tsp_dt
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

  !write(*,*) gas_siggrey(gas_nr), gas_cap(1,gas_nr)

end subroutine analytic_source
