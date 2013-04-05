subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use particlemod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2, x3, x4, srcren
  real*8 :: specint

  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)
  if(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)
     if (tsp_time<=gas_theav*pc_day) then
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1) 
              x4 = 1d0/gas_wl(ig)
              gas_exsource(ig,ir)=gas_srcmax*(x4-x3)/(x2-x1)
           enddo
        enddo
        do ir = gas_nheav+1, gas_nr
           do ig = 1, gas_ng
              gas_exsource(ig,ir)=0d0
           enddo
        enddo
     else
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              gas_exsource(ig,ir)=0d0
           enddo
        enddo
        !Ryan W.: Temporary fix (?) for over generation,
        !resetting prt_ns
        if(gas_sigcoef==0d0) then
           prt_ns = 1
        endif
     endif
     !
  elseif(gas_srctype=='strt') then
     !Linear source profile
     do ir = 1, gas_nr
        srcren = gas_srcmax*(gas_rarr(gas_nr+1)- &
             0.5d0*(gas_rarr(ir)+gas_rarr(ir+1)))/ & 
             (gas_rarr(gas_nr+1)-gas_rarr(1))
        do ig = 1, gas_ng
           x3 = 1d0/gas_wl(ig+1) 
           x4 = 1d0/gas_wl(ig)
           gas_exsource(ig,ir)=srcren*(x4-x3)/(x2-x1)
        enddo
     enddo
  elseif(gas_srctype=='manu') then
     !Manufactured Source (for gas_grptype='line')
     stop 'analytic_source: manu not yet available'
  else
     stop 'analytic_source: gas_srctype invalid'
  endif

end subroutine analytic_source
