subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use inputparmod
  use manufacmod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2, x3, x4, srcren
  real*8 :: specint
  !
  real*8 :: xx3, xx4
  real*8 :: aleff1 = 0.5d0

  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)

  gas_emitex = 0d0

  if(gas_srctype=='none') then
    return
  elseif(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)!{{{
     if (tsp_texp<=(in_tfirst+gas_theav)*pc_day) then
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1) 
              x4 = 1d0/gas_wl(ig)
              gas_emitex(ig,ir)=gas_srcmax*(x4-x3)/(x2-x1)
              if(gas_isvelocity) then
                 gas_emitex(ig,ir)=gas_emitex(ig,ir)/tsp_texp**3
              endif
           enddo
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!
        enddo
        do ir = gas_nheav+1, gas_nr
           do ig = 1, gas_ng
              gas_emitex(ig,ir)=0d0
           enddo
        enddo
     else
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              gas_emitex(ig,ir)=0d0
           enddo
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
        do ig = 1, gas_ng
           x3 = 1d0/gas_wl(ig+1) 
           x4 = 1d0/gas_wl(ig)
           gas_emitex(ig,ir)=srcren*(x4-x3)/(x2-x1)
        enddo
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!
!-- no temp source for strt (matsrc=0.0)
!--
     enddo!}}}
  elseif(gas_srctype=='manu') then
     !!{{{
!
!-- radiation source
     call generate_manuradsrc(in_totmass,in_sigcoef,tsp_texp,tsp_dt)
!
!-- temperature source
     call generate_manutempsrc(in_totmass,in_sigcoef,tsp_texp,tsp_dt)
!     
     !}}}
  else
     stop 'analytic_source: gas_srctype invalid'
  endif

  !write(*,*) gas_emitex(1,1), gas_emitex(2,1), gas_emitex(3,1)  
  !write(*,*) gas_emitex(1,:)
  !write(*,*)
  !write(*,*) gas_emitex(2,:)
  !write(*,*) gas_siggrey(gas_nr), gas_cap(1,gas_nr)

end subroutine analytic_source
