subroutine fleck_factor

  use inputparmod
  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------

  integer :: ir
  real*8 :: Um, beta

  if(.not.in_isbdf2.or.tsp_it==1) then
 !-- calculating Fleck factor: 
     do ir = 1, gas_nr
        Um = gas_vals2(ir)%bcoef*gas_temp(ir)
        if(gas_temp(ir)<=0d0.or.gas_vals2(ir)%bcoef==0d0) then
           beta = 0d0
        else
           beta = 4.0*gas_vals2(ir)%ur/Um
        endif
        gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_siggrey(ir))
     enddo
  else
!-- calculating BDF-2 modified Fleck factor
     do ir = 1, gas_nr
        Um = gas_vals2(ir)%bcoef*gas_temp(ir)
        if(gas_temp(ir)<=0d0.or.gas_vals2(ir)%bcoef==0d0) then
           beta = 0d0
        elseif(gas_temp(ir)==gas_tempold(ir)) then
           beta = 4.0*gas_vals2(ir)%ur/Um
        elseif(gas_siggrey(ir)<=0d0.or.gas_siggreyold(ir)<=0d0) then
           beta = 4.0*gas_vals2(ir)%ur/Um
        elseif(tsp_it>1) then
           beta = (pc_acoef/gas_vals2(ir)%bcoef)*&
                (4d0*gas_temp(ir)**3+(gas_temp(ir)**4) &
                *log(gas_cap(ir)/gas_capold(ir))/ &
                (gas_temp(ir)-gas_temphist(ir,tsp_it-1)))
        else
           beta = 4.0*gas_vals2(ir)%ur/Um
        endif
        gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_siggrey(ir))
     enddo
  endif

end subroutine fleck_factor
