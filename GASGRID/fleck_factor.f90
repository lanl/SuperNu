subroutine fleck_factor(dtempfrac)

  use inputparmod
  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------
  real*8, intent(in) :: dtempfrac
  integer :: ir
  real*8 :: Um, beta, beta2, dlogsig

!-- calculating modified Fleck factor
  do ir = 1, gas_nr
     Um = gas_vals2(ir)%bcoef*gas_temp(ir)
     if(gas_temp(ir)<=0d0.or.gas_vals2(ir)%bcoef==0d0) then
        beta = 0d0
     elseif(gas_siggrey(ir)<=0d0.or.gas_siggreyold(ir)<=0d0) then
        beta = 4.0*gas_vals2(ir)%ur/Um
     else
        dlogsig = log(gas_siggrey(ir)/gas_siggreyold(ir))/(gas_temp(ir)-dtempfrac*gas_temp(ir))
        if(tsp_it==1) then
           beta2 = min((pc_acoef*in_tempradinit**4-gas_vals2(ir)%ur)*dlogsig/gas_vals2(ir)%bcoef,0d0)
        else
           beta2 = min((gas_vals2(ir)%eraddens-gas_vals2(ir)%ur)*dlogsig/gas_vals2(ir)%bcoef,0d0)
        endif

        beta = 4.0*gas_vals2(ir)%ur/Um-beta2

     endif
     gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_siggrey(ir))

  enddo

end subroutine fleck_factor
