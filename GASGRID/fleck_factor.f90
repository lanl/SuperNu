subroutine fleck_factor

  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------

  integer :: ir
  real*8 :: Um, beta

  !Calculating Fleck factor: 
  do ir = 1, gas_nr
     Um = gas_vals2(ir)%bcoef*gas_temp(ir)
     if(gas_temp(ir)<=0d0.or.gas_vals2(ir)%bcoef==0d0) then
        beta = 0d0
     else
        beta = 4.0*gas_vals2(ir)%ur/Um
     endif
     gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_siggrey(ir))
  enddo

end subroutine fleck_factor
