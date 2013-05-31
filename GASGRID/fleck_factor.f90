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
     Um = gas_vals2(ir)%bcoef*gas_vals2(ir)%tempkev
     beta = 4.0*gas_vals2(ir)%ur/Um
     gas_fcoef(ir) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_sigmap(ir))
  enddo

end subroutine fleck_factor
