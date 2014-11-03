subroutine fleck_factor(tempalt,siggreyalt)

  use inputparmod
  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------
  real*8, intent(in) :: tempalt(gas_nx,gas_ny,gas_nz)
  real*8, intent(in) :: siggreyalt(gas_nx,gas_ny,gas_nz)
  integer :: i,j,k
  real*8 :: Um, beta, beta2, dlogsig

!-- init (necessary for domain decomposition
  gas_fcoef = 0d0

!-- calculating modified Fleck factor
  do k=1,gas_nz
  do j=1,gas_ny
  do i=1,gas_nx
     Um = gas_vals2(i,j,k)%bcoef*gas_temp(i,j,k)
     if(gas_temp(i,j,k)<=0d0.or.gas_vals2(i,j,k)%bcoef==0d0) then
        beta = 0d0
     elseif(gas_siggrey(i,j,k)<=0d0.or.siggreyalt(i,j,k)<=0d0) then
        beta = 4.0*gas_vals2(i,j,k)%ur/Um
     else
        if(.not.in_ismodimc) then
           beta2 = 0d0
        else
!-- convert from per gram
           dlogsig = log(gas_siggrey(i,j,k)/(siggreyalt(i,j,k)*gas_vals2(i,j,k)%rho))/&
              (gas_temp(i,j,k)-tempalt(i,j,k))
           beta2 = min(0d0,(gas_vals2(i,j,k)%eraddens-gas_vals2(i,j,k)%ur)*&
              dlogsig/gas_vals2(i,j,k)%bcoef)
        endif
        beta = 4.0*gas_vals2(i,j,k)%ur/Um - beta2

     endif
     gas_fcoef(i,j,k) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_siggrey(i,j,k))
  enddo !i
  enddo !j
  enddo !k

end subroutine fleck_factor
