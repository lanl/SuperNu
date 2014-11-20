subroutine fleck_factor(tempalt,siggreyalt)

  use inputparmod
  use gasmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------
  real*8, intent(in) :: tempalt(gas_ncell)
  real*8, intent(in) :: siggreyalt(gas_ncell)
  integer :: i
  real*8 :: Um, beta, beta2, dlogsig

!-- init (necessary for domain decomposition
  gas_fcoef = 0d0

!-- calculating modified Fleck factor
  do i=1,gas_ncell
     Um = gas_bcoef(i)*gas_temp(i)
     if(gas_temp(i)<=0d0.or.gas_bcoef(i)==0d0) then
        beta = 0d0
     elseif(gas_capgrey(i)<=0d0.or.siggreyalt(i)<=0d0) then
        beta = 4.0*gas_ur(i)/Um
     else
        if(.not.in_ismodimc) then
           beta2 = 0d0
        else
!-- convert from per gram
           dlogsig = log(gas_capgrey(i)/(siggreyalt(i)*gas_rho(i)))/&
              (gas_temp(i)-tempalt(i))
           beta2 = min(0d0,(gas_eraddens(i)-gas_ur(i))*&
              dlogsig/gas_bcoef(i))
        endif
        beta = 4.0*gas_ur(i)/Um - beta2

     endif
     gas_fcoef(i) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*gas_capgrey(i))
  enddo !i

end subroutine fleck_factor
