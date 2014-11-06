subroutine fleck_factor(tempalt,siggreyalt)

  use inputparmod
  use gasgridmod
  use timestepmod
  use physconstmod
  implicit none

  !-----------------------
  !This subroutine computes the Fleck factor
  !-----------------------
  real*8, intent(in) :: tempalt(dd_ncell)
  real*8, intent(in) :: siggreyalt(dd_ncell)
  integer :: i
  real*8 :: Um, beta, beta2, dlogsig

!-- init (necessary for domain decomposition
  dd_fcoef = 0d0

!-- calculating modified Fleck factor
  do i=1,dd_ncell
     Um = dd_bcoef(i)*dd_temp(i)
     if(dd_temp(i)<=0d0.or.dd_bcoef(i)==0d0) then
        beta = 0d0
     elseif(dd_siggrey(i)<=0d0.or.siggreyalt(i)<=0d0) then
        beta = 4.0*dd_ur(i)/Um
     else
        if(.not.in_ismodimc) then
           beta2 = 0d0
        else
!-- convert from per gram
           dlogsig = log(dd_siggrey(i)/(siggreyalt(i)*dd_rho(i)))/&
              (dd_temp(i)-tempalt(i))
           beta2 = min(0d0,(dd_eraddens(i)-dd_ur(i))*&
              dlogsig/dd_bcoef(i))
        endif
        beta = 4.0*dd_ur(i)/Um - beta2

     endif
     dd_fcoef(i) = 1.0/(1.0+tsp_alpha*beta*pc_c*tsp_dt*dd_siggrey(i))
  enddo !i

end subroutine fleck_factor
