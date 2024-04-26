subroutine rprocess_fitfrm_source(it,nt,tcenter)
  use gasmod
  use miscmod, only:binsrch
  use rprocmod, only:v_grid
  use gridmod, only:grd_xarr, grd_nx
  implicit none
  integer,intent(in) :: it, nt
  real*8,intent(in) :: tcenter
!--------------------------------------------------
!- Calculate heating rate from analytic fit formula (interpolated).
!- See Rosswog & Korobkin (2023) for details
!--------------------------------------------------
  integer,parameter :: num_rad_types = 4
  integer :: icol, i
  real*8 :: helparr(gas_ncell)
  real*8 :: therm_fracs(gas_ncell,num_rad_types)

!-- reset material source
  gas_matsrc = 1d10
  gas_matsrc = gas_matsrc * gas_rho * gas_ye * (1d0 - grd_xarr / grd_xarr(grd_nx))

!-- initialize helper and thermalization fraction array
  helparr = 0d0
  therm_fracs = 0d0

!-- set the helper array
  where(gas_rho>0d0) helparr = 2d0/(tcenter*gas_rho)

end subroutine rprocess_fitfrm_source
