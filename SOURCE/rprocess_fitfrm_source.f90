subroutine rprocess_fitfrm_source(it,nt,tcenter)
  use physconstmod, only:pc_msun,pc_pi,pc_c
  use gasmod, only:gas_ye0,gas_mass,gas_rho,gas_matsrc,gas_ncell
  use rprocmod, only:heating_rate,v_grid
  use gridmod, only:grd_vol
  implicit none
  integer,intent(in) :: it, nt
  real*8,intent(in) :: tcenter
!--------------------------------------------------
!- Calculate heating rate from analytic fit formula (interpolated).
!- See Rosswog & Korobkin (2023) for details
!--------------------------------------------------
  integer,parameter :: num_rad_types = 4
  integer :: icol, i
  real*8 :: helparr(gas_ncell),vexp(gas_ncell)
  real*8 :: therm_fracs(gas_ncell,num_rad_types)
  real*8, parameter :: M_ref = 0.05d0*pc_msun  !< reference mass when computing vexp

!-- reset material source
  vexp = 0.1d0 ! default value of the expansion velocity
  if (tcenter > 0) &
     where(gas_mass > 0) &
        vexp = (3*M_ref/(4*pc_pi) * grd_vol/gas_mass)**(1d0/3d0)/tcenter
  vexp = vexp/pc_c
  vexp = max(v_grid(1), vexp)
  vexp = min(v_grid(size(v_grid)), vexp)

  where(gas_rho > 0d0) &
     gas_matsrc = gas_rho*heating_rate(vexp, gas_ye0, tcenter)

end subroutine rprocess_fitfrm_source
