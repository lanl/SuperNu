subroutine analytic_initial

  use gridmod
  use gasgridmod
  use inputparmod
  use physconstmod
  use timestepmod
  use manufacmod
  use profiledatamod
  implicit none

  real*8 :: trad(grd_nx,grd_ny,grd_nz)

!###############################################
! This subroutines attributes radiation energy to
! each cell and group depeding on user specification
! of gas_srctype
!###############################################
!
  grd_evolinit = 0d0
!
!-- initial radiation energy
  trad = in_tempradinit
!
!-- map radiation temperature to grd_evolinit
  grd_evolinit = pc_acoef*trad**4 * grd_vol
!--
!
!-- source specific initial conditions (overrides gas_inittyp)
!-- currently only supplying nonzero for gas_srctype=manu
  if(gas_srctype=='none') then
     if(gas_opacanaltype=='pick') then
!-- tstd initial energy profile currently approximation
        stop 'analytic_initial: gas_opacanaltype==pick not implemented'
     else
        return
     endif
  elseif(gas_srctype=='heav') then
     return
  elseif(gas_srctype=='strt') then
     return
  elseif(gas_srctype=='manu') then
     call init_manuprofile
  else
     stop 'analytic_initial: invalid gas_srctype'
  endif

end subroutine analytic_initial
