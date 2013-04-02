subroutine analytic_source

  use gasgridmod
  use physconstmod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2
  real*8 :: specint

  if(gas_srctype=='heav') then
     !Grey Heaviside source (uniform source sphere)

  elseif(gas_srctype=='strt') then
     !Grey linear source profile

  elseif(gas_srctype=='manu') then
     !Manufactured Source (for gas_grptype='line')
     stop 'analytic_source: manu not yet available'
  else
     stop 'analytic_source: gas_srctype invalid'
  endif

end subroutine analytic_source
