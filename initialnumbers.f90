subroutine initialnumbers

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine computes the distribution of initial particle energy
  !before the first time step.  A fraction of the total initial particle
  !number is given to each cell based on the amount of inital radiative
  !energy profile.
!##################################################
  
  integer :: ir, ig
  integer :: nvolinittot, nvolinitapp
  real*8 :: wl0, mu0, Ep0, r0
  real*8 :: help
  real*8 :: exsumg, rrcenter, etotinit,denom2
  real*8 :: x1,x2,x3,x4
  !
  nvolinitapp = 0
  etotinit = 0d0
  !
  gas_vals2(:)%eraddens=0d0
  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)

  if(gas_isvelocity) then
     help = gas_velout*tsp_texp
  else
     help = gas_l0+gas_lr
  endif
                
     
     do ir = 1, gas_nr
        nvolinit(ir)=nint(gas_vals2(ir)%eraddens*gas_vals2(ir)%volr*help**3 &
             *nvolinittot/etotinit)
        nvolinitapp = nvolinitapp+nvolinit(ir)
     enddo


end subroutine initialnumbers
