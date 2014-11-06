subroutine emission_probability

  use gasgridmod
  use physconstmod
  use miscmod, only:specint
  implicit none

!-----------------------
  !multigroup volume emission probabilities
!-----------------------

  integer :: i,ig
  real*8 :: x1, x2

!-- init
  dd_emitprob = 0d0

!-- Calculating grouped volume emission probabilities:
  if(gas_opacanaltype=='pick') then
     do i=1,dd_ncell
        dd_emitprob(1,i) = gas_ppick(1)*dd_cap(1,i)/dd_siggrey(i)
        dd_emitprob(2,i) = gas_ppick(2)*dd_cap(2,i)/dd_siggrey(i)
!       dd_emitprob(3:gas_ng,i) = 0d0  !-- not necessary
     enddo !i
  else
     if(gas_ng==1) then
        dd_emitprob = 1d0
     else
        do i=1,dd_ncell
           do ig=1,gas_ng
              x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*dd_temp(i))
              x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*dd_temp(i))
              if(dd_siggrey(i)<=0d0) then
!                dd_emitprob(ig,i) = 0d0  !-- not necessary
              else
                 dd_emitprob(ig,i) = 15d0*specint(x1,x2,3)*dd_cap(ig,i)/ &
                      (dd_siggrey(i)*pc_pi**4)
              endif
           enddo !ig
        enddo !i
     endif
  endif

end subroutine emission_probability
