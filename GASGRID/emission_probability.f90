subroutine emission_probability

  use gasgridmod
  use physconstmod
  implicit none

!-----------------------
  !multigroup volume emission probabilities
!-----------------------

  integer :: ir, ig
  real*8 :: x1, x2
  real*8 :: specint

  !Calculating grouped volume emission probabilities:
  if(gas_isanalgrp.and.gas_grptype=='pick') then
     do ir = 1, gas_nr
        gas_emitprobg(1,ir) = gas_ppick(1)*gas_sigmapg(1,ir)/gas_siggrey(ir)
        gas_emitprobg(2,ir) = gas_ppick(2)*gas_sigmapg(2,ir)/gas_siggrey(ir)
        do ig = 3, gas_ng
           gas_emitprobg(ig,ir) = 0d0
        enddo
     enddo
  else
     if(gas_ng==1) then
        do ir = 1, gas_nr
           gas_emitprobg(1,ir) = 1d0
        enddo
     else
        do ir = 1, gas_nr
           do ig = 1, gas_ng
              x1 = (pc_h*pc_c/(pc_ev*gas_wl(ig+1)))/(1d3*gas_vals2(ir)%tempkev)
              x2 = (pc_h*pc_c/(pc_ev*gas_wl(ig)))/(1d3*gas_vals2(ir)%tempkev)
              gas_emitprobg(ig,ir) = 15d0*specint(x1,x2,3)*gas_sigmapg(ig,ir)/ &
                   (gas_siggrey(ir)*pc_pi**4)              
           enddo
        enddo
     endif
  endif

end subroutine emission_probability
