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
  if(gas_opacanaltype=='pick') then
     do ir = 1, gas_nr
        gas_emitprob(1,ir) = gas_ppick(1)*gas_cap(1,ir)/gas_siggrey(ir)
        gas_emitprob(2,ir) = gas_ppick(2)*gas_cap(2,ir)/gas_siggrey(ir)
        do ig = 3, gas_ng
           gas_emitprob(ig,ir) = 0d0
        enddo
     enddo
  else
     if(gas_ng==1) then
        do ir = 1, gas_nr
           gas_emitprob(1,ir) = 1d0
        enddo
     else
        do ir = 1, gas_nr
           do ig = 1, gas_ng
              x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_vals2(ir)%temp)
              x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_vals2(ir)%temp)
              gas_emitprob(ig,ir) = 15d0*specint(x1,x2,3)*gas_cap(ig,ir)/ &
                   (gas_siggrey(ir)*pc_pi**4)              
           enddo
        enddo
     endif
  endif

end subroutine emission_probability
