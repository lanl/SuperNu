subroutine leakage_opacity

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC lumped leakage opacities.
!##################################################
  logical :: missive = .false.
  logical :: lhelp
  integer :: ir, ig
  real*8 :: help
  real*8 :: thelp, ppl, ppr, specval, speclump
!
!-- setting vel-space helper
  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- calculating leakage opacities
  gas_opacleak = 0d0
  do ir = 1, gas_nr
!
!-- initializing Planck integral
     speclump = 0d0
     do ig = 1, gas_ng
!-- finding lumpable groups
        if(gas_cap(ig,ir)*gas_drarr(ir)*thelp>=prt_taulump) then
!-- summing lumpable Planck function integrals
           speclump = speclump + gas_siggrey(ir)*gas_emitprob(ig,ir)/&
                gas_cap(ig,ir)
        endif
     enddo
!-- lumping opacity
     do ig = 1, gas_ng
        if(gas_cap(ig,ir)*gas_drarr(ir)*thelp>=prt_taulump) then
!
!-- obtaining spectral weight
           specval = gas_siggrey(ir)*gas_emitprob(ig,ir)/&
                gas_cap(ig,ir)
!
!-- case differentiation
           if(ir==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir-1)+ &
                 gas_sig(ir-1))*gas_drarr(ir-1)*thelp<prt_tauddmc
           endif
!
!-- calculating inward leakage opacity
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir)+gas_sig(ir))*gas_drarr(ir)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(1,ir)=gas_opacleak(1,ir)+(specval/speclump)*&
                   1.5d0*ppl*(thelp*gas_rarr(ir))**2/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
           else
!-- DDMC interior
              help = ((gas_sig(ir)+gas_cap(ig,ir))*gas_drarr(ir)+&
                   (gas_sig(ir-1)+gas_cap(ig,ir-1))*gas_drarr(ir-1))*thelp
              gas_opacleak(1,ir)=gas_opacleak(1,ir)+(specval/speclump)*&
                   2.0d0*(thelp*gas_rarr(ir))**2/ &
                   (help*3d0*gas_vals2(ir)%vol/pc_pi4)
           endif
!
!-- case differentation
           if(ir==gas_nr) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir+1)+ &
                 gas_sig(ir+1))*gas_drarr(ir+1)*thelp<prt_tauddmc
           endif
!
!-- calculating outward leakage opacity
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir)+gas_sig(ir))*gas_drarr(ir)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(2,ir)=gas_opacleak(2,ir)+(specval/speclump)*&
                   1.5d0*ppr*(thelp*gas_rarr(ir+1))**2/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
           else
!-- DDMC interior
              help = ((gas_sig(ir)+gas_cap(ig,ir))*gas_drarr(ir)+&
                   (gas_sig(ir+1)+gas_cap(ig,ir+1))*gas_drarr(ir+1))*thelp
              gas_opacleak(2,ir)=gas_opacleak(2,ir)+(specval/speclump)*&
                   2.0d0*(thelp*gas_rarr(ir+1))**2/ &
                   (help*3d0*gas_vals2(ir)%vol/pc_pi4)
           endif
        endif
     enddo
  enddo
  

end subroutine leakage_opacity
