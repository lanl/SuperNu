subroutine leakage_opacity2

  use inputparmod
  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC 2D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: iz, ir, ig
  real*8 :: thelp, help
  real*8 :: speclump, specval
  real*8 :: ppl, ppr
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)
!
!-- setting vel-space helper
  if(in_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!
!-- calculating leakage opacities
  gas_opacleak = 0d0
  do iz=1,gas_nz
  do ir=1,gas_nr
!
!-- initializing Planck integral
     speclump = 0d0
     do ig = 1, gas_ng
!-- finding lumpable groups
        if(gas_cap(ig,ir,iz)*min(gas_drarr(ir),gas_dzarr(iz))*thelp >= &
             prt_taulump) then
!-- summing lumpable Planck function integrals
           speclump = speclump + gas_capgrey(ir,iz)*gas_emitprob(ig,ir,iz)/&
                gas_cap(ig,ir,iz)
        endif
     enddo !ig
!-- lumping opacity
     do ig = 1, gas_ng
        if(gas_cap(ig,ir,iz)*min(gas_drarr(ir),gas_dzarr(iz))*thelp >= &
             prt_taulump) then
!
!-- obtaining spectral weight
           specval = gas_capgrey(ir,iz)*gas_emitprob(ig,ir,iz)/&
                gas_cap(ig,ir,iz)

!
!-- calculating inward leakage opacity
           if(ir==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir-1,iz)+ &
                 gas_sig(ir-1,iz))*min(gas_drarr(ir-1),gas_dzarr(iz)) * &
                 thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir,iz)+gas_sig(ir,iz))*gas_drarr(ir)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(1,ir,iz)=gas_opacleak(1,ir,iz)+(specval/speclump)*&
                   ppl*(thelp*gas_rarr(ir))/ &
                   ((gas_rarr(ir+1)**2-gas_rarr(ir)**2)*thelp**2)
           else
!-- DDMC interior
              help = ((gas_sig(ir,iz)+gas_cap(ig,ir,iz))*gas_drarr(ir)+&
                   (gas_sig(ir-1,iz)+gas_cap(ig,ir-1,iz))*gas_drarr(ir-1))*thelp
              gas_opacleak(1,ir,iz)=gas_opacleak(1,ir,iz)+(specval/speclump)*&
                   (4d0/3d0)*(thelp*gas_rarr(ir))/ &
                   (help*(gas_rarr(ir+1)**2-gas_rarr(ir)**2)*thelp**2)
           endif

!
!-- calculating outward leakage opacity
           if(ir==gas_nr) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir+1,iz)+ &
                   gas_sig(ir+1,iz))*min(gas_drarr(ir+1),gas_dzarr(iz)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir,iz)+gas_sig(ir,iz))*gas_drarr(ir)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(2,ir,iz)=gas_opacleak(2,ir,iz)+(specval/speclump)*&
                   ppr*(thelp*gas_rarr(ir+1))/ &
                   ((gas_rarr(ir+1)**2-gas_rarr(ir)**2)*thelp**2)
           else
!-- DDMC interior
              help = ((gas_sig(ir,iz)+gas_cap(ig,ir,iz))*gas_drarr(ir)+&
                   (gas_sig(ir+1,iz)+gas_cap(ig,ir+1,iz))*gas_drarr(ir+1))*thelp
              gas_opacleak(2,ir,iz)=gas_opacleak(2,ir,iz)+(specval/speclump)*&
                   (4d0/3d0)*(thelp*gas_rarr(ir+1))/ &
                   (help*(gas_rarr(ir+1)**2-gas_rarr(ir)**2)*thelp**2)
           endif

!
!-- calculating downward leakage opacity
           if(iz==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir,iz-1)+ &
                   gas_sig(ir,iz-1))*min(gas_drarr(ir),gas_dzarr(iz-1)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir,iz)+gas_sig(ir,iz))*gas_dzarr(iz)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(3,ir,iz)=gas_opacleak(3,ir,iz)+(specval/speclump)*&
                   0.5d0*ppl/(thelp*gas_dzarr(iz))
           else
!-- DDMC interior
              help = ((gas_sig(ir,iz)+gas_cap(ig,ir,iz))*gas_dzarr(iz)+&
                   (gas_sig(ir,iz-1)+gas_cap(ig,ir,iz-1))*gas_dzarr(iz-1))*thelp
              gas_opacleak(3,ir,iz)=gas_opacleak(3,ir,iz)+(specval/speclump)*&
                   (2d0/3d0)/(help*gas_dzarr(iz)*thelp)
           endif

!
!-- calculating upward leakage opacity
           if(iz==gas_nz) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,ir,iz+1)+ &
                   gas_sig(ir,iz+1))*min(gas_drarr(ir),gas_dzarr(iz+1)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,ir,iz)+gas_sig(ir,iz))*gas_dzarr(iz)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(4,ir,iz)=gas_opacleak(4,ir,iz)+(specval/speclump)*&
                   0.5d0*ppr/(thelp*gas_dzarr(iz))
           else
!-- DDMC interior
              help = ((gas_sig(ir,iz)+gas_cap(ig,ir,iz))*gas_dzarr(iz)+&
                   (gas_sig(ir,iz+1)+gas_cap(ig,ir,iz+1))*gas_dzarr(iz+1))*thelp
              gas_opacleak(4,ir,iz)=gas_opacleak(4,ir,iz)+(specval/speclump)*&
                   (2d0/3d0)/(help*gas_dzarr(iz)*thelp)
           endif
        endif
     enddo !ig
  enddo !ir
  enddo !iz


end subroutine leakage_opacity2
