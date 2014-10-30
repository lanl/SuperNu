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
  logical :: lhelp
  integer :: i,j,k, ig
  real*8 :: help
  real*8 :: thelp, ppl, ppr, specval, speclump
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
!
!-- setting vel-space helper
  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- init (necessary for domain decomposition)
  gas_opacleak = 0d0

!
!-- calculating leakage opacities
  do k = 1, gas_nz
  do j = 1, gas_ny
  do i = 1, gas_nx
!
!-- initializing Planck integral
     speclump = 0d0
     do ig = 1, gas_ng
!-- finding lumpable groups
        if(gas_cap(ig,i,j,k)*dx(i)*thelp>=prt_taulump) then
!-- summing lumpable Planck function integrals
           speclump = speclump + gas_siggrey(i,j,k)*gas_emitprob(ig,i,j,k)/&
                gas_cap(ig,i,j,k)
        endif
     enddo !ig
!-- lumping opacity
     do ig = 1, gas_ng
        if(gas_cap(ig,i,j,k)*dx(i)*thelp>=prt_taulump) then
!
!-- obtaining spectral weight
           specval = gas_siggrey(i,j,k)*gas_emitprob(ig,i,j,k)/&
                gas_cap(ig,i,j,k)
!
!
!-- calculating left leakage opacity
           if(i==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i-1,j,k)+ &
                 gas_sig(i-1,j,k))*dx(i-1)*thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dx(i)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(1,i,j,k)=gas_opacleak(1,i,j,k)+(specval/speclump)*&
                   1.5d0*ppl*(thelp*gas_xarr(i))**2/ &
                   (3d0*gas_vals2(i,j,k)%vol/pc_pi4)
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dx(i)+&
                   (gas_sig(i-1,j,k)+gas_cap(ig,i-1,j,k))*dx(i-1))*thelp
              gas_opacleak(1,i,j,k)=gas_opacleak(1,i,j,k)+(specval/speclump)*&
                   2.0d0*(thelp*gas_xarr(i))**2/ &
                   (help*3d0*gas_vals2(i,j,k)%vol/pc_pi4)
           endif
!
!
!-- calculating right leakage opacity
           if(i==gas_nx) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i+1,j,k)+ &
                 gas_sig(i+1,j,k))*dx(i+1)*thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dx(i)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(2,i,j,k)=gas_opacleak(2,i,j,k)+(specval/speclump)*&
                   1.5d0*ppr*(thelp*gas_xarr(i+1))**2/ &
                   (3d0*gas_vals2(i,j,k)%vol/pc_pi4)
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dx(i)+&
                   (gas_sig(i+1,j,k)+gas_cap(ig,i+1,j,k))*dx(i+1))*thelp
              gas_opacleak(2,i,j,k)=gas_opacleak(2,i,j,k)+(specval/speclump)*&
                   2.0d0*(thelp*gas_xarr(i+1))**2/ &
                   (help*3d0*gas_vals2(i,j,k)%vol/pc_pi4)
           endif
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k
  

end subroutine leakage_opacity
