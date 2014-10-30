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
  integer :: i,j,k, ig
  real*8 :: thelp, help
  real*8 :: speclump, specval
  real*8 :: ppl, ppr
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dx2(l)= gas_xarr(l+1)**2-gas_xarr(l)**2
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
  do k=1,gas_nz
  do j=1,gas_ny
  do i=1,gas_nx
!
!-- initializing Planck integral
     speclump = 0d0
     do ig = 1, gas_ng
!-- finding lumpable groups
        if(gas_cap(ig,i,j,k)*min(dx(i),dy(j))*thelp >= &
             prt_taulump) then
!-- summing lumpable Planck function integrals
           speclump = speclump + gas_siggrey(i,j,k)*gas_emitprob(ig,i,j,k)/&
                gas_cap(ig,i,j,k)
        endif
     enddo !ig
!-- lumping opacity
     do ig = 1, gas_ng
        if(gas_cap(ig,i,j,k)*min(dx(i),dy(j))*thelp >= &
             prt_taulump) then
!
!-- obtaining spectral weight
           specval = gas_siggrey(i,j,k)*gas_emitprob(ig,i,j,k)/&
                gas_cap(ig,i,j,k)

!
!-- calculating inward leakage opacity
           if(i==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i-1,j,k)+ &
                 gas_sig(i-1,j,k))*min(dx(i-1),dy(j)) * &
                 thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dx(i)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(1,i,j,k)=gas_opacleak(1,i,j,k)+(specval/speclump)*&
                   ppl*(thelp*gas_xarr(i))/(dx2(i)*thelp**2)
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dx(i)+&
                   (gas_sig(i-1,j,k)+gas_cap(ig,i-1,j,k))*dx(i-1))*thelp
              gas_opacleak(1,i,j,k)=gas_opacleak(1,i,j,k)+(specval/speclump)*&
                   (4d0/3d0)*(thelp*gas_xarr(i))/(help*dx2(i)*thelp**2)
           endif

!
!-- calculating outward leakage opacity
           if(i==gas_nx) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i+1,j,k)+ &
                   gas_sig(i+1,j,k))*min(dx(i+1),dy(j)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dx(i)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(2,i,j,k)=gas_opacleak(2,i,j,k)+(specval/speclump)*&
                   ppr*(thelp*gas_xarr(i+1))/(dx2(i)*thelp**2)
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dx(i)+&
                   (gas_sig(i+1,j,k)+gas_cap(ig,i+1,j,k))*dx(i+1))*thelp
              gas_opacleak(2,i,j,k)=gas_opacleak(2,i,j,k)+(specval/speclump)*&
                   (4d0/3d0)*(thelp*gas_xarr(i+1))/(help*dx2(i)*thelp**2)
           endif

!
!-- calculating downward leakage opacity
           if(j==1) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i,j-1,k)+ &
                   gas_sig(i,j-1,k))*min(dx(i),dy(j-1)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dy(j)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(3,i,j,k)=gas_opacleak(3,i,j,k)+(specval/speclump)*&
                   0.5d0*ppl/(thelp*dy(j))
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dy(j)+&
                   (gas_sig(i,j-1,k)+gas_cap(ig,i,j-1,k))*dy(j-1))*thelp
              gas_opacleak(3,i,j,k)=gas_opacleak(3,i,j,k)+(specval/speclump)*&
                   (2d0/3d0)/(help*dy(j)*thelp)
           endif

!
!-- calculating upward leakage opacity
           if(j==gas_ny) then
              lhelp = .true.
           else
              lhelp = (gas_cap(ig,i,j+1,k)+ &
                   gas_sig(i,j+1,k))*min(dx(i),dy(j+1)) * &
                   thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (gas_cap(ig,i,j,k)+gas_sig(i,j,k))*dy(j)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              gas_opacleak(4,i,j,k)=gas_opacleak(4,i,j,k)+(specval/speclump)*&
                   0.5d0*ppr/(thelp*dy(j))
           else
!-- DDMC interior
              help = ((gas_sig(i,j,k)+gas_cap(ig,i,j,k))*dy(j)+&
                   (gas_sig(i,j+1,k)+gas_cap(ig,i,j+1,k))*dy(j+1))*thelp
              gas_opacleak(4,i,j,k)=gas_opacleak(4,i,j,k)+(specval/speclump)*&
                   (2d0/3d0)/(help*dy(j)*thelp)
           endif
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity2
