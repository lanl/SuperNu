subroutine leakage_opacity1

  use gridmod
  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC 1D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: i,j,k, ig
  real*8 :: help
  real*8 :: thelp, ppl, ppr, specval, speclump
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
!
!-- setting vel-space helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- init (necessary for domain decomposition)
  grd_opacleak = 0d0

!
!-- calculating leakage opacities
  do k = 1, grd_nz
  do j = 1, grd_ny
  do i = 1, grd_nx
!
!-- initializing Planck integral
     speclump = 0d0
     do ig = 1, gas_ng
!-- finding lumpable groups
        if(grd_cap(ig,i,j,k)*dx(i)*thelp>=prt_taulump) then
!-- summing lumpable Planck function integrals
           speclump = speclump + grd_siggrey(i,j,k)*grd_emitprob(ig,i,j,k)/&
                grd_cap(ig,i,j,k)
        endif
     enddo !ig
!-- lumping opacity
     do ig = 1, gas_ng
        if(grd_cap(ig,i,j,k)*dx(i)*thelp>=prt_taulump) then
!
!-- obtaining spectral weight
           specval = grd_siggrey(i,j,k)*grd_emitprob(ig,i,j,k)/&
                grd_cap(ig,i,j,k)
!
!
!-- calculating left leakage opacity
           if(i==1) then
              lhelp = .true.
           else
              lhelp = (grd_cap(ig,i-1,j,k)+ &
                 grd_sig(i-1,j,k))*dx(i-1)*thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval/speclump)*&
                   1.5d0*ppl*(thelp*grd_xarr(i))**2/ &
                   (3d0*grd_vol(i,j,k)/pc_pi4)
           else
!-- DDMC interior
              help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                   (grd_sig(i-1,j,k)+grd_cap(ig,i-1,j,k))*dx(i-1))*thelp
              grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval/speclump)*&
                   2.0d0*(thelp*grd_xarr(i))**2/ &
                   (help*3d0*grd_vol(i,j,k)/pc_pi4)
           endif
!
!
!-- calculating right leakage opacity
           if(i==grd_nx) then
              lhelp = .true.
           else
              lhelp = (grd_cap(ig,i+1,j,k)+ &
                 grd_sig(i+1,j,k))*dx(i+1)*thelp<prt_tauddmc
           endif
!
           if(lhelp) then
!-- DDMC interface
              help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
              ppr = 4d0/(3d0*help+6d0*pc_dext)
              grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval/speclump)*&
                   1.5d0*ppr*(thelp*grd_xarr(i+1))**2/ &
                   (3d0*grd_vol(i,j,k)/pc_pi4)
           else
!-- DDMC interior
              help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                   (grd_sig(i+1,j,k)+grd_cap(ig,i+1,j,k))*dx(i+1))*thelp
              grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval/speclump)*&
                   2.0d0*(thelp*grd_xarr(i+1))**2/ &
                   (help*3d0*grd_vol(i,j,k)/pc_pi4)
           endif
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k
  

end subroutine leakage_opacity1
