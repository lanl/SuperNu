subroutine leakage_opacity2

  use miscmod
  use gridmod
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
! real*8 :: x1, x2
  real*8 :: speclump, specval
  real*8 :: specarr(grd_ng)
  real*8 :: ppl, ppr
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2-grd_xarr(l)**2
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
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
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
!!
!!-- initializing Planck integral
!     speclump = 0d0
!     do ig=1,grd_ng
!!-- finding lumpable groups
!        if(grd_cap(ig,i,j,k)*min(dx(i),dy(j))*thelp < prt_taulump) cycle
!!-- summing lumpable Planck function integrals
!        x1 = pc_h*pc_c/(grd_wl(ig+1)*pc_kb*grd_temp(i,j,k))
!        x2 = pc_h*pc_c/(grd_wl(ig)*pc_kb*grd_temp(i,j,k))
!        specarr(ig) = specint(x1,x2,3,10)
!        speclump = speclump + specarr(ig)
!     enddo !ig
!     speclump = 1d0/speclump
!
!-- initializing Planck integral vectorized
     specarr = specint3(pc_h*pc_c/(grd_wl*pc_kb*grd_temp(i,j,k)),grd_ng+1)
     help = min(dx(i),dy(j))*thelp
     speclump = 1d0/sum(specarr,grd_cap(:,i,j,k)*help>=prt_taulump)
!-- lumping opacity
     do ig=1,grd_ng
        if(grd_cap(ig,i,j,k)*min(dx(i),dy(j))*thelp < prt_taulump) cycle
!
!-- obtaining spectral weight
        specval = specarr(ig)
!
!-- calculating inward leakage opacity
        if(i==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i-1,j,k)+ &
              grd_sig(i-1,j,k))*min(dx(i-1),dy(j))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                ppl*(thelp*grd_xarr(i))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i-1,j,k)+grd_cap(ig,i-1,j,k))*dx(i-1))*thelp
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating outward leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i+1,j,k)+ &
                grd_sig(i+1,j,k))*min(dx(i+1),dy(j)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                ppr*(thelp*grd_xarr(i+1))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i+1,j,k)+grd_cap(ig,i+1,j,k))*dx(i+1))*thelp
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i+1))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating downward leakage opacity
        if(j==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j-1,k)+ &
                grd_sig(i,j-1,k))*min(dx(i),dy(j-1)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dy(j)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                0.5d0*ppl/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dy(j)+&
                (grd_sig(i,j-1,k)+grd_cap(ig,i,j-1,k))*dy(j-1))*thelp
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!
!-- calculating upward leakage opacity
        if(j==grd_ny) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j+1,k)+ &
                grd_sig(i,j+1,k))*min(dx(i),dy(j+1)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dy(j)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                0.5d0*ppr/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dy(j)+&
                (grd_sig(i,j+1,k)+grd_cap(ig,i,j+1,k))*dy(j+1))*thelp
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity2
