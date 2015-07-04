subroutine leakage_opacity2

  use miscmod
  use groupmod
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
  integer :: icnb(4) !neighbor cells
  real*8 :: thelp, dist, help
  real*8 :: speclump, specval
  real*8 :: specarr(grp_ng)
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

!
!-- calculating leakage opacities
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
!
!-- work distribution
     if(l<grd_idd1) cycle
     if(l>grd_idd1+grd_ndd-1) cycle
!
!-- zero
     grd_opacleak(:,l) = 0d0
!
!-- neighbors
     icnb(1) = grd_icell(max(i-1,1),j,k)      !left neighbor
     icnb(2) = grd_icell(min(i+1,grd_nx),j,k) !right neighbor
     icnb(3) = grd_icell(i,max(j-1,1),k)      !left neighbor
     icnb(4) = grd_icell(i,min(j+1,grd_ny),k) !right neighbor
!
!-- distance
     dist = min(dx(i),dy(j))*thelp
!
!-- initializing Planck integral vectorized
     specarr = specintv(1d0/grd_temp(l),0)
     speclump = sum(specarr, grd_cap(:,l)*dist>=prt_taulump .and. &
       (grd_sig(l) + grd_cap(:,l))*dist >= prt_tauddmc)
     if(speclump>0d0) then
        speclump = 1d0/speclump
     else
        speclump = 0d0
     endif
!-- lumping opacity
     do ig=1,grp_ng
        if(grd_cap(ig,l)*dist < prt_taulump) cycle
        if((grd_sig(l) + grd_cap(ig,l))*dist < prt_tauddmc) cycle
!
!-- obtaining spectral weight
        specval = specarr(ig)
!
!-- calculating inward leakage opacity
        if(i==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(1))+ &
              grd_sig(icnb(1)))*min(dx(i-1),dy(j))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(1,l)=grd_opacleak(1,l)+(specval*speclump)*&
                ppl*(thelp*grd_xarr(i))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(1))+grd_cap(ig,icnb(1)))*dx(i-1))*thelp
           grd_opacleak(1,l)=grd_opacleak(1,l)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating outward leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(2))+ &
                grd_sig(icnb(2)))*min(dx(i+1),dy(j)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(2,l)=grd_opacleak(2,l)+(specval*speclump)*&
                ppr*(thelp*grd_xarr(i+1))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(2))+grd_cap(ig,icnb(2)))*dx(i+1))*thelp
           grd_opacleak(2,l)=grd_opacleak(2,l)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i+1))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating downward leakage opacity
        if(j==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(3))+ &
                grd_sig(icnb(3)))*min(dx(i),dy(j-1)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dy(j)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(3,l)=grd_opacleak(3,l)+(specval*speclump)*&
                0.5d0*ppl/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dy(j)+&
                (grd_sig(icnb(3))+grd_cap(ig,icnb(3)))*dy(j-1))*thelp
           grd_opacleak(3,l)=grd_opacleak(3,l)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!
!-- calculating upward leakage opacity
        if(j==grd_ny) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(4))+ &
                grd_sig(icnb(4)))*min(dx(i),dy(j+1)) * &
                thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dy(j)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(4,l)=grd_opacleak(4,l)+(specval*speclump)*&
                0.5d0*ppr/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dy(j)+&
                (grd_sig(icnb(4))+grd_cap(ig,icnb(4)))*dy(j+1))*thelp
           grd_opacleak(4,l)=grd_opacleak(4,l)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity2
