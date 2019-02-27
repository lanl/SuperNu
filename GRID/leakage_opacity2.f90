!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
!See LANL_COPYING and LANL_README for details of LANL copyright assertion.
subroutine leakage_opacity2

  use miscmod
  use groupmod
  use gridmod
  use timestepmod
  use transportmod
  use physconstmod
  implicit none
!##################################################
  !This subroutine computes
  !DDMC 2D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: i,j,k, ig, khelp,igemitmax
  integer :: icnb(6) !neighbor cells
  real*8 :: thelp, dist, help, emitmax
  real*8 :: speclump, caplump, doplump, specval
  real*8 :: specarr(grp_ng)
  real*8 :: ppl, ppr
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy,xm,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2-grd_xarr(l)**2
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
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
     grd_opaclump(:,l) = 0d0
!
!-- neighbors
     icnb(1) = grd_icell(max(i-1,1),j,k)      !left neighbor
     icnb(2) = grd_icell(min(i+1,grd_nx),j,k) !right neighbor
     icnb(3) = grd_icell(i,max(j-1,1),k)      !left neighbor
     icnb(4) = grd_icell(i,min(j+1,grd_ny),k) !right neighbor
!-- azimuthal neighbors are cyclic
     if(k==1) then
      icnb(5) = grd_icell(i,j,grd_nz) !left neighbor
     else
      icnb(5) = grd_icell(i,j,k-1) !left neighbor
     endif
     if(k==grd_nz) then
      icnb(6) = grd_icell(i,j,1) !right neighbor
     else
      icnb(6) = grd_icell(i,j,k+1) !right neighbor
     endif
!
!-- distance
     dist = min(dx(i),dy(j),xm(i)*dz(k))*thelp
!
!-- initializing Planck integral vectorized
     call specintv(grd_tempinv(l),grp_ng,specarr)
     speclump = sum(specarr, grd_cap(:,l)*dist>=trn_taulump .and. &
       (grd_sig(l) + grd_cap(:,l))*dist >= trn_tauddmc)
     if(speclump>0d0) then
        speclump = 1d0/speclump
     else
        speclump = 0d0
     endif
     grd_opaclump(7,l) = speclump
!
!-- caplump
     caplump = 0d0
     emitmax = 0d0
     igemitmax = 0
     do ig=1,grp_ng
        if(grd_cap(ig,l)*dist < trn_taulump) cycle
        if((grd_sig(l) + grd_cap(ig,l))*dist < trn_tauddmc) cycle
        help = specarr(ig)*grd_cap(ig,l)
        caplump = caplump + help
        if(help > emitmax) then
           emitmax = help
           igemitmax = ig
        endif
     enddo
!-- doplump
     doplump = 0d0
     if(grd_isvelocity) then
        do ig=1,grp_ng-1
           if(grd_cap(ig,l)*dist < trn_taulump) cycle
           if(grd_cap(ig+1,l)*dist >= trn_taulump) cycle
           if((grd_sig(l) + grd_cap(ig,l))*dist < trn_tauddmc) cycle
           help = dopspeccalc(grd_tempinv(l),ig) / (pc_c*tsp_t)
           doplump = doplump + help
        enddo
     endif
!-- store regrouped data
     grd_opaclump(8,l) = caplump
     grd_opaclump(9,l) = igemitmax
     grd_opaclump(10,l) = doplump
!
!-- lumping opacity
     do ig=1,grp_ng
        if(grd_cap(ig,l)*dist < trn_taulump) cycle
        if((grd_sig(l) + grd_cap(ig,l))*dist < trn_tauddmc) cycle
!
!-- obtaining spectral weight
        specval = specarr(ig)
!
!-- calculating inward leakage opacity
        if(i==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(1))+ &
              grd_sig(icnb(1)))*min(dx(i-1),dy(j), &
              xm(i-1)*dz(k))*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(1,l)=grd_opaclump(1,l)+(specval*speclump)*&
                ppl*(thelp*grd_xarr(i))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(1))+grd_cap(ig,icnb(1)))*dx(i-1))*thelp
           grd_opaclump(1,l)=grd_opaclump(1,l)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating outward leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(2))+ &
              grd_sig(icnb(2)))*min(dx(i+1),dy(j), &
              xm(i+1)*dz(k))*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(2,l)=grd_opaclump(2,l)+(specval*speclump)*&
                ppr*(thelp*grd_xarr(i+1))/(dx2(i)*thelp**2)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(2))+grd_cap(ig,icnb(2)))*dx(i+1))*thelp
           grd_opaclump(2,l)=grd_opaclump(2,l)+(specval*speclump)*&
                (4d0/3d0)*(thelp*grd_xarr(i+1))/(help*dx2(i)*thelp**2)
        endif

!
!-- calculating downward leakage opacity
        if(j==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(3))+ &
              grd_sig(icnb(3)))*min(dx(i),dy(j-1), &
              xm(i)*dz(k))*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dy(j)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(3,l)=grd_opaclump(3,l)+(specval*speclump)*&
                0.5d0*ppl/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dy(j)+&
                (grd_sig(icnb(3))+grd_cap(ig,icnb(3)))*dy(j-1))*thelp
           grd_opaclump(3,l)=grd_opaclump(3,l)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!
!-- calculating upward leakage opacity
        if(j==grd_ny) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(4))+ &
              grd_sig(icnb(4)))*min(dx(i),dy(j+1), &
              xm(i)*dz(k))*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dy(j)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(4,l)=grd_opaclump(4,l)+(specval*speclump)*&
                0.5d0*ppr/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dy(j)+&
                (grd_sig(icnb(4))+grd_cap(ig,icnb(4)))*dy(j+1))*thelp
           grd_opaclump(4,l)=grd_opaclump(4,l)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!-- 0 azimuthal leakage opacities if nz=1
        if(grd_nz==1) cycle
!
!-- calculating k->k-1 leakage opacity
        if(k==1) then
           khelp = grd_nz
        else
           khelp = k-1
        endif
        lhelp = (grd_cap(ig,icnb(5))+ &
           grd_sig(icnb(5)))*min(dx(i),dy(j), &
           xm(i)*dz(khelp))*thelp<trn_tauddmc

        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*xm(i) * &
              dz(k)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(5,l)=grd_opaclump(5,l)+(specval*speclump)*&
              0.5d0*ppl/(xm(i)*dz(k)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dz(k) + &
              (grd_sig(icnb(5))+grd_cap(ig,icnb(5)))*dz(khelp))
           grd_opaclump(5,l)=grd_opaclump(5,l)+(specval*speclump)*&
              2.0d0/(3d0*xm(i)**2*dz(k)*help*thelp**2)
        endif

!
!-- calculating k->k+1 leakage opacity
        if(k==grd_nz) then
           khelp = 1
        else
           khelp = k+1
        endif
        lhelp = (grd_cap(ig,icnb(6))+ &
           grd_sig(icnb(6)))*min(dx(i),dy(j), &
           xm(i)*dz(khelp))*thelp<trn_tauddmc

        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*xm(i) * &
              dz(k)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(6,l)=grd_opaclump(6,l)+(specval*speclump)*&
              0.5d0*ppr/(xm(i)*dz(k)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dz(k) + &
              (grd_sig(icnb(6))+grd_cap(ig,icnb(6)))*dz(khelp))
           grd_opaclump(6,l)=grd_opaclump(6,l)+(specval*speclump)*&
              2.0d0/(3d0*xm(i)**2*dz(k)*help*thelp**2)
        endif

     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity2
! vim: fdm=marker
