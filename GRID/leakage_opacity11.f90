!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine leakage_opacity11

  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use transportmod
  use physconstmod
  implicit none
!##################################################
  !This subroutine computes
  !DDMC 1D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: i,j,k, ig,igemitmax
  real*8 :: thelp, dist, help, emitmax
  real*8 :: speclump, caplump, doplump, specval
  real*8 :: specarr(grp_ng)
  real*8 :: ppl, ppr
  integer :: icnb(2) !neighbor cell pointers
!-- statement functions
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
!
!-- distance
     dist = dx(i)*thelp
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
              grd_sig(icnb(1)))*dx(i-1)*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(1,l)=grd_opaclump(1,l)+(specval*speclump)*&
                1.5d0*ppl*(thelp*grd_xarr(i))**2/ &
                (3d0*grd_vol(l)/pc_pi4)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(1))+grd_cap(ig,icnb(1)))*dx(i-1))*thelp
           grd_opaclump(1,l)=grd_opaclump(1,l)+(specval*speclump)*&
                2.0d0*(thelp*grd_xarr(i))**2/ &
                (help*3d0*grd_vol(l)/pc_pi4)
        endif

!
!-- calculating outward leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,icnb(2))+ &
              grd_sig(icnb(2)))*dx(i+1)*thelp<trn_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,l)+grd_sig(l))*dx(i)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           grd_opaclump(2,l)=grd_opaclump(2,l)+(specval*speclump)*&
                1.5d0*ppr*(thelp*grd_xarr(i+1))**2/ &
                (3d0*grd_vol(l)/pc_pi4)
        else
!-- DDMC interior
           help = ((grd_sig(l)+grd_cap(ig,l))*dx(i)+&
                (grd_sig(icnb(2))+grd_cap(ig,icnb(2)))*dx(i+1))*thelp
           grd_opaclump(2,l)=grd_opaclump(2,l)+(specval*speclump)*&
                2.0d0*(thelp*grd_xarr(i+1))**2/ &
                (help*3d0*grd_vol(l)/pc_pi4)
        endif
     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity11
! vim: fdm=marker
