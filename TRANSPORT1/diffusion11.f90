subroutine diffusion11(ptcl,ig,isvacant,icell,specarr)

  use groupmod
  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  use totalsmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  integer,intent(inout) :: ig
  logical,intent(inout) :: isvacant
  integer,intent(inout) :: icell(3)
  real*8,intent(inout) :: specarr(grp_ng)
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  real*8,parameter :: deleff=0.38d0
!
  integer :: iig, iiig
  logical :: lhelp
  integer,external :: binsrch,emitgroup
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus, pa
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: mfphelp, ppl, ppr
  real*8 :: opacleak(2)
  real*8 :: probleak(2)
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer :: glumps(grp_ng)
  real*8 :: dtinv, tempinv, capgreyinv
  real*8 :: help

  integer,pointer :: ix
  real*8,pointer :: r, mu, e, e0, wl
!-- statement function
  integer :: l
  real*8 :: dx,dx3
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx3(l) = grd_xarr(l+1)**3 - grd_xarr(l)**3

  ix => ptcl%ix
  r => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl
!
!-- shortcuts
  dtinv = 1d0/tsp_dt
  tempinv = 1d0/grd_temp(ix,1,1)
  capgreyinv = max(1d0/grd_capgrey(ix,1,1),0d0) !catch nans

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!
!-- lump testing ---------------------------------------------
  glump = 0
  gunlump = grp_ng
  glumps = 0
!
!-- find lumpable groups
  if(grd_cap(ig,ix,1,1)*dx(ix)*thelp >= prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ix,1,1)*dx(ix)*thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
  endif
! write(0,*) prt_ipart,prt_istep,glump,ig,ix

!
!-- only do this if needed
  if(glump>0 .and. .not.all(icell==[ix,1,1])) then
     icell = [ix,1,1]
     specarr = specintv(tempinv) !this is slow!
  endif

!
  if(glump==0) then
     forall(iig=1:grp_ng) glumps(iig)=iig
  endif

!
!-- lumping
  speclump = 0d0
  do iig=1,glump
     iiig = glumps(iig)
     specig = specarr(iiig)
     speclump = speclump + specig
  enddo
  if(speclump>0d0) then
     speclump = 1d0/speclump
  else
     speclump = 0d0
  endif

!write(0,*) impi,glump,speclump
!
  emitlump = 0d0
  caplump = 0d0
!-- calculate lumped values
  if(speclump>0d0) then
     if(glump==grp_ng) then!{{{
        emitlump = 1d0
        caplump = grd_capgrey(ix,1,1)
     else
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- emission lump
           emitlump = emitlump + specig*capgreyinv*grd_cap(iiig,ix,1,1)
!-- Planck x-section lump
           caplump = caplump + specig*grd_cap(iiig,ix,1,1)*speclump
        enddo
        emitlump = min(emitlump,1d0)
     endif
!-- leakage opacities
     opacleak = grd_opacleak(:2,ix,1,1)
!!}}}
!-- calculating unlumped values
  else
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ix,1,1)!{{{
     caplump = grd_cap(ig,ix,1,1)
!-- inward
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,ix-1,1,1)+ &
          grd_sig(ix-1,1,1))*dx(ix-1)*thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ix,1,1)+grd_sig(ix,1,1))*dx(ix)*thelp
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= 1.5d0*ppl*(thelp*grd_xarr(ix))**2/ &
             (thelp**3*dx3(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,1,1)+grd_cap(ig,ix,1,1))*dx(ix)+&
             (grd_sig(ix-1,1,1)+grd_cap(ig,ix-1,1,1))*dx(ix-1))*thelp
        opacleak(1)=2.0d0*(thelp*grd_xarr(ix))**2/ &
             (mfphelp*thelp**3*dx3(ix))
     endif
!
!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        lhelp = (grd_cap(ig,ix+1,1,1)+ &
           grd_sig(ix+1,1,1))*dx(ix+1)*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ix,1,1)+grd_sig(ix,1,1))*dx(ix)*thelp
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
             (thelp**3*dx3(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ix,1,1)+grd_cap(ig,ix,1,1))*dx(ix)+&
             (grd_sig(ix+1,1,1)+grd_cap(ig,ix+1,1,1))*dx(ix+1))*thelp
        opacleak(2)=2.0d0*(thelp*grd_xarr(ix+1))**2/ &
             (mfphelp*thelp**3*dx3(ix))
     endif!}}}
  endif
!
!-------------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
       (1d0-emitlump)*(1d0-grd_fcoef(ix,1,1))*caplump
  if(prt_isddmcanlog) then
     denom = denom+grd_fcoef(ix,1,1)*caplump
  endif

  r1 = rand()
  prt_tlyrand = prt_tlyrand+1
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     grd_eraddens(ix,1,1) = grd_eraddens(ix,1,1)+e*ddmct*dtinv
  else
     grd_edep(ix,1,1) = grd_edep(ix,1,1)+e*(1d0-exp(-grd_fcoef(ix,1,1) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ix,1,1)*caplump*dx(ix)*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ix,1,1)*caplump)
        grd_eraddens(ix,1,1)= &
             grd_eraddens(ix,1,1)+e* &
             (1d0-exp(-grd_fcoef(ix,1,1)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        grd_eraddens(ix,1,1) = grd_eraddens(ix,1,1)+e*ddmct*dtinv
     endif
     e = e*exp(-grd_fcoef(ix,1,1)*caplump*pc_c*ddmct)
!!}}}
  endif


!-- updating particle time
  ptcl%t = ptcl%t+ddmct


!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     prt_done = .true.
     grd_numcensus(ix,1,1)=grd_numcensus(ix,1,1)+1
     return
  endif


!-- otherwise, perform event
  r1 = rand()
  prt_tlyrand = prt_tlyrand+1
  help = 1d0/denom

!-- leak probability
  probleak = opacleak*help

!-- absorption probability
  if(prt_isddmcanlog) then
     pa = grd_fcoef(ix,1,1)*caplump*help
  else
     pa = 0d0
  endif

!-- absorption sample
  if(r1<pa) then
     isvacant = .true.
     prt_done = .true.
     grd_edep(ix,1,1) = grd_edep(ix,1,1)+e

!-- left leakage sample
  elseif (r1>=pa .and. r1<pa+probleak(1)) then
!{{{
!-- checking if at inner bound
     if (ix == 1) then
        stop 'diffusion11: non-physical inward leakage'

!-- sample adjacent group (assumes aligned ig bounds)
     else

        if(speclump<=0d0) then
           iiig = ig
        else
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           denom2 = 0d0
           help = 1d0/opacleak(1)
           do iig=1,glump
              iiig = glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,ix-1,1,1)+ &
                   grd_sig(ix-1,1,1))*dx(ix-1)*thelp<prt_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ix,1,1)+grd_sig(ix,1,1))*dx(ix)*thelp
                 ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppl*(thelp*grd_xarr(ix))**2/ &
                      (thelp**3*dx3(ix))
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ix,1,1)+grd_cap(iiig,ix,1,1))*dx(ix)+&
                      (grd_sig(ix-1,1,1)+grd_cap(iiig,ix-1,1,1))*dx(ix-1))*thelp
                 resopacleak = 2.0d0*(thelp*grd_xarr(ix))**2/ &
                      (mfphelp*thelp**3*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif


        if((grd_sig(ix-1,1,1)+grd_cap(iiig,ix-1,1,1))*dx(ix-1) &
             *thelp >= prt_tauddmc) then
           ix = ix-1
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
           ig = iiig
        else
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!
!-- method changed to IMC
           ptcl%itype = 1
           grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
!
!-- location set right bound of left cell
           r = grd_xarr(ix)
!-- current particle cell set to 1 left
           ix = ix-1
!
!-- particl angle sampled from isotropic b.c. inward
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = -max(r1,r2)
!
!-- doppler and aberration corrections
           if(grd_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
              help = 1d0/(1.0-r*mu*cinv)
              tot_evelo=tot_evelo+e*(1d0 - help)
!
              e = e*help
              e0 = e0*help
              wl = wl*(1.0-r*mu*cinv)
           endif
!
!-- group reset
           ig = iiig
!
        endif

     endif!}}}


!-- right leakage sample
  elseif (r1>=pa+probleak(1) .and. r1<pa+sum(probleak)) then
!!{{{
!-- checking if at outer bound
     if (ix == grd_nx) then
        isvacant = .true.
        prt_done = .true.
        tot_eout = tot_eout+e
!-- outbound luminosity tally
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        r2 = rand()
        prt_tlyrand = prt_tlyrand+1
        mu = max(r1,r2)
        if(speclump<=0d0) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl=1d0/(r1*grp_wlinv(ig+1) + (1d0-r1)*grp_wlinv(ig))
!-- changing from comoving frame to observer frame
           if(grd_isvelocity) then
              help = 1d0+mu*grd_xarr(grd_nx+1)*cinv
              wl = wl/help
           else
              help = 1d0
           endif
           iiig = binsrch(wl,flx_wl,flx_ng+1,0)
           if(iiig>flx_ng.or.iiig<1) then
              if(iiig>flx_ng) then
                 iiig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 iiig=1
                 wl=flx_wl(1)
              endif
           endif
           flx_luminos(iiig,1,1) = flx_luminos(iiig,1,1) + e*dtinv*help
           flx_lumdev(iiig,1,1) = flx_lumdev(iiig,1,1) + (e*dtinv*help)**2
           flx_lumnum(iiig,1,1) = flx_lumnum(iiig,1,1) + 1
        else
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig=glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              mfphelp = (grd_cap(iiig,ix,1,1)+grd_sig(ix,1,1))*dx(ix)*thelp
              ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
                   (thelp**3*dx3(ix))
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))
!-- changing from comoving frame to observer frame
           if(grd_isvelocity) then
              help = 1d0+mu*grd_xarr(grd_nx+1)*cinv
              wl = wl/help
           else
              help = 1d0
           endif
!-- obtaining lab frame flux group
           iiig = binsrch(wl,flx_wl,flx_ng+1,0)
           if(iiig>flx_ng.or.iiig<1) then
              if(iiig>flx_ng) then
                 iiig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 iiig=1
                 wl=flx_wl(1)
              endif
           endif
           flx_luminos(iiig,1,1) = flx_luminos(iiig,1,1) + e*dtinv*help
           flx_lumdev(iiig,1,1) = flx_lumdev(iiig,1,1) + (e*dtinv*help)**2
           flx_lumnum(iiig,1,1) = flx_lumnum(iiig,1,1) + 1
        endif
!
!
     else
!
!-- sample adjacent group (assumes aligned ig bounds)
        if(speclump<=0d0) then
           iiig = ig
        else
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig = glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,ix+1,1,1)+ &
                   grd_sig(ix+1,1,1))*dx(ix+1)*thelp<prt_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ix,1,1)+grd_sig(ix,1,1))*dx(ix)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
                      (thelp**3*dx3(ix))
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ix,1,1)+grd_cap(iiig,ix,1,1))*dx(ix)+&
                      (grd_sig(ix+1,1,1)+grd_cap(iiig,ix+1,1,1))*dx(ix+1))*thelp
                 resopacleak = 2.0d0*(thelp*grd_xarr(ix+1))**2/ &
                      (mfphelp*thelp**3*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif


        if((grd_sig(ix+1,1,1)+grd_cap(iiig,ix+1,1,1))*dx(ix+1) &
             *thelp >= prt_tauddmc) then
!
           ix = ix+1
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!--
!
        else
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!
!-- method changed to IMC
           ptcl%itype = 1
           grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
!
!-- location set left bound of right cell
           r = grd_xarr(ix+1)
!-- current particle cell set 1 right
           ix = ix+1
!
!--  particl angle sampled from isotropic b.c. outward
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           if(grd_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
              help = 1d0/(1.0-r*mu*cinv)
              tot_evelo=tot_evelo+e*(1d0 - help)
!
              e = e*help
              e0 = e0*help
              wl = wl*(1.0-r*mu*cinv)
           endif
!
!-- group reset
           ig = iiig
!
        endif
     endif!!}}}


!-- effective scattering sample
  else
!{{{
     if(glump==grp_ng) stop 'diffusion11: effective scattering with glump==ng'

     r1 = rand()
     prt_tlyrand = prt_tlyrand+1

     if(glump==0) then
        iiig = emitgroup(r1,ix,1,1)
     else
        denom2 = 1d0-emitlump
        denom2 = 1d0/denom2
        denom3 = 0d0
        do iig=glump+1,grp_ng
           iiig=glumps(iig)
           if(all(icell==[ix,1,1])) then
              help = specarr(iiig)*grd_cap(iiig,ix,1,1)*capgreyinv
           else
              help = specint0(tempinv,iiig)*grd_cap(iiig,ix,1,1)*capgreyinv
           endif
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
!
     ig = iiig
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))

     if((grd_sig(ix,1,1)+grd_cap(ig,ix,1,1))*dx(ix) &
          *thelp >= prt_tauddmc) then
        ptcl%itype = 2
     else
        ptcl%itype = 1
        grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
!-- direction sampled isotropically           
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        mu = 1.0-2.0*r1
!-- position sampled uniformly
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        r = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
        r = min(r,grd_xarr(ix+1))
        r = max(r,grd_xarr(ix))
!-- doppler and aberration corrections
        if(grd_isvelocity) then
           mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-r*mu*cinv)
           tot_evelo=tot_evelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1.0-r*mu*cinv)
        endif
     endif
!}}}
  endif

end subroutine diffusion11
