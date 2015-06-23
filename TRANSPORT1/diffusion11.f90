pure subroutine diffusion11(ptcl,ptcl2,rndstate,edep,eraddens,totevelo,icspec,specarr,ierr)

  use randommod
  use miscmod
  use groupmod
  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens
  real*8,intent(inout) :: totevelo
  integer,intent(inout) :: icspec
  real*8,intent(inout) :: specarr(grp_ng)
  integer,intent(out) :: ierr
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
  real*8 :: dist, help

  integer,pointer :: ix, ic, ig
  integer,parameter :: iy=1, iz=1
  real*8,pointer :: r, mu, e, e0, wl
!-- statement function
  integer :: l
  real*8 :: dx,dx3
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx3(l) = grd_xarr(l+1)**3 - grd_xarr(l)**3

  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  r => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
!
!-- shortcuts
  dtinv = 1d0/tsp_dt
  tempinv = 1d0/grd_temp(ic)
  capgreyinv = max(1d0/grd_capgrey(ic),0d0) !catch nans

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
  dist = dx(ix)*thelp
  if(grd_cap(ig,ic)*dist >= prt_taulump) then
     do iig=1,grp_ng
        if(grd_cap(iig,ic)*dist >= prt_taulump .and. &
             (grd_sig(ic) + grd_cap(iig,ic))*dist >= prt_tauddmc) then
           glump=glump+1
           glumps(glump)=iig
        else
           glumps(gunlump)=iig
           gunlump=gunlump-1
        endif
     enddo
  endif
! write(0,*) ptcl2%ipart,ptcl2%istep,glump,ig,ix
!
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic))*dist < prt_tauddmc) then
     ierr = 100
     return
  endif

!
!-- only do this if needed
  if(glump>0 .and. icspec/=ic) then
     icspec = ic
     specarr = specintv(tempinv,0) !this is slow!
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
        caplump = grd_capgrey(ic)
     else
        do iig=1,glump
           iiig = glumps(iig)
           specig = specarr(iiig)
!-- emission lump
           emitlump = emitlump + specig*capgreyinv*grd_cap(iiig,ic)
!-- Planck x-section lump
           caplump = caplump + specig*grd_cap(iiig,ic)*speclump
        enddo
        emitlump = min(emitlump,1d0)
     endif
!-- leakage opacities
     opacleak = grd_opacleak(:2,ic)
!!}}}
!-- calculating unlumped values
  else
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ic)!{{{
     caplump = grd_cap(ig,ic)
!-- inward
     if(ix/=1) l = grd_icell(ix-1,iy,iz)
     if(ix==1) then
        opacleak(1) = 0d0
     elseif((grd_cap(ig,l)+ &
          grd_sig(l))*dx(ix-1)*thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(1)= 1.5d0*ppl*(thelp*grd_xarr(ix))**2/ &
             (thelp**3*dx3(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix-1))*thelp
        opacleak(1)=2.0d0*(thelp*grd_xarr(ix))**2/ &
             (mfphelp*thelp**3*dx3(ix))
     endif
!
!-- outward
     if(ix==grd_nx) then
        lhelp = .true.
     else
        l = grd_icell(ix+1,iy,iz)
        lhelp = (grd_cap(ig,l)+ &
           grd_sig(l))*dx(ix+1)*thelp<prt_tauddmc
     endif
!
     if(lhelp) then
!-- DDMC interface
        mfphelp = (grd_cap(ig,ic)+grd_sig(ic))*dx(ix)*thelp
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleak(2)=1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
             (thelp**3*dx3(ix))
     else
!-- DDMC interior
        mfphelp = ((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)+&
             (grd_sig(l)+grd_cap(ig,l))*dx(ix+1))*thelp
        opacleak(2)=2.0d0*(thelp*grd_xarr(ix+1))**2/ &
             (mfphelp*thelp**3*dx3(ix))
     endif!}}}
  endif
!
!-------------------------------------------------------------
!

!-- calculate time to census or event
  denom = sum(opacleak) + &
       (1d0-emitlump)*(1d0-grd_fcoef(ic))*caplump
  if(prt_isddmcanlog) then
     denom = denom+grd_fcoef(ic)*caplump
  endif

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     eraddens = e*ddmct*dtinv
  else
     edep = e*(1d0-exp(-grd_fcoef(ic) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*dx(ix)*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        eraddens  = e*ddmct*dtinv
     endif
     e = e*exp(-grd_fcoef(ic)*caplump*pc_c*ddmct)
!!}}}
  endif


!-- updating particle time
  ptcl%t = ptcl%t+ddmct


!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
     return
  endif


!-- otherwise, perform event
  call rnd_r(r1,rndstate)
  help = 1d0/denom

!-- leak probability
  probleak = opacleak*help

!-- absorption probability
  if(prt_isddmcanlog) then
     pa = grd_fcoef(ic)*caplump*help
  else
     pa = 0d0
  endif

!-- absorption sample
  if(r1<pa) then
     ptcl2%isvacant = .true.
     ptcl2%done = .true.
     edep = e

!-- left leakage sample
  elseif (r1>=pa .and. r1<pa+probleak(1)) then
!{{{
!-- checking if at inner bound
     if(ix==1) then
!       stop 'diffusion11: non-physical inward leakage'
        ierr = 101
        return

!-- sample adjacent group (assumes aligned ig bounds)
     else

        l = grd_icell(ix-1,iy,iz)
        if(speclump<=0d0) then
           iiig = ig
        else
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(1)
           do iig=1,glump
              iiig = glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,l)+ &
                   grd_sig(l))*dx(ix-1)*thelp<prt_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ic)+grd_sig(ic))*dx(ix)*thelp
                 ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppl*(thelp*grd_xarr(ix))**2/ &
                      (thelp**3*dx3(ix))
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dx(ix)+&
                      (grd_sig(l)+grd_cap(iiig,l))*dx(ix-1))*thelp
                 resopacleak = 2.0d0*(thelp*grd_xarr(ix))**2/ &
                      (mfphelp*thelp**3*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif
!
!-- sampling wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!
!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix-1) &
             *thelp < prt_tauddmc) then
           ptcl2%itype = 1
!
!-- location set right bound of left cell
           r = grd_xarr(ix)
!
!-- particle angle sampled from isotropic b.c. inward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
!
!-- doppler and aberration corrections
           if(grd_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
              help = 1d0/(1.0-r*mu*cinv)
              totevelo = totevelo+e*(1d0 - help)
!
              e = e*help
              e0 = e0*help
              wl = wl*(1.0-r*mu*cinv)
           endif
        endif
!
!-- update particle
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
        ig = iiig

     endif!}}}


!-- right leakage sample
  elseif (r1>=pa+probleak(1) .and. r1<pa+sum(probleak)) then
!!{{{
!-- checking if at outer bound
     if(ix==grd_nx) then
        ptcl2%isvacant = .true.!{{{
        ptcl2%done = .true.
        ptcl2%lflux = .true.
!-- outbound luminosity tally
        call rnd_r(r1,rndstate)
        call rnd_r(r2,rndstate)
        mu = max(r1,r2)
        if(speclump<=0d0) then
           call rnd_r(r1,rndstate)
           wl=1d0/(r1*grp_wlinv(ig+1) + (1d0-r1)*grp_wlinv(ig))
        else
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig=glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic))*dx(ix)*thelp
              ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
                   (thelp**3*dx3(ix))
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
           call rnd_r(r1,rndstate)
           wl=1d0/(r1*grp_wlinv(iiig+1) + (1d0-r1)*grp_wlinv(iiig))
        endif
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
           help = 1d0+mu*grd_xarr(grd_nx+1)*cinv
!-- velocity effects accounting
           totevelo = totevelo+e*(1d0 - help)
           wl = wl/help
           e = e*help
           e0 = e0*help
        endif
        ig = binsrch(wl,flx_wl,flx_ng+1,.false.)
!
!!}}}
     else

        l = grd_icell(ix+1,iy,iz)
!-- sample adjacent group (assumes aligned ig bounds)
        if(speclump<=0d0) then
           iiig = ig
        else
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig = glumps(iig)
              specig = specarr(iiig)
!-- calculating resolved leakage opacities
              if((grd_cap(iiig,l)+ &
                   grd_sig(l))*dx(ix+1)*thelp<prt_tauddmc) then
!-- DDMC interface
                 mfphelp = (grd_cap(iiig,ic)+grd_sig(ic))*dx(ix)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleak = 1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
                      (thelp**3*dx3(ix))
              else
!-- IMC interface
                 mfphelp = ((grd_sig(ic)+grd_cap(iiig,ic))*dx(ix)+&
                      (grd_sig(l)+grd_cap(iiig,l))*dx(ix+1))*thelp
                 resopacleak = 2.0d0*(thelp*grd_xarr(ix+1))**2/ &
                      (mfphelp*thelp**3*dx3(ix))
              endif
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
        endif

!-- sampling wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))

!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix+1) &
             *thelp < prt_tauddmc) then
!
           ptcl2%itype = 1
!
!-- location set left bound of right cell
           r = grd_xarr(ix+1)
!
!-- particle angle sampled from isotropic b.c. outward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
!
!-- doppler and aberration corrections
           if(grd_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
              help = 1d0/(1.0-r*mu*cinv)
              totevelo = totevelo+e*(1d0 - help)
!
              e = e*help
              e0 = e0*help
              wl = wl*(1.0-r*mu*cinv)
           endif
        endif
!
!-- update particle
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
        ig = iiig
!
     endif!!}}}


!-- effective scattering sample
  else
!{{{
     if(glump==grp_ng) then
!       stop 'diffusion11: effective scattering with glump==ng'
        ierr = 102
        return
     endif

     call rnd_r(r1,rndstate)

     if(glump==0) then
        iiig = emitgroup(r1,ic)
        if(iiig>grp_ng) then
           ierr = 103
           return
        endif
     else
        denom2 = 1d0-emitlump
        denom2 = 1d0/denom2
        denom3 = 0d0
        do iig=glump+1,grp_ng
           iiig=glumps(iig)
           if(icspec==ic) then
              help = specarr(iiig)*grd_cap(iiig,ic)*capgreyinv
           else
              help = specint0(tempinv,iiig)*grd_cap(iiig,ic)*capgreyinv
           endif
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
!
     ig = iiig
     call rnd_r(r1,rndstate)
     wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))

     if((grd_sig(ic)+grd_cap(ig,ic))*dx(ix) &
          *thelp >= prt_tauddmc) then
        ptcl2%itype = 2
     else
        ptcl2%itype = 1
!-- direction sampled isotropically           
        call rnd_r(r1,rndstate)
        mu = 1.0-2.0*r1
!-- position sampled uniformly
        call rnd_r(r1,rndstate)
        r = (r1*grd_xarr(ix+1)**3 + (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
        r = min(r,grd_xarr(ix+1))
        r = max(r,grd_xarr(ix))
!-- doppler and aberration corrections
        if(grd_isvelocity) then
           mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-r*mu*cinv)
           totevelo = totevelo+e*(1d0 - help)
!
           e = e*help
           e0 = e0*help
           wl = wl*(1.0-r*mu*cinv)
        endif
     endif
!}}}
  endif

end subroutine diffusion11
