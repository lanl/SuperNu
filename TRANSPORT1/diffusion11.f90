pure subroutine diffusion11(ptcl,ptcl2,cache,rndstate,edep,eraddens,totevelo,ierr)

  use randommod
  use miscmod
  use groupmod
  use gridmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(grp_t_cache),target,intent(inout) :: cache
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
!
  integer :: iig, iiig
  logical :: lhelp
  real*8 :: r1, r2, thelp
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus, pa
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, caplump
  real*8 :: specig
  real*8 :: mfphelp, ppl, ppr
  real*8 :: opacleak(2)
  real*8 :: probleak(2)
  real*8 :: resopacleak
  integer :: glump, gunlump
  integer,pointer :: glumps(:)
  logical,pointer :: llumps(:)
  real*8,pointer :: tempinv, capgreyinv
  real*8,pointer :: speclump
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

  tempinv => cache%tempinv
  capgreyinv => cache%capgreyinv
  speclump => cache%speclump
  glumps => cache%glumps
  llumps => cache%llumps

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0

!
!-- set expansion helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
  dist = dx(ix)*thelp

!
!-- update cache
  if(ic/=cache%ic) then
     cache%ic = ic!{{{
     cache%istat = 0 !specarr is not cached yet
     tempinv = 1d0/grd_temp(ic)
     capgreyinv = max(1d0/grd_capgrey(ic),0d0) !catch nans

!
!-- lump testing ---------------------------------------------
     glump = 0
     gunlump = grp_ng
     glumps = 0
!
!-- find lumpable groups
     speclump = grd_opaclump(0,ic)
     if(speclump==0d0) then
        glump=0
     else
        do iig=1,grp_ng
           if(grd_cap(iig,ic)*dist >= prt_taulump .and. &
                (grd_sig(ic) + grd_cap(iig,ic))*dist >= prt_tauddmc) then
              llumps(iig) = .false.
              glump=glump+1
              glumps(glump)=iig
           else
              llumps(iig) = .true.
              glumps(gunlump)=iig
              gunlump=gunlump-1
           endif
        enddo
     endif
!
!-- calculate lumped values
     if(glump==grp_ng) then
        emitlump = 1d0
        caplump = grd_capgrey(ic)
     else
!-- Planck x-section lump
        caplump = grd_opaclump(-1,ic)*speclump
        emitlump = grd_opaclump(-1,ic)*capgreyinv
        emitlump = min(emitlump,1d0)
     endif
!
!-- save
     cache%nlump = glump
     cache%emitlump = emitlump
     cache%caplump = caplump
!}}}
  endif !cache%ic /= ic

!
!-- in lump?
  if(grd_cap(ig,ic)*dist >= prt_taulump) then
     glump = cache%nlump
  else
     glump = 0
  endif
!
!-- sanity check
  if((grd_sig(ic) + grd_cap(ig,ic))*dist < prt_tauddmc) then
     ierr = 100
     return
  endif
!
!-- retrieve from cache
  if(glump>0) then
     emitlump = cache%emitlump
     caplump = cache%caplump
  else
!-- outside the lump
     emitlump = specint0(tempinv,ig)*capgreyinv*grd_cap(ig,ic)
     caplump = grd_cap(ig,ic)
  endif




  if(glump>0) then
!-- leakage opacities
     opacleak = grd_opaclump(1:2,ic)
!-- calculating unlumped values
  else
!{{{
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
  denom = 1d0/denom

  call rnd_r(r1,rndstate)
  tau = abs(log(r1)*denom*cinv)
  tcensus = tsp_t+tsp_dt-ptcl%t
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  if(prt_isddmcanlog) then
     eraddens = e*ddmct*tsp_dtinv
  else
     edep = e*(1d0-exp(-grd_fcoef(ic) &!{{{
          *caplump*pc_c*ddmct))
     if(grd_fcoef(ic)*caplump*dx(ix)*thelp>1d-6) then
        help = 1d0/(grd_fcoef(ic)*caplump)
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic)*caplump*pc_c*ddmct))* &
             help*cinv*tsp_dtinv
     else
        eraddens  = e*ddmct*tsp_dtinv
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
  if (tcensus < tau) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
     return
  endif


!-- otherwise, perform event
  call rnd_r(r1,rndstate)

!-- leak probability
  probleak = opacleak*denom

!-- absorption probability
  if(prt_isddmcanlog) then
     pa = grd_fcoef(ic)*caplump*denom
  else
     pa = 0d0
  endif

!-- update specarr cache only when necessary. this is slow
  if(r1>=pa .and. r1<pa+sum(probleak) .and. speclump>0d0 .and. &
        iand(cache%istat,2)==0) then
     cache%istat = cache%istat + 2
     call specintw(tempinv,cache%specarr,mode=0,mask=.not.llumps)
  endif

!-- absorption sample
  if(r1<pa) then
     ptcl2%idist = -1
     ptcl2%isvacant = .true.
     ptcl2%done = .true.
     edep = e

!-- left leakage sample
  elseif (r1>=pa .and. r1<pa+probleak(1)) then
     ptcl2%idist = -3
!{{{
!-- checking if at inner bound
     if(ix==1) then
!       stop 'diffusion11: non-physical inward leakage'
        ierr = 101
        return

!-- sample adjacent group (assumes aligned ig bounds)
     else

        l = grd_icell(ix-1,iy,iz)
        if(glump==0) then
           iiig = ig
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(1)
           do iig=1,glump
              iiig = glumps(iig)
              specig = cache%specarr(iiig)
              !specig = specint0(tempinv,iiig)
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
!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix-1)*thelp < prt_tauddmc) then
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set right bound of left cell
           r = grd_xarr(ix)
!-- particle angle sampled from isotropic b.c. inward
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
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
     ptcl2%idist = -4
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
        if(glump==0) then
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig=glumps(iig)
              specig = cache%specarr(iiig)
              !specig = specint0(tempinv,iiig)
!-- calculating resolved leakage opacities
              mfphelp = (grd_cap(iiig,ic)+grd_sig(ic))*dx(ix)*thelp
              ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
              resopacleak = 1.5d0*ppr*(thelp*grd_xarr(ix+1))**2/ &
                   (thelp**3*dx3(ix))
              denom2 = denom2+specig*resopacleak*speclump*help
              if(denom2>r1) exit
           enddo
           ig = iiig
        endif
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/(r1*grp_wlinv(ig+1) + (1d0-r1)*grp_wlinv(ig))
!-- changing from comoving frame to observer frame
        if(grd_isvelocity) then
           help = 1d0+mu*grd_xarr(grd_nx+1)*cinv
!-- velocity effects accounting
           totevelo = totevelo+e*(1d0 - help)
           wl = wl/help
           e = e*help
           e0 = e0*help
           ig = binsrch(wl,flx_wl,flx_ng+1,.false.)
        endif
!
!!}}}
     else

        l = grd_icell(ix+1,iy,iz)
!-- sample adjacent group (assumes aligned ig bounds)
        if(glump==0) then
           iiig = ig
        else
!-- sample group
           call rnd_r(r1,rndstate)
           denom2 = 0d0
           help = 1d0/opacleak(2)
           do iig=1,glump
              iiig = glumps(iig)
              specig = cache%specarr(iiig)
              !specig = specint0(tempinv,iiig)
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

!-- method changes to IMC
        if((grd_sig(l)+grd_cap(iiig,l))*dx(ix+1)*thelp < prt_tauddmc) then
!
           ptcl2%itype = 1
!-- sampling wavelength
           call rnd_r(r1,rndstate)
           wl = 1d0/(r1*grp_wlinv(iiig+1)+(1d0-r1)*grp_wlinv(iiig))
!-- location set left bound of right cell
           r = grd_xarr(ix+1)
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
     ptcl2%idist = -2
!{{{
     if(glump==grp_ng) then
!       stop 'diffusion11: effective scattering with glump==ng'
        ierr = 102
        return
     endif

     if(glump==0) then
!-- sample group
        if(cache%emitlump<.99d0 .or. trn_nolumpshortcut) then
           call rnd_r(r1,rndstate)
           iiig = emitgroup(r1,ic)
!-- don't sample, it will end up in the lump anyway
        else
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(-2,ic))
           return
        endif
     else
!
!-- update specarr cache. this is slow
        if(iand(cache%istat,1)==0) then
           cache%istat = cache%istat + 1
           call specintw(tempinv,cache%specarr,mode=0,mask=llumps)
        endif
!
        call rnd_r(r1,rndstate)
        denom2 = 1d0/(1d0-emitlump)
        denom3 = 0d0
        do iig=grp_ng,glump+1,-1
           iiig = glumps(iig)
           help = cache%specarr(iiig)*grd_cap(iiig,ic)*capgreyinv
           denom3 = denom3 + help*denom2
           if(denom3>r1) exit
        enddo
     endif
     ig = iiig

     if((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)*thelp < prt_tauddmc) then
        ptcl2%itype = 1
!-- sample wavelength
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig) + r1*grp_wlinv(ig+1))
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
