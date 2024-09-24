! © 2023. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
! Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
! Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
! National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
! The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
! copies to the public, perform publicly and display publicly, and to permit others to do so.
!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
pure subroutine transport11(ptcl,ptcl2,rndstate,edep,eraddens,eamp,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use physconstmod
  use particlemod
  use transportmod
  use inputparmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  type(packet2),target,intent(inout) :: ptcl2
  type(rnd_t),intent(inout) :: rndstate
  real*8,intent(out) :: edep, eraddens, eamp
  real*8,intent(inout) :: totevelo
  integer,intent(out) :: ierr
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c

  logical :: lout
  integer :: ixnext
  real*8 :: r1, r2, thelp,thelpinv
  real*8 :: db, dcol, dcen, dthm, ddop
  real*8 :: darr(5)
  real*8 :: elabfact
  real*8 :: xold, muold
! real*8 :: x1, x2, xx0
  real*8 :: help
!-- distance out of physical reach
  real*8 :: far
  real*8 :: emitlump

  integer,pointer :: ix,ic,ig
  integer,parameter :: iy=1, iz=1
  real*8,pointer :: x, mu, e, e0, wl, d
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

  ix => ptcl2%ix
  ic => ptcl2%ic
  ig => ptcl2%ig
  d => ptcl2%dist
  x => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl

!-- no error by default
  ierr = 0
!-- init
  edep = 0d0
  eraddens = 0d0
  eamp = 0d0
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     elabfact = 1d0/(1d0 + mu*x*cinv)
     thelp = tsp_t
  else
     elabfact = 1d0
     thelp = 1d0
  endif
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*tsp_dt*thelpinv)

!-- census distance
  dcen = abs(pc_c*(tsp_t1-ptcl%t)*thelpinv)
!
!-- boundary distances
  if(ix==1 .or. mu>=-sqrt(1d0-(grd_xarr(ix)/x)**2)) then
!-- outer boundary
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
     ixnext = ix+1
  else
!-- inner boundary
     db = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
     ixnext = ix-1
  endif
!-- sanity check
  if(db/=db) then
!    stop 'transport11: db/=db'
     ierr = 1
     return
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = -log(r1)*thelpinv/grd_sig(ic)
  else
     dthm = far
  endif
!
!-- effective collision distance
  if(grd_cap(ig,ic)<=0d0) then
     dcol = far
  elseif(trn_isimcanlog) then
!-- calculating dcol for analog MC
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/grd_cap(ig,ic)
  elseif(grd_fcoef(ic)<1d0.and.grd_fcoef(ic)>=0d0) then
     call rnd_r(r1,rndstate)
     dcol = -log(r1)*thelpinv/&
          ((1d0-grd_fcoef(ic))*grd_cap(ig,ic))
  else
     dcol = far
  endif
!
!-- Doppler shift distance
  if(grd_isvelocity .and. ig < grp_ng) then
     ddop = -pc_c*log(wl*grp_wlinv(ig+1))
!     ddop = pc_c*(grp_wl(ig+1)-wl)/wl
!     if(ddop<0d0) ddop = 0d0
  else
     ddop = far
  endif
!
!--finding minimum distance
  darr = [dcen,dcol,dthm,ddop,db]
  ptcl2%idist = minloc(darr,dim=1)
  d = minval(darr)
  if(any(darr/=darr)) then
     ierr = 3
     return
  endif
  if(d<0d0) then
     ierr = 4
     return
  endif

!-- updating position, angle
  xold = x
  x = sqrt(x**2 + d**2 + 2d0*d*x*mu)
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (xold*mu+d)/x
  endif
!-- updating energy, wavelength
  if(grd_isvelocity) then
     totevelo = totevelo + e*(1d0-exp(-d*cinv))
     e = e * exp(-d*cinv)
     wl = wl * exp(d*cinv)
!     wl = wl*(1d0+d*cinv)
  endif

!
!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d

!-- tallying energy densities
  if(trn_isimcanlog) then
!-- analog energy density
     eraddens = e * &
          d*thelp*cinv*tsp_dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ic)*grd_cap(ig,ic)*dx(ix)*thelp>1d-6) then
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic) * &
             grd_cap(ig,ic)*d*thelp)) / &
             (grd_fcoef(ic) * &
             grd_cap(ig,ic)*pc_c*tsp_dt)
     else
!-- analog energy density
        eraddens = e * &
             d*thelp*cinv*tsp_dtinv
     endif
!-- depositing nonanalog absorbed energy
     edep = e* &
          (1d0-exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
          d*thelp))
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic) * &
          d*thelp)
  endif

!
!-- updating transformation factors
  if(grd_isvelocity) then
     elabfact = 1d0/(1d0 + mu*x*cinv)
  endif

!
!-- census
  if(d == dcen) then
     ptcl2%stat = 'cens'
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
     call rnd_r(r1,rndstate)
     mu = 1d0-2d0*r1
  elseif(d==db) then
     lout = mu>=0d0.and.ix==grd_nx
     if(lout) then
!-- changing to observer frame (energy is in cmf now)
        if (grd_isvelocity) then
           help = 1d0/elabfact
           totevelo=totevelo+e*(1d0 - help)
           e = e*help
           e0 = e0*help
           wl = wl*exp(1d0-help) !elabfact
           mu = (mu+x*cinv)*elabfact
           mu = min(mu,1d0)
        endif
!-- observer time correction
        ptcl%t=ptcl%t-mu*x*thelp*cinv
!-- ending particle
        ptcl2%stat = 'flux'
        return
     endif
  endif

!
!-- Thomson scatter
  if(d == dthm) then
!
!-- effective collision
  elseif(d == dcol) then
     call rnd_r(r1,rndstate)
!-- checking if analog
     if(trn_isimcanlog.and.r1<=grd_fcoef(ic)) then
!-- effective absorption
        ptcl2%stat = 'dead'
!-- adding comoving energy to deposition energy
        edep = e
        return
     else
!-- effective scattering
!-- redistributing wavelength
        emitlump = grd_opaclump(8,ic)/grd_capgrey(ic)
        if(grp_ng==1) then
        elseif(emitlump<.99d0 .or. trn_nolumpshortcut .or. in_puretran) then
           call rnd_r(r1,rndstate)
           ig = emitgroup(r1,ic)
!-- checking if DDMC in new group
           if(((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)*thelp >= trn_tauddmc) &
              .and. .not.in_puretran) ptcl2%itype = 2
!-- don't sample, it will end up in the lump anyway
        else
           ptcl2%itype = 2
!-- always put this in the single most likely group
           ig = nint(grd_opaclump(9,ic))
        endif
!-- uniformly in new group
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
     endif

!
!-- outer radial bound
  elseif(d==db .and. ixnext>ix) then

     l = grd_icell(ix+1,iy,iz)!{{{
     if((grd_sig(l)+grd_cap(ig,l)) * &
          dx(ix+1)*thelp < trn_tauddmc &
          .or. in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix+1)
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        help = ddmc_emiss_bc(dx(ix+1)*thelp, grd_fcoef(l), &
                grd_cap(ig,l), grd_sig(l), pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5*abs(mu))) then
           ptcl2%itype = 2
!-- update
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
           x = grd_xarr(ix+1)
        endif
     endif!}}}

!
!-- inner radial bound
  elseif(d==db) then

     l = grd_icell(ix-1,iy,iz)!{{{
     if((grd_sig(l)+grd_cap(ig,l))*dx(ix-1) &
          *thelp < trn_tauddmc .or. &
          in_puretran) then
!-- IMC in adjacent cell
        x = grd_xarr(ix)
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
     else
!-- DDMC in adjacent cell
        help = ddmc_emiss_bc(dx(ix-1)*thelp, grd_fcoef(l), &
                grd_cap(ig,l), grd_sig(l), pc_dext)
!-- sampling
        call rnd_r(r1,rndstate)
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl2%itype = 2
           ix = ix-1
           ic = grd_icell(ix,iy,iz)
        else
!-- resampling direction
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
           x = grd_xarr(ix)
        endif
     endif!}}}

!
!-- Doppler shift
  elseif(d == ddop) then
     if(.not.grd_isvelocity) then
!       stop 'transport11: ddop and no velocity'
        ierr = 16
        return
     endif
!-- shifting group
     if(ig<grp_ng) then
        ig = ig+1
        wl = grp_wl(ig) !*elabfact
     else
        ierr = 18
     endif

!-- check if ddmc region
     if (((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)* &
          thelp >= trn_tauddmc) &
          .and. .not.in_puretran) then
        ptcl2%itype = 2
     endif
  else
!    stop 'transport11: invalid distance'
     ierr = 17
     return
  endif

end subroutine transport11
! vim: fdm=marker
