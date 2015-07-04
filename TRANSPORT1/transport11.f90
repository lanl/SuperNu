pure subroutine transport11(ptcl,ptcl2,rndstate,edep,eraddens,eamp,totevelo,ierr)

  use randommod
  use miscmod
  use gridmod
  use groupmod
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
!
  real*8 :: r1, r2, thelp,thelpinv
  real*8 :: db, dcol, dcen, dthm, ddop
  real*8 :: darr(5)
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: xold, P, muold
! real*8 :: x1, x2, xx0
  real*8 :: help
  real*8 :: ppl, ppr
!-- distance out of physical reach
  real*8 :: far

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
  if(grd_isvelocity) then
     siglabfact = 1d0 - mu*x*cinv
     dcollabfact = tsp_t*(1d0-mu*x*cinv)
     thelp = tsp_t
  else
     siglabfact = 1d0
     dcollabfact = 1d0
     thelp = 1d0
  endif
  thelpinv = 1d0/thelp

!-- distance longer than distance to census
  far = 2d0*abs(pc_c*tsp_dt*thelpinv)

!
!== DISTANCE CALCULATIONS
!
!-- distance to boundary = db
  if (ix == 1) then
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  elseif (mu < -sqrt(1d0-(grd_xarr(ix)/x)**2)) then
     db = abs(sqrt(grd_xarr(ix)**2-(1d0-mu**2)*x**2)+mu*x)
  else
     db = abs(sqrt(grd_xarr(ix+1)**2-(1d0-mu**2)*x**2)-mu*x)
  endif
!-- sanity check
  if(db/=db) then
!    stop 'transport11: db/=db'
     ierr = 1
     return
  endif
!
!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(grd_cap(ig,ic)>0d0) then
        call rnd_r(r1,rndstate)
        dcol = abs(log(r1)/(grd_cap(ig,ic)*dcollabfact))
     else
        dcol = far
     endif
  else
     if((1d0-grd_fcoef(ic))*grd_cap(ig,ic)>0d0) then
        call rnd_r(r1,rndstate)
        dcol = abs(log(r1)/((1d0-grd_fcoef(ic))*grd_cap(ig,ic)*dcollabfact))
     else
        dcol = far
     endif
  endif
!
!-- distance to Thomson-type collision = dthm
  if(grd_sig(ic)>0d0) then
     call rnd_r(r1,rndstate)
     dthm = abs(log(r1)/(grd_sig(ic)*dcollabfact))
  else
     dthm = far
  endif
!
!-- distance to census = dcen
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- distance to Doppler shift = ddop
  if(grd_isvelocity.and.ig<grp_ng) then
     ddop = abs(pc_c*(1d0-wl*grp_wlinv(ig+1))-x*mu)
  else
     ddop = far
  endif
!
!-- minimum distance = d
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
!
!== END OF DISTANCE CALCULATIONS
!
!-- position, angle, time update  
  xold = x
! x = sqrt((1d0-mu**2)*x**2+(d+x*mu)**2)
  x = sqrt(x**2 + d**2 + 2d0*d*x*mu)
!
  ptcl%t = ptcl%t + thelp*d*cinv
  muold = mu
  if(x==0d0) then
     mu = 1d0
  else
     mu = (xold*mu+d)/x
  endif

!-- transformation factor set
  if(grd_isvelocity) then
     elabfact = 1d0 - muold*xold*cinv
  else
     elabfact = 1d0
  endif
  !calculating energy deposition and density
  !
  if(.not.prt_isimcanlog) then
     edep = e*(1d0-exp(-grd_fcoef(ic) &
          *grd_cap(ig,ic)*siglabfact*d*thelp))*elabfact
     !--
     if(grd_fcoef(ic)*grd_cap(ig,ic)*dx(ix)*thelp>1d-6) then     
        eraddens = e* &
             (1d0-exp(-grd_fcoef(ic)*siglabfact*grd_cap(ig,ic)*d*thelp))* &
             elabfact/(grd_fcoef(ic)*siglabfact*grd_cap(ig,ic)*pc_c*tsp_dt)
     else
        eraddens = e* &
             elabfact*d*dcollabfact*cinv*tsp_dtinv
     endif
     !--
!     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic)*d*dcollabfact)
     e = e*exp(-grd_fcoef(ic)*grd_cap(ig,ic)*siglabfact*d*thelp)

  else
     !
     eraddens = e* &
          elabfact*d*dcollabfact*cinv*tsp_dtinv
  endif

!-- transformation factor reset
  if(grd_isvelocity) then
     elabfact = 1d0 - mu*x*cinv
  else
     elabfact = 1d0
  endif

  !
  if(d == ddop) then !group shift
!-- redshifting!{{{
     if(ig<grp_ng) then
        ig = ig+1
!-- lab frame wavelength
        wl = (grp_wl(ig)+1d-6*(grp_wl(ig+1)-grp_wl(ig)))*(1d0-mu*x*cinv)
     else
        call rnd_r(r1,rndstate)
        wl=1d0/(r1*grp_wlinv(grp_ng+1) + (1d0-r1)*grp_wlinv(grp_ng))  !reciprocal sample
        wl = wl*(1d0-mu*x*cinv)
     endif
!-- check if ddmc region
     if (((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)* &
          thelp >= prt_tauddmc) &
          .and. .not.in_puretran) then
        ptcl2%itype = 2
        if(grd_isvelocity) then
!-- velocity effects accounting
           totevelo=totevelo+e*x*mu*cinv
!
           e = e*(1d0-x*mu*cinv)
           e0 = e0*(1d0-x*mu*cinv)
           wl = wl/(1d0-x*mu*cinv)
        endif
     endif
!!}}}
  elseif (d == dthm) then  !physical scattering (Thomson-type)
     !!{{{
     call rnd_r(r1,rndstate)
     mu = 1d0-2d0*r1
     if(abs(mu)<0.0000001d0) then
        mu = 0.0000001d0
     endif
     if(grd_isvelocity) then
        mu = (mu+x*cinv)/(1d0+x*mu*cinv)
!-- velocity effects accounting
        help = 1d0/(1d0-mu*x*cinv)
        totevelo=totevelo+e*(1d0-elabfact*help)
!
        e = e*elabfact*help
!        e0 = e0*elabfact/(1d0-mu*x*cinv)
        wl = wl*(1d0-mu*x*cinv)/elabfact
     endif
     !
     !!}}}
  elseif(d == dcol) then  !fictitious scattering with implicit capture
     !!{{{
     call rnd_r(r1,rndstate)
     if(r1<=grd_fcoef(ic).and.prt_isimcanlog) then
        ptcl2%isvacant = .true.
        ptcl2%done = .true.
        edep = e*elabfact
!-- velocity effects accounting
        totevelo = totevelo+e*(1d0-elabfact)
!
     else
        call rnd_r(r1,rndstate)
        mu = 1d0-2d0*r1
        if(abs(mu)<0.0000001d0) then
           mu = 0.0000001d0
        endif
        if(grd_isvelocity) then
           mu = (mu+x*cinv)/(1d0+x*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1d0-mu*x*cinv)
           totevelo = totevelo+e*(1d0-elabfact*help)
!
           e = e*elabfact*help
!           wl = wl*(1d0-mu*x*cinv)/elabfact
        endif
!
!-- sample wavelength
        call rnd_r(r1,rndstate)
        if(grp_ng>1) ig = emitgroup(r1,ic)
! sampling comoving wavelength in group
        call rnd_r(r1,rndstate)
        wl = 1d0/((1d0-r1)*grp_wlinv(ig)+r1*grp_wlinv(ig+1))
!-- converting comoving wavelength to lab frame wavelength
        if(grd_isvelocity) wl = wl*(1d0-x*mu*cinv)
!
        if (((grd_sig(ic)+grd_cap(ig,ic))*dx(ix)*thelp >= prt_tauddmc) &
             .and. .not.in_puretran) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo = totevelo+e*x*mu*cinv
!
              e = e*(1d0-x*mu*cinv)
              e0 = e0*(1d0-x*mu*cinv)
           endif
        else
           ptcl%icorig = ic
        endif
     endif
     !!}}}
  elseif(d==db .and. mu>=0d0) then   !------ outwards boundary crossing ----
     if(ix/=grd_nx) l = grd_icell(ix+1,iy,iz)!{{{
     if(ix == grd_nx) then
!           if(ig/=1) then
        ptcl2%isvacant = .true.
        ptcl2%done = .true.
        ptcl2%lflux = .true.
!
!-- Checking if DDMC region right
     elseif (((grd_sig(l)+grd_cap(ig,l))*dx(ix+1) &
          *thelp >= prt_tauddmc) &
          .and. .not.in_puretran) then
        call rnd_r(r1,rndstate)
        if(grd_isvelocity) then
           mu = (mu-x*cinv)/(1d0-x*mu*cinv)
        endif
        help = (grd_cap(ig,l)+grd_sig(l))*dx(ix+1)*thelp
        ppl = 4d0/(3d0*help+6d0*pc_dext)
        P = ppl*(1d0+1.5*abs(mu))
!--
        if (r1 < P) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo=totevelo+e*(1d0-elabfact)
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
!-- update
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = -max(r1,r2)
           if(grd_isvelocity) then
              mu = (mu+x*cinv)/(1d0+x*mu*cinv)
           endif
           x = grd_xarr(ix+1)
        endif
     ! End of check
     else
        x = grd_xarr(ix+1)
        ix = ix+1
        ic = grd_icell(ix,iy,iz)
     endif!}}}
  elseif(d==db .and. mu<0d0) then   !------ inwards boundary crossing ----
     if(ix/=1) l = grd_icell(ix-1,iy,iz)!{{{
     if(ix==1) then
        l = grd_icell(ix+1,iy,iz)
        if (((grd_sig(l)+grd_cap(ig,l))*dx(ix+1) &
             *thelp >= prt_tauddmc) &
             .and. .not.in_puretran) then
           call rnd_r(r1,rndstate)
           if(grd_isvelocity) then
              mu = (mu-x*cinv)/(1d0-x*mu*cinv)
           endif
           help = (grd_cap(ig,l)+grd_sig(l))*dx(ix+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           P = ppl*(1d0+1.5*abs(mu))
           if (r1 < P) then
              ptcl2%itype = 2
              if(grd_isvelocity) then
!-- velocity effects accounting
                 totevelo=totevelo+e*(1d0-elabfact)
!
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
!-- update
              ix = ix+1
              ic = grd_icell(ix,iy,iz)
           else
              call rnd_r(r1,rndstate)
              call rnd_r(r2,rndstate)
              mu = -max(r1,r2)
              if(grd_isvelocity) then
                 mu = (mu+x*cinv)/(1d0+x*mu*cinv)
              endif
              x = grd_xarr(ix+1)
           endif
        else
           x = grd_xarr(ix+1)
           ix = ix+1
           ic = grd_icell(ix,iy,iz)
        endif
     elseif (((grd_sig(l)+grd_cap(ig,l))*dx(ix-1) &
          *thelp >= prt_tauddmc) &
          .and. .not.in_puretran) then
        call rnd_r(r1,rndstate)
!
!-- amplification factor
        if(grd_isvelocity) then
           mu = (mu-x*cinv)/(1d0-x*mu*cinv)
!
!             e0=e0*(1d0+2d0*min(0.055d0*prt_tauddmc,1d0)*x*cinv)
!             e = e*(1d0+2d0*min(0.055d0*prt_tauddmc,1d0)*x*cinv)
              if(.not.in_trn_noamp .and. mu<0d0) then
                 help = 1d0/abs(mu)
                 help = min(100d0, help) !-- truncate singularity
!
!-- velocity effects accounting
                 totevelo = totevelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv
!
!-- apply the excess (higher than factor 2d0) to the energy deposition
                 eamp = e*2d0*0.55d0*max(0d0,help-2d0)*x*cinv
!-- apply limited correction to the particle
                 help = min(2d0,help)
                 e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
                 e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*x*cinv)
              endif
!--
        endif
        help = (grd_cap(ig,l)+grd_sig(l))*dx(ix-1)*thelp
        ppr = 4d0/(3d0*help+6d0*pc_dext)
        P = ppr*(1d0+1.5*abs(mu))
!--
        if (r1 < P) then
           ptcl2%itype = 2
           if(grd_isvelocity) then
!-- velocity effects accounting
              totevelo = totevelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
!-- update
           ix = ix-1
           ic = grd_icell(ix,iy,iz)
        else
           call rnd_r(r1,rndstate)
           call rnd_r(r2,rndstate)
           mu = max(r1,r2)
           if(grd_isvelocity) then
              mu = (mu+x*cinv)/(1d0+x*mu*cinv)
           endif
           x = grd_xarr(ix)
        endif
     ! End of check
     else
        x = grd_xarr(ix)
        ix = ix-1
        ic = grd_icell(ix,iy,iz)
     endif!}}}
  elseif(d == dcen) then
     ptcl2%done = .true.
     ptcl2%lcens = .true.
!     toterad = toterad + e*elabfact
!
  endif

end subroutine transport11
