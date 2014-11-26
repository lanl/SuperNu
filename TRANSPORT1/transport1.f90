subroutine transport1(ptcl,isvacant)

  use gridmod
  use totalsmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
  use fluxmod
  implicit none
!
  type(packet),target,intent(inout) :: ptcl
  logical,intent(inout) :: isvacant
!##################################################
!This subroutine passes particle parameters as input and modifies
!them through one IMC transport event (Fleck&Cummings, 1971).  If
!the puretran boolean is set to false, this routine couples to the
!analogous DDMC diffusion routine through the advance.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
!
  integer :: ig, iig, g, binsrch
  real*8 :: r1, r2, thelp,thelpinv
  real*8 :: db, dcol, dcen, dthm, ddop, d
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, P, denom2, told, muold
! real*8 :: x1, x2, xx0
  real*8 :: dtinv
  real*8 :: help
  real*8 :: ppl, ppr

  integer,pointer :: ix
  real*8,pointer :: r, mu, e, e0, wl
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

  ix => ptcl%ix
  r => ptcl%x
  mu => ptcl%mu
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl
!
!-- shortcut
  dtinv = 1d0/tsp_dt

  if(grd_isvelocity) then
     siglabfact = 1.0d0 - mu*r*cinv
     dcollabfact = tsp_t*(1d0-mu*r*cinv)
     thelp = tsp_t
  else
     siglabfact = 1d0
     dcollabfact = 1d0
     thelp = 1d0
  endif
  thelpinv = 1d0/thelp

!
!-- calculating current group (rev. 120)
  if(grd_isvelocity) then
     g = binsrch(wl/(1.0d0-r*mu*cinv),grd_wl,grd_ng+1,in_ng)
  else
     g = binsrch(wl,grd_wl,grd_ng+1,in_ng)
  endif
  if(g>grd_ng.or.g<1) then
     !particle out of wlgrid energy bound
     if(g>grd_ng) then
        g=grd_ng
        if(grd_isvelocity) then
           wl=grd_wl(grd_ng+1)*(1.0d0-r*mu*cinv)
        else
           wl=grd_wl(grd_ng+1)
        endif
     elseif(g<1) then
        g=1
        if(grd_isvelocity) then
           wl=grd_wl(1)*(1.0d0-r*mu*cinv)
        else
           wl=grd_wl(1)
        endif
     else
        write(*,*) 'domain leak!!'
        prt_done = .true.
        isvacant = .true.
     endif
  endif
!
!== DISTANCE CALCULATIONS
!
!-- distance to boundary = db
  if (ix == 1) then
     db = abs(sqrt(grd_xarr(ix+1)**2-(1.0-mu**2)*r**2)-mu*r)
  elseif (mu < -sqrt(1.0d0-(grd_xarr(ix)/r)**2)) then
     db = abs(sqrt(grd_xarr(ix)**2-(1.0d0-mu**2)*r**2)+mu*r)
  else
     db = abs(sqrt(grd_xarr(ix+1)**2-(1.0d0-mu**2)*r**2)-mu*r)
  endif
!-- sanity check
  if(db/=db) stop 'transport1: db/=db'
!
!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(grd_cap(g,ix,1,1)>0d0) then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        dcol = abs(log(r1)/(grd_cap(g,ix,1,1)*dcollabfact))
     else
        dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
     endif
  else
     if((1.0d0-grd_fcoef(ix,1,1))*grd_cap(g,ix,1,1)>0.0d0) then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        dcol = abs(log(r1)/((1.0d0-grd_fcoef(ix,1,1))*grd_cap(g,ix,1,1)*dcollabfact))
     else
        dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
     endif
  endif
!
!-- distance to Thomson-type collision = dthm
  if(grd_sig(ix,1,1)>0.0d0) then
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     dthm = abs(log(r1)/(grd_sig(ix,1,1)*dcollabfact))
  else
     dthm = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- distance to census = dcen
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- distance to Doppler shift = ddop
  if(grd_isvelocity.and.g<grd_ng) then
!      r1 = rand()
!      prt_tlyrand=prt_tlyrand+1
!      ddop = pc_c*tsp_t*(grd_wl(g+1)-grd_wl(g))*abs(log(r1))/(grd_wl(g)*dcollabfact)
!     wl = r1*grd_wl(g)+(1d0-r1)*grd_wl(g+1) !uniform sample
!      wl=1d0/(r1/grd_wl(g+1) + (1d0-r1)/grd_wl(g))  !reciprocal sample
!      wl=wl*(1d0-mu*r*cinv)
!      ddop = pc_c*(1d0-mu*r*cinv)*(1d0-wl/(1d0-mu*r*cinv*grd_wl(g+1)))
!     ddop = pc_c*(1d0-mu*r*cinv)*(1d0-&
!          grd_wl(g)*log(grd_wl(g+1)/grd_wl(g))/(grd_wl(g+1)-grd_wl(g)))
!     write(*,*) pc_c*(wl/grd_wl(g+1)-1d0)+r*mu
     ddop = abs(pc_c*(1d0-wl/grd_wl(g+1))-r*mu)
  else
     ddop = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- minimum distance = d
  d = min(dcol,dthm,db,dcen,ddop)
  if(d<0d0) stop 'transport1: negative distance'
!
!== END OF DISTANCE CALCULATIONS
!
!-- position, angle, time update  
  rold = r
  r = sqrt((1.0d0-mu**2)*r**2+(d+r*mu)**2)
!  r = sqrt(r**2+d**2+2d0*d*r*mu)
!
  told = ptcl%t
  ptcl%t = ptcl%t + thelp*d*cinv
  muold = mu
  mu = (rold*mu+d)/r

!-- transformation factor set
  if(grd_isvelocity) then
     elabfact = 1.0d0 - muold*rold*cinv
  else
     elabfact = 1d0
  endif
  !calculating energy deposition and density
  !
  if(.not.prt_isimcanlog) then
        grd_edep(ix,1,1)=grd_edep(ix,1,1)+e*(1.0d0-exp(-grd_fcoef(ix,1,1) &
             *grd_cap(g,ix,1,1)*siglabfact*d*thelp))*elabfact
     !--
     if(grd_fcoef(ix,1,1)*grd_cap(g,ix,1,1)*dx(ix)*thelp>1d-6) then     
        grd_eraddens(ix,1,1) = grd_eraddens(ix,1,1)+e* &
             (1.0d0-exp(-grd_fcoef(ix,1,1)*siglabfact*grd_cap(g,ix,1,1)*d*thelp))* &
             elabfact/(grd_fcoef(ix,1,1)*siglabfact*grd_cap(g,ix,1,1)*pc_c*tsp_dt)
     else
        grd_eraddens(ix,1,1) = grd_eraddens(ix,1,1)+e* &
             elabfact*d*dcollabfact*cinv*dtinv
     endif
     !--
!     e = e*exp(-grd_fcoef(ix)*grd_cap(g,ix)*d*dcollabfact)
     e = e*exp(-grd_fcoef(ix,1,1)*grd_cap(g,ix,1,1)*siglabfact*d*thelp)

  else
     !
     grd_eraddens(ix,1,1) = grd_eraddens(ix,1,1)+e* &
          elabfact*d*dcollabfact*cinv*dtinv
  endif

!-- transformation factor reset
  if(grd_isvelocity) then
     elabfact = 1.0d0 - mu*r*cinv
  else
     elabfact = 1d0
  endif

  !
  if(d == ddop) then !group shift
!     r1 = rand()!{{{
!     prt_tlyrand=prt_tlyrand+1
!-- redshifting
     if(g<grd_ng) then
        g = g+1
!-- lab frame wavelength
!     wl = r1*grd_wl(g)+(1d0-r1)*grd_wl(g+1) !uniform sample
!        wl=1d0/(r1/grd_wl(g+1) + (1d0-r1)/grd_wl(g))  !reciprocal sample
!        wl = wl*(1d0-mu*r*cinv)
        wl = (grd_wl(g)+1d-6*(grd_wl(g+1)-grd_wl(g)))*(1d0-mu*r*cinv)
     else
        r1 = rand()
        prt_tlyrand=prt_tlyrand+1
!     wl = r1*grd_wl(grd_ng)+(1d0-r1)*grd_wl(grd_ng+1) !uniform sample
        wl=1d0/(r1/grd_wl(grd_ng+1) + (1d0-r1)/grd_wl(grd_ng))  !reciprocal sample
        wl = wl*(1d0-mu*r*cinv)
!        wl = grd_wl(grd_ng+1)*(1d0-mu*r*cinv)
     endif
!-- check if ddmc region
     if (((grd_sig(ix,1,1)+grd_cap(g,ix,1,1))*dx(ix)* &
          thelp >= prt_tauddmc) &
          .and.(in_puretran.eqv..false.)) then
        ptcl%itype = 2
        grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
        if(grd_isvelocity) then
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*r*mu*cinv
!
           e = e*(1.0-r*mu*cinv)
           e0 = e0*(1.0-r*mu*cinv)
           wl = wl/(1.0-r*mu*cinv)
        endif
     else
        ptcl%itype = 1
     endif
!!}}}
  elseif (d == dthm) then  !physical scattering (Thomson-type)
     !!{{{
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu = 1.0-2.0*r1
     if(abs(mu)<0.0000001d0) then
        mu = 0.0000001d0
     endif
     if(grd_isvelocity) then
        mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
        help = 1d0/(1.0-mu*r*cinv)
        tot_evelo=tot_evelo+e*(1d0-elabfact*help)
!
        e = e*elabfact*help
!        e0 = e0*elabfact/(1.0-mu*r*cinv)
        wl = wl*(1.0-mu*r*cinv)/elabfact
     endif
     !
     !!}}}
  elseif (d == dcol) then  !fictitious scattering with implicit capture
     !!{{{
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     if(r1<=grd_fcoef(ix,1,1).and.prt_isimcanlog) then
        isvacant = .true.
        prt_done = .true.
        grd_edep(ix,1,1) = grd_edep(ix,1,1) + e*elabfact
!-- velocity effects accounting
        tot_evelo = tot_evelo+e*(1d0-elabfact)
!
     else
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        mu = 1.0-2.0*r1
        if(abs(mu)<0.0000001d0) then
           mu = 0.0000001d0
        endif
        if(grd_isvelocity) then
           mu = (mu+r*cinv)/(1.0+r*mu*cinv)
!-- velocity effects accounting
           help = 1d0/(1.0-mu*r*cinv)
           tot_evelo=tot_evelo+e*(1d0-elabfact*help)
!
           e = e*elabfact*help
!           wl = wl*(1.0-mu*r*cinv)/elabfact
           
        endif
!
        denom2 = 0.0
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        do ig = 1, grd_ng
           iig = ig
           if ((r1>=denom2).and.(r1<denom2+grd_emitprob(ig,ix,1,1))) exit
           denom2 = denom2+grd_emitprob(ig,ix,1,1)
        enddo
        g = iig
        !(rev 121): calculating radiation energy tally per group
        !grd_eraddens(ix)=grd_eraddens(ix)+e*elabfact
        !-------------------------------------------------------
        ! sampling comoving wavelength in group
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        wl = 1d0/((1d0-r1)/grd_wl(g)+r1/grd_wl(g+1))
        !wl = (1d0-r1)*grd_wl(g)+r1*grd_wl(g+1)
        !wl = 0.5d0*(grd_wl(g)+grd_wl(g+1))
        !
        ! sampling sub-group Planck function:
!         x1 = pc_h*pc_c/(grd_wl(g+1)*pc_kb*grd_temp(ix))
!         x2 = pc_h*pc_c/(grd_wl(g)*pc_kb*grd_temp(ix))
!         if (x2<pc_plkpk) then
!            bmax = x2**3/(exp(x2)-1d0)
!         elseif (x1>pc_plkpk) then
!            bmax = x1**3/(exp(x1)-1d0)
!         else
!            bmax = pc_plkpk
!         endif
!         r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!         r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!         xx0 = (1d0-r1)*x1+r1*x2
!         do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!            r1 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            r2 = rand()
!                 prt_tlyrand = prt_tlyrand+1
!            xx0 = (1d0-r1)*x1+r1*x2
!         enddo
!         wl = pc_h*pc_c/(xx0*pc_kb*grd_temp(ix))
        !
        !
        if(grd_isvelocity) then
!-- converting comoving wavelength to lab frame wavelength
           wl = wl*(1.0-r*mu*cinv)
        endif
        if (((grd_sig(ix,1,1)+grd_cap(g,ix,1,1))*dx(ix)* &
             thelp >= prt_tauddmc) &
             .and.(in_puretran.eqv..false.)) then
           ptcl%itype = 2
           grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo = tot_evelo+e*r*mu*cinv
!
              e = e*(1.0-r*mu*cinv)
              e0 = e0*(1.0-r*mu*cinv)
              wl = wl/(1.0-r*mu*cinv)
           endif
        else
           ptcl%itype = 1
        endif
     endif
     !!}}}
  elseif (d == db) then   !------boundary crossing ----
     if (mu>=0.0d0) then!{{{
        if (ix == grd_nx) then
!           if(g/=1) then
           isvacant = .true.
           prt_done = .true.
!
!-- retrieve lab frame flux group
           g = binsrch(wl,flx_wl,flx_ng+1,0)
!
!-- check group bounds
           if(g>flx_ng.or.g<1) then
              if(g>flx_ng) then
                 g=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 g=1
                 wl=flx_wl(1)
              endif
           endif
!
!-- outbound luminosity tally
           tot_eout = tot_eout+e
           flx_luminos(g,1,1) = flx_luminos(g,1,1)+e*dtinv
           flx_lumdev(g,1,1) = flx_lumdev(g,1,1)+(e*dtinv)**2
           flx_lumnum(g,1,1) = flx_lumnum(g,1,1)+1
        ! Checking if DDMC region right
        elseif (((grd_sig(ix+1,1,1)+grd_cap(g,ix+1,1,1))*dx(ix+1) &
             *thelp >= prt_tauddmc) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(grd_isvelocity) then
              mu = (mu-r*cinv)/(1.0-r*mu*cinv)
           endif
           help = (grd_cap(g,ix+1,1,1)+grd_sig(ix+1,1,1))*dx(ix+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
           P = ppl*(1.0+1.5*abs(mu))
!--
           if (r1 < P) then
              ptcl%itype = 2
              grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
              if(grd_isvelocity) then
!-- velocity effects accounting
                 tot_evelo=tot_evelo+e*(1d0-elabfact)
!
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              ix = ix+1
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = -max(r1,r2)
              if(grd_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
              r = grd_xarr(ix+1)
           endif
        ! End of check
        else
           r = grd_xarr(ix+1)
           ix = ix+1
        endif
     else
        if (ix==1) then
           if (((grd_sig(ix+1,1,1)+grd_cap(g,ix+1,1,1))*dx(ix+1) &
                *thelp >= prt_tauddmc) &
                .and.(in_puretran.eqv..false.)) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(grd_isvelocity) then
                 mu = (mu-r*cinv)/(1.0-r*mu*cinv)
              endif
              help = (grd_cap(g,ix+1,1,1)+grd_sig(ix+1,1,1))*dx(ix+1)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              P = ppl*(1.0+1.5*abs(mu))
              if (r1 < P) then
                 ptcl%itype = 2
                 grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
                 if(grd_isvelocity) then
!-- velocity effects accounting
                    tot_evelo=tot_evelo+e*(1d0-elabfact)
!
                    e = e*elabfact
                    e0 = e0*elabfact
                    wl = wl/elabfact
                 endif
                 ix = ix+1
              else
                 r1 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 r2 = rand()
                 prt_tlyrand = prt_tlyrand+1
                 mu = -max(r1,r2)
                 if(grd_isvelocity) then
                    mu = (mu+r*cinv)/(1.0+r*mu*cinv)
                 endif
                 r = grd_xarr(ix+1)
              endif
           else
              r = grd_xarr(ix+1)
              ix = ix+1
           endif
        elseif (((grd_sig(ix-1,1,1)+grd_cap(g,ix-1,1,1))*dx(ix-1) &
             *thelp >= prt_tauddmc) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(grd_isvelocity) then
              mu = (mu-r*cinv)/(1.0-r*mu*cinv)
!-- amplification
!
!              e0=e0*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
!              e = e*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
              if(mu<0d0) then
!-- velocity effects accounting
                 help = 1d0/abs(mu)
                 tot_evelo = tot_evelo-e*2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv
!
                 e0 = e0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
                 e = e*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
              endif
!--
           endif
           help = (grd_cap(g,ix-1,1,1)+grd_sig(ix-1,1,1))*dx(ix-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           P = ppr*(1.0+1.5*abs(mu))
!--
           if (r1 < P) then
              ptcl%itype = 2
              grd_methodswap(ix,1,1)=grd_methodswap(ix,1,1)+1
              if(grd_isvelocity) then
!-- velocity effects accounting
                 tot_evelo = tot_evelo+e*(1d0-elabfact)
!
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              ix = ix-1
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = max(r1,r2)
              if(grd_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
              r = grd_xarr(ix)
           endif
        ! End of check
        else
           r = grd_xarr(ix)
           ix = ix-1
        endif
     endif!}}}
  elseif (d == dcen) then
     prt_done = .true.
     grd_numcensus(ix,1,1) = grd_numcensus(ix,1,1)+1
!     tot_erad = tot_erad + e*elabfact
!
  endif

end subroutine transport1
