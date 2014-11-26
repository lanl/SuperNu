subroutine transport2(ptcl,isvacant)

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
  logical,intent(inout) :: isvacant
!##################################################
  !This subroutine passes particle parameters as input and modifies
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch

  logical :: loutx,louty
  integer :: ig, imu, ihelp
  real*8 :: elabfact, dirdotu, mu0
  real*8 :: dtinv, thelp, thelpinv, help
  real*8 :: dcen,dcol,dthm,dbx,dby,ddop,d
  real*8 :: rold, zold, omold
  real*8 :: r1, r2, denom2

  integer,pointer :: ix,iy
  real*8,pointer :: x,y,mu,om,e,e0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)

  ix => ptcl%ix
  iy => ptcl%iy
  x => ptcl%x
  y => ptcl%y
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
  wl => ptcl%wl
!
!-- shortcut
  dtinv = 1d0/tsp_dt
!
!-- setting vel-grid helper variables
  if(grd_isvelocity) then
!-- calculating initial transformation factors
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
     thelp = tsp_t
  else
     dirdotu = 0d0
     elabfact = 1d0
     thelp = 1d0
  endif
!
!-- inverting vel-grid factor
  thelpinv = 1d0/thelp
!
!-- looking up initial group
  ig = binsrch(wl/elabfact,grd_wl,grd_ng+1,in_ng)
!-- checking group bounds
  if(ig>grd_ng.or.ig<1) then
     if(ig>grd_ng) then
        ig = grd_ng
     elseif(ig<1) then
        ig = 1
     else
        write(0,*) ig,grd_ng,wl,elabfact
        stop 'transport2 (1): particle group invalid'
     endif
  endif

!-- census distance
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- boundary distances
!-- to x-bound
  if(abs(mu)==1d0) then
!-- making greater than dcen
     dbx = 2d0*pc_c*tsp_dt*thelpinv
  else
     if(abs(sin(om))<grd_xarr(ix)/x .and. cos(om)<0d0) then
!-- inner boundary
        dbx = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbx = 0d0
     else
!-- outer boundary
        dbx = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
     endif
  endif
  if(dbx/=dbx) stop 'transport2: dbx nan'

!-- to y-bound
  if(mu>0d0) then
     dby = (grd_yarr(iy+1)-y)/mu
  elseif(mu<0d0) then
     dby = (grd_yarr(iy)-y)/mu
  else
!-- making greater than dcen
     dby = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- Thomson scattering distance
  if(grd_sig(ix,iy,1)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ix,iy,1))
  else
!-- making greater than dcen
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
  if(dthm/=dthm) stop 'transport2: dthm nan'
!
!-- effective collision distance
  if(grd_cap(ig,ix,iy,1)<=0d0) then
!-- making greater than dcen
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*grd_cap(ig,ix,iy,1))
  elseif(grd_fcoef(ix,iy,1)<1d0.and.grd_fcoef(ix,iy,1)>=0d0) then
     r1 = rand()
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-grd_fcoef(ix,iy,1))*grd_cap(ig,ix,iy,1))
  else
!-- making greater than dcen
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
  if(dcol/=dcol) stop 'transport2: dthm nan'
!
!-- Doppler shift distance
  if(grd_isvelocity.and.ig<grd_ng) then
     ddop = pc_c*(elabfact-wl/grd_wl(ig+1))
     if(ddop<0d0) then
        ddop = 2d0*pc_c*tsp_dt*thelpinv
     endif
  else
!-- making greater than dcen
     ddop = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- finding minimum distance
  d = min(dcen,dbx,dby,dthm,dcol,ddop)
  if(d<0d0) then
     write(*,*) dcen,dby,dbx,dthm,dcol,ddop
     stop 'transport2: negative distance'
  endif
!
!-- updating position
  rold = x
  zold = y
  x = sqrt(rold**2+(1d0-mu**2)*d**2+2d0*rold*sqrt(1d0-mu**2)*d*cos(om))
  y = zold+mu*d
!-- updating azimuthal direction
  omold = om
  if(om>pc_pi.and.om<pc_pi2) then
     om = pc_pi2 + &
          atan2(-sqrt(max(x**2-(rold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-mu**2))
  else
     om = atan2(sqrt(max(x**2-(rold*cos(omold)+d*sqrt(1d0-mu**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-mu**2))
  endif
  if(om/=om) then
     write(*,*) d, x, rold, omold, om, mu
     stop 'transport2: om is nan'
  endif
!
!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d
!
!-- updating transformation factors
  if(grd_isvelocity) then
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
  endif

!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     grd_eraddens(ix,iy,1)=grd_eraddens(ix,iy,1)+e*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(grd_fcoef(ix,iy,1)*grd_cap(ig,ix,iy,1)* &
          min(dx(ix),dy(iy))*thelp>1d-6) then
        grd_eraddens(ix,iy,1) = grd_eraddens(ix,iy,1)+e* &
             (1.0d0-exp(-grd_fcoef(ix,iy,1)*elabfact* &
             grd_cap(ig,ix,iy,1)*d*thelp))* &
             elabfact/(grd_fcoef(ix,iy,1)*elabfact * &
             grd_cap(ig,ix,iy,1)*pc_c*tsp_dt)
     else
!-- analog energy density
        grd_eraddens(ix,iy,1)=grd_eraddens(ix,iy,1)+e*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     grd_edep(ix,iy,1)=grd_edep(ix,iy,1)+e* &
          (1d0-exp(-grd_fcoef(ix,iy,1)*grd_cap(ig,ix,iy,1)* &
          elabfact*d*thelp))*elabfact
     if(grd_edep(ix,iy,1)/=grd_edep(ix,iy,1)) then
        stop 'transport2: invalid energy deposition'
     endif
!-- reducing particle energy
     e = e*exp(-grd_fcoef(ix,iy,1)*grd_cap(ig,ix,iy,1) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!
!-- census
  if(d==dcen) then
     prt_done = .true.
     grd_numcensus(ix,iy,1) = grd_numcensus(ix,iy,1)+1
     return
  endif

!-- common manipulations for collisions
  if(d==dthm.or.d==dcol) then
!-- resampling direction
     r1 = rand()
     mu = 1d0 - 2d0*r1
     r1 = rand()
     om = pc_pi2*r1
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- calculating transformation factors
        dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
        om = atan2(sqrt(1d0-mu**2)*sin(om), &
             sqrt(1d0-mu**2)*cos(om)+x/pc_c)
!-- transforming to lab:
!-- y-cosine
        mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
        if(mu>1d0) then
           mu = 1d0
        elseif(mu<-1d0) then
           mu = -1d0
        endif
!-- azimuthal direction angle
        if(om<0d0) om=om+pc_pi2
!-- recalculating dirdotu
        dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     endif
  elseif(any([dbx,dby]==d)) then
!-- checking if escapted domain
     loutx = d==dbx.and.(cos(om)>=0d0.and.ix==grd_nx)
     louty = d==dby.and.((mu>=0d0.and.iy==grd_ny).or.(mu<0.and.iy==1))
     if(loutx.or.louty) then
!-- ending particle
        isvacant = .true.
        prt_done = .true.
!-- retrieving lab frame flux group, polar bin
        imu = binsrch(mu,flx_mu,flx_nmu+1,0)
        ig = binsrch(wl,flx_wl,flx_ng+1,0)
!-- checking group bounds
        if(ig>flx_ng.or.ig<1) then
           if(ig>flx_ng) then
              ig=flx_ng
           else
              ig=1
           endif
        endif
!-- tallying outbound luminosity
        flx_luminos(ig,imu,1) = flx_luminos(ig,imu,1)+e*dtinv
        flx_lumdev(ig,imu,1) = flx_lumdev(ig,imu,1)+(e*dtinv)**2
        flx_lumnum(ig,imu,1) = flx_lumnum(ig,imu,1)+1
        return
     endif
  endif

!
!-- Thomson scatter
  if(d==dthm) then
!-- checking velocity dependence
     if(grd_isvelocity) then
!-- wavelength
        wl = wl*(1d0-dirdotu*cinv)/elabfact
        help = elabfact/(1d0-dirdotu*cinv)
!-- velocity effects accounting
        tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
        e = e*help
        e0 = e0*help
     endif

!
!-- effective collision
  elseif(d==dcol) then
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ix,iy,1)) then
!-- effective absorption
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,1) = grd_edep(ix,iy,1) + e*elabfact
        return
     else
!-- effective scattering
!-- redistributing wavelength
        denom2 = 0d0
        r1 = rand()
        do ig = 1, grd_ng
           if ((r1>=denom2).and.(r1<denom2+grd_emitprob(ig,ix,iy,1))) exit
           denom2 = denom2+grd_emitprob(ig,ix,iy,1)
        enddo
!-- uniformly in new group
        r1 = rand()
        wl = 1d0/((1d0-r1)/grd_wl(ig)+r1/grd_wl(ig+1))
!-- transforming to lab
        if(grd_isvelocity) then
           wl = wl*(1d0-dirdotu*cinv)
           help = elabfact/(1d0-dirdotu*cinv)
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-help)
!
!-- energy weight
           e = e*help
           e0 = e0*help
        endif
!-- checking if DDMC in new group
        if((grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1)) * &
             min(dx(ix),dy(iy))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- transforming to cmf
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*dirdotu*cinv
!
              e = e*(1d0-dirdotu*cinv)
              e0 = e0*(1d0-dirdotu*cinv)
              wl = wl/(1d0-dirdotu*cinv)
           endif
        endif
     endif

!
!-- x-bound
  elseif(d==dbx) then

     if(cos(om)>=0d0) then
        ihelp = 1
        x = grd_xarr(ix+1)
     else
        if(ix==1) stop &
             'transport2_gamgrey: cos(om)<0 and ix=1'
        ihelp = -1
        x = grd_xarr(ix)
     endif

     if((grd_cap(ig,ix+ihelp,iy,1)+grd_sig(ix+ihelp,iy,1)) * &
          min(dx(ix+ihelp),dy(iy))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        ix = ix+ihelp
     else
!-- DDMC in adjacent cell
        if(grd_isvelocity) then
           om = atan2(sqrt(1d0-mu**2)*sin(om), &
                sqrt(1d0-mu**2)*cos(om)-x*cinv)
           if(om<0d0) om=om+pc_pi2
           mu = (mu-y*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
        endif
!-- x-cosine
        mu0 = sqrt(1d0-mu**2)*cos(om)
        help = (grd_cap(ig,ix+ihelp,iy,1)+grd_sig(ix+ihelp,iy,1)) * &
             dx(ix+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rand()
        if (r1 < help*(1d0+1.5d0*abs(mu0))) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           ix = ix + ihelp
        else
           r1 = rand()
           r2 = rand()
           mu0 = -ihelp*max(r1,r2)
           r1 = rand()
!-- resampling y-cosine
           mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
!-- resampling azimuthal
           om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
           if(om<0d0) om=om+pc_pi2
           if(grd_isvelocity) then
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              om = atan2(sqrt(1d0-mu**2)*sin(om), &
                   sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming mu to lab
              mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- y-bound
  elseif(d==dby) then

     if(mu>=0d0) then
        ihelp = 1
        y = grd_yarr(iy+1)
     else
        ihelp = -1
        y = grd_yarr(iy)
     endif

     if((grd_cap(ig,ix,iy+ihelp,1)+grd_sig(ix,iy+ihelp,1)) * &
          min(dx(ix),dy(iy+ihelp))*thelp < prt_tauddmc &
          .or.in_puretran) then
!-- IMC in adjacent cell
        iy = iy+ihelp
     else
        if(grd_isvelocity) then
!-- transforming y-cosine to cmf
           mu = (mu-y*cinv)/elabfact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
        endif
        help = (grd_cap(ig,ix,iy+ihelp,1)+grd_sig(ix,iy+ihelp,1)) * &
             dy(iy+ihelp)*thelp
        help = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
        r1 = rand()
        if (r1 < help*(1d0+1.5d0*abs(mu))) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
           if(grd_isvelocity) then
!-- velocity effects accounting
              tot_evelo=tot_evelo+e*(1d0-elabfact)
!
              e = e*elabfact
              e0 = e0*elabfact
              wl = wl/elabfact
           endif
           iy = iy + ihelp
        else
           r1 = rand()
           r2 = rand()
           mu = -ihelp*max(r1,r2)
!-- resampling azimuthal
           r1 = rand()
           om = pc_pi2*r1
           if(grd_isvelocity) then
              dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
              om = atan2(sqrt(1d0-mu**2)*sin(om), &
                   sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming mu to lab
              mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
!-- transforming om to lab
              if(om<0d0) om=om+pc_pi2
           endif
        endif
     endif

!
!-- Doppler shift
  elseif(d==ddop) then
     if(.not.grd_isvelocity) stop 'transport2: ddop and no velocity'
     if(ig<grd_ng) then
!-- shifting group
        ig = ig+1
        wl = (grd_wl(ig)+1d-6*(grd_wl(ig+1)-grd_wl(ig)))*elabfact
     else
!-- resampling wavelength in highest group
        r1 = rand()
        wl=1d0/(r1/grd_wl(grd_ng+1) + (1d0-r1)/grd_wl(grd_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if ((grd_sig(ix,iy,1)+grd_cap(ig,ix,iy,1)) * &
          min(dx(ix),dy(iy))*thelp >= prt_tauddmc &
          .and..not.in_puretran) then
        ptcl%itype = 2
        grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
        if(grd_isvelocity) then
!-- velocity effects accounting
           tot_evelo=tot_evelo+e*(1d0-elabfact)
!
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  else
     stop 'transport2: invalid distance'
  endif

end subroutine transport2
