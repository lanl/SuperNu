subroutine transport2(ptcl,isvacant)

  use gridmod
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
  !them through one IMC transport event.  If
  !the puretran boolean is set to false, this routine couples to the
  !corresponding DDMC diffusion routine.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch

  integer :: ig, imu
  real*8 :: dirdotu, mu0
  real*8 :: dtinv, elabfact, thelp, thelpinv 
  real*8 :: dcen,dcol,dthm,db,dbr,dbz,ddop,d
  real*8 :: rold, zold, omold, ppl, ppr, help
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
!
!-- calculating distance to census:
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%t)*thelpinv)
!
!-- calculating distance to boundary:
!-- to x-bound
  if(abs(mu)==1d0) then
!-- making greater than dcen
     dbr = 2d0*pc_c*tsp_dt*thelpinv
  else
     if(abs(sin(om))<grd_xarr(ix)/x .and. cos(om)<0d0) then
!-- inner boundary
        dbr = abs(x*cos(om)/sqrt(1d0-mu**2) &
             +sqrt(((cos(om)*x)**2-x**2+grd_xarr(ix)**2)/(1d0-mu**2)))
     elseif(abs(grd_xarr(ix+1)-x)<1d-15*x .and. cos(om)>0d0) then
!-- on outer boundary moving out
        dbr = 0d0
     else
!-- outer boundary
        dbr = -x*cos(om)/sqrt(1d0-mu**2) &
             + sqrt(((cos(om)*x)**2 + grd_xarr(ix+1)**2-x**2)/(1d0-mu**2))
     endif
  endif
  if(dbr/=dbr) stop 'transport2: dbr nan'

!-- to y-bound
  if(mu>0d0) then
     dbz = (grd_yarr(iy+1)-y)/mu
  elseif(mu<0d0) then
     dbz = (grd_yarr(iy)-y)/mu
  else
!-- making greater than dcen
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  endif

!-- finding minim boundary distance
  db = min(dbr,dbz)
!
!-- calculating distance to Thomson scattering:
  if(grd_sig(ix,iy,1)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*grd_sig(ix,iy,1))
  else
!-- making greater than dcen
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
  if(dthm/=dthm) stop 'transport2: dthm nan'
!
!-- calculating distance to effective collision:
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
!-- calculating distance to Doppler shift
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
  d = min(dcen,db,dthm,dcol,ddop)
  if(d<0d0) then
     write(*,*) dcen,dbz,dbr,dthm,dcol,ddop
     stop 'transport2: negative distance'
  endif
!
!-- using min distance to stream particle to event location
  rold = x
  zold = y
!-- updating position
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

!-- updating time
  ptcl%t = ptcl%t + thelp*cinv*d
!
!-- updating transformation factors
  if(grd_isvelocity) then
     dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
     elabfact = 1d0 - dirdotu*cinv
  else
     dirdotu = 0d0
     elabfact = 1d0
  endif

!
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
!-- distance to census
  if(d==dcen) then
!-- censusing particle
     prt_done = .true.
     grd_numcensus(ix,iy,1) = grd_numcensus(ix,iy,1)+1

!
!-- distance to Thomson scatter
  elseif(d==dthm) then
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
!-- transforming to cmf, then to lab:
!-- wavelength
        wl = wl*(1d0-dirdotu*cinv)/elabfact
!-- energy weight
        e = e*elabfact/(1d0-dirdotu*cinv)
        e0 = e0*elabfact/(1d0-dirdotu*cinv)
     endif

!
!-- distance to effective collision
  elseif(d==dcol) then
!-- sampling
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=grd_fcoef(ix,iy,1)) then
!-- effective absorption:
!-- ending particle
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        grd_edep(ix,iy,1) = grd_edep(ix,iy,1) + e*elabfact
     else
!-- effectively scattered:
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
                sqrt(1d0-mu**2)*cos(om)+x*cinv)
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
!-- transforming to cmf, then to lab:
!-- energy weight
           e = e*elabfact/(1d0-dirdotu*cinv)
           e0 = e0*elabfact/(1d0-dirdotu*cinv)
        endif
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
        endif
!-- checking if DDMC in new group
        if((grd_cap(ig,ix,iy,1)+grd_sig(ix,iy,1)) * &
             min(dx(ix),dy(iy))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%itype = 2
           grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
!-- transforming to cmf
           if(grd_isvelocity) then
              e = e*(1d0-dirdotu*cinv)
              e0 = e0*(1d0-dirdotu*cinv)
              wl = wl/(1d0-dirdotu*cinv)
           endif
        endif
     endif

!
!-- distance to y-boundary
  elseif(d==dbz) then
     if(mu>=0d0) then
!-- checking if particle escapes top
        if(iy == grd_ny) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame flux group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           ig = binsrch(wl,flx_wl,flx_ng+1,0)
!-- checking group bounds
           if(ig>flx_ng.or.ig<1) then
              if(ig>flx_ng) then
                 ig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 ig=1
                 wl=flx_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           flx_luminos(ig,imu,1) = flx_luminos(ig,imu,1)+e*dtinv
           flx_lumdev(ig,imu,1) = flx_lumdev(ig,imu,1)+(e*dtinv)**2
           flx_lumnum(ig,imu,1) = flx_lumnum(ig,imu,1)+1
!-- checking if above cell is DDMC
        elseif((grd_cap(ig,ix,iy+1,1)+grd_sig(ix,iy+1,1)) * &
             min(dx(ix),dy(iy+1))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming y-cosine to cmf
           if(grd_isvelocity) then
              mu = (mu-y*cinv)/elabfact
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
           endif
           help = (grd_cap(ig,ix,iy+1,1)+grd_sig(ix,iy+1,1))*dy(iy+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppl*(1d0+1.5d0*abs(mu))) then
              ptcl%itype = 2
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
              if(grd_isvelocity) then
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              iy = iy + 1
           else
!-- resampling y-cosine
              r1 = rand()
              r2 = rand()
              mu = -max(r1,r2)
!-- resampling azimuthal
              r1 = rand()
              om = pc_pi2*r1
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
              endif
!-- fix position at boundary
              y = grd_yarr(iy+1)
           endif
        else
!-- IMC in upper cell
           y = grd_yarr(iy+1)
           iy = iy+1
        endif
!-- mu<0
     else
!-- checking if particle escapes bottom
        if(iy == 1) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame flux group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           ig = binsrch(wl,flx_wl,flx_ng+1,0)
!-- checking group bounds
           if(ig>flx_ng.or.ig<1) then
              if(ig>flx_ng) then
                 ig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 ig=1
                 wl=flx_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           flx_luminos(ig,imu,1) = flx_luminos(ig,imu,1)+e*dtinv
           flx_lumdev(ig,imu,1) = flx_lumdev(ig,imu,1)+(e*dtinv)**2
           flx_lumnum(ig,imu,1) = flx_lumnum(ig,imu,1)+1
!-- checking if lower cell is DDMC
        elseif((grd_cap(ig,ix,iy-1,1)+grd_sig(ix,iy-1,1)) * &
             min(dx(ix),dy(iy-1))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming y-cosine to cmf
           if(grd_isvelocity) then
              mu = (mu-y*cinv)/elabfact
              if(mu>1d0) then
                 mu = 1d0
              elseif(mu<-1d0) then
                 mu = -1d0
              endif
           endif
           help = (grd_cap(ig,ix,iy-1,1)+grd_sig(ix,iy-1,1)) * &
                dy(iy-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppr*(1d0+1.5d0*abs(mu))) then
              ptcl%itype = 2
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
              if(grd_isvelocity) then
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              iy = iy - 1
           else
!-- resampling y-cosine
              r1 = rand()
              r2 = rand()
              mu = max(r1,r2)
!-- resampling azimuthal
              r1 = rand()
              om = pc_pi2*r1
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
              endif
!-- fix position at boundary
              y = grd_yarr(iy)
           endif
        else
!-- IMC in lower cell
           y = grd_yarr(iy)
           iy = iy-1
        endif
     endif

!
!-- distance to x-boundary
  elseif(d==dbr) then
     if(cos(om)>=0d0) then
!-- checking if particle escapes at outer radius
        if(ix == grd_nx) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame flux group and polar bin
           imu = binsrch(mu,flx_mu,flx_nmu+1,0)
           ig = binsrch(wl,flx_wl,flx_ng+1,0)
!-- checking group bounds
           if(ig>flx_ng.or.ig<1) then
              if(ig>flx_ng) then
                 ig=flx_ng
                 wl=flx_wl(flx_ng+1)
              else
                 ig=1
                 wl=flx_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           flx_luminos(ig,imu,1) = flx_luminos(ig,imu,1)+e*dtinv
           flx_lumdev(ig,imu,1) = flx_lumdev(ig,imu,1)+(e*dtinv)**2
           flx_lumnum(ig,imu,1) = flx_lumnum(ig,imu,1)+1
!-- checking if outer cell is DDMC
        elseif((grd_cap(ig,ix+1,iy,1)+grd_sig(ix+1,iy,1)) * &
             min(dx(ix+1),dy(iy))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming x-cosine to cmf
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
           help = (grd_cap(ig,ix+1,iy,1)+grd_sig(ix+1,iy,1))*dx(ix+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppl*(1d0+1.5d0*abs(mu0))) then
              ptcl%itype = 2
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
              if(grd_isvelocity) then
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              ix = ix + 1
           else
!-- resampling direction
              r1 = rand()
              r2 = rand()
              mu0 = -max(r1,r2)
              r1 = rand()
              mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
              if(om<0d0) om = om+pc_pi2
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
              endif
              x = grd_xarr(ix+1)
           endif
        else
!-- IMC in outer cell
           x = grd_xarr(ix+1)
           ix = ix + 1
        endif
!-- cos(om)<0
     else
        if(ix==1) then
           write(*,*) om, omold, x, rold, db
           stop 'transport2: cos(om)<0 and ix=1'
        endif
        if((grd_cap(ig,ix-1,iy,1)+grd_sig(ix-1,iy,1)) * &
             min(dx(ix-1),dy(iy))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming x-cosine to cmf
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
           help = (grd_cap(ig,ix-1,iy,1)+grd_sig(ix-1,iy,1))*dx(ix-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppr*(1d0+1.5d0*abs(mu0))) then
              ptcl%itype = 2
              grd_methodswap(ix,iy,1)=grd_methodswap(ix,iy,1)+1
              if(grd_isvelocity) then
                 e = e*elabfact
                 e0 = e0*elabfact
                 wl = wl/elabfact
              endif
              ix = ix - 1
           else
!-- resampling direction
              r1 = rand()
              r2 = rand()
              mu0 = max(r1,r2)
              r1 = rand()
              mu = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(1d0-mu0**2)*sin(pc_pi2*r1),mu0)
              if(om<0d0) om = om+pc_pi2
              if(grd_isvelocity) then
                 dirdotu = mu*y+sqrt(1d0-mu**2)*cos(om)*x
                 om = atan2(sqrt(1d0-mu**2)*sin(om), &
                      sqrt(1d0-mu**2)*cos(om)+x*cinv)
!-- transforming y-axis direction cosine to lab
                 mu = (mu+y*cinv)/(1d0+dirdotu*cinv)
                 if(mu>1d0) then
                    mu = 1d0
                 elseif(mu<-1d0) then
                    mu = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(om<0d0) om=om+pc_pi2
              endif
              x = grd_xarr(ix)
           endif
        else
!-- IMC in inner cell
           x = grd_xarr(ix)
           ix = ix - 1
        endif
     endif

!
!-- distance to doppler shift
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
           e = e*elabfact
           e0 = e0*elabfact
           wl = wl/elabfact
        endif
     endif
  endif


end subroutine transport2
