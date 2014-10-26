subroutine transport2(ptcl,isvacant)

  use gasgridmod
  use timestepmod
  use physconstmod
  use particlemod
  use inputparmod
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

  integer :: g, ig, iig
  real*8 :: dirdotu, azidotu, mu0
  real*8 :: dtinv, elabfact, thelp, thelpinv 
  real*8 :: dcen,dcol,dthm,db,dbr,dbz,ddop,d
  real*8 :: rold, zold, omold, ppl, ppr, help
  real*8 :: r1, r2, denom2

  integer,pointer :: zr,zz
  real*8,pointer :: r,z,xi,om,ep,ep0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dy
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dy(l) = gas_yarr(l+1) - gas_yarr(l)

  zr => ptcl%zsrc
  zz => ptcl%iy
  r => ptcl%rsrc
  z => ptcl%y
  xi => ptcl%musrc
  om => ptcl%om
  ep => ptcl%esrc
  ep0 => ptcl%ebirth
  wl => ptcl%wlsrc
!
!-- shortcut
  dtinv = 1d0/tsp_dt
!
!-- setting vel-grid helper variables
  if(in_isvelocity) then
!-- calculating initial transformation factors
     dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
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
  g = binsrch(wl/elabfact,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
  if(g>gas_ng.or.g<1) then
     if(g==gas_ng+1) then
        g = gas_ng
     elseif(g==0) then
        g = 1
     else
        stop 'transport2 (1): particle group invalid'
     endif
  endif
!
!-- calculating distance to census:
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%tsrc)*thelpinv)
!
!-- calculating distance to boundary:
!-- to r-bound
  if(abs(xi)==1d0) then
!-- making greater than dcen
     dbr = 2d0*pc_c*tsp_dt*thelpinv
  else
     if(abs(sin(om))<gas_xarr(zr)/r.and.cos(om)<0d0) then
!-- inner boundary
        dbr = abs(r*cos(om)/sqrt(1d0-xi**2) &
             +sqrt(((cos(om)*r)**2-r**2+gas_xarr(zr)**2)/(1d0-xi**2)))
        if(dbr/=dbr) then
           write(*,*) ((cos(om)*r)**2-r**2+gas_xarr(zr)**2)/(1d0-xi**2)
           stop 'transport2: invalid inner dbr'
        endif
     else
!-- outer boundary
        dbr = -r*cos(om)/sqrt(1d0-xi**2) &
             +sqrt(((cos(om)*r)**2+gas_xarr(zr+1)**2-r**2)/(1d0-xi**2))
        if(dbr/=dbr) then
           write(*,*) ((cos(om)*r)**2-r**2+gas_xarr(zr+1)**2)/(1d0-xi**2)
           stop 'transport2: invalid outer dbr'
        endif
     endif
  endif

!-- to z-bound
  if(xi>0d0) then
     dbz = (gas_yarr(zz+1)-z)/xi
     if(dbz<0d0) then
        write(*,*) zz, gas_yarr(zz+1), z, xi
        stop 'upward dbz'
     endif
  elseif(xi<0d0) then
     dbz = (gas_yarr(zz)-z)/xi
     if(dbz<0d0) stop 'downward dbz'
  else
!-- making greater than dcen
     dbz = 2d0*pc_c*tsp_dt*thelpinv
  endif

!-- finding minim boundary distance
  db = min(dbr,dbz)
!
!-- calculating distance to Thomson scattering:
  if(gas_sig(zr,zz,1)>0d0) then
     r1 = rand()
     dthm = -log(r1)*thelpinv/(elabfact*gas_sig(zr,zz,1))
  else
!-- making greater than dcen
     dthm = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- calculating distance to effective collision:
  if(gas_cap(g,zr,zz,1)<=0d0) then
!-- making greater than dcen
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  elseif(prt_isimcanlog) then
!-- calculating dcol for analog MC
     r1 = rand()
     dcol = -log(r1)*thelpinv/(elabfact*gas_cap(g,zr,zz,1))
  elseif(gas_fcoef(zr,zz,1)<1d0.and.gas_fcoef(zr,zz,1)>=0d0) then
     r1 = rand()
     dcol = -log(r1)*thelpinv/&
          (elabfact*(1d0-gas_fcoef(zr,zz,1))*gas_cap(g,zr,zz,1))
  else
!-- making greater than dcen
     dcol = 2d0*pc_c*tsp_dt*thelpinv
  endif
!
!-- calculating distance to Doppler shift
  if(in_isvelocity.and.g<gas_ng) then
     ddop = pc_c*(elabfact-wl/gas_wl(g+1))
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
  if(any((/dcen,dbz,dbr,dthm,dcol,ddop/)<0d0)) then
     write(*,*) dcen,dbz,dbr,dthm,dcol,ddop
     stop 'transport2: negative distance'
  endif
!
!-- using min distance to stream particle to event location
  rold = r
  zold = z
!-- updating position
  r = sqrt(rold**2+(1d0-xi**2)*d**2+2d0*rold*sqrt(1d0-xi**2)*d*cos(om))
  z = zold+xi*d
!-- updating azimuthal direction
  omold = om
  if(om>pc_pi.and.om<pc_pi2) then
     om = pc_pi2 + &
          atan2(-sqrt(max(r**2-(rold*cos(omold)+d*sqrt(1d0-xi**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-xi**2))
  else
     om = atan2(sqrt(max(r**2-(rold*cos(omold)+d*sqrt(1d0-xi**2))**2,0d0)), &
          rold*cos(omold)+d*sqrt(1d0-xi**2))
  endif
  if(om/=om) then
     write(*,*) d, r, rold, omold, om, xi
     stop 'transport2: om is nan'
  endif

!-- updating time
  ptcl%tsrc = ptcl%tsrc + thelp*cinv*d
!
!-- updating transformation factors
  if(in_isvelocity) then
     dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
     elabfact = 1d0 - dirdotu*cinv
  else
     dirdotu = 0d0
     elabfact = 1d0
  endif

!
!-- tallying energy densities
  if(prt_isimcanlog) then
!-- analog energy density
     gas_eraddens(zr,zz,1)=gas_eraddens(zr,zz,1)+ep*elabfact* &
          d*thelp*cinv*dtinv
  else
!-- nonanalog energy density
     if(gas_fcoef(zr,zz,1)*gas_cap(g,zr,zz,1)* &
          min(dx(zr),dy(zz))*thelp>1d-6) then
        gas_eraddens(zr,zz,1) = gas_eraddens(zr,zz,1)+ep* &
             (1.0d0-exp(-gas_fcoef(zr,zz,1)*elabfact* &
             gas_cap(g,zr,zz,1)*d*thelp))* &
             elabfact/(gas_fcoef(zr,zz,1)*elabfact * &
             gas_cap(g,zr,zz,1)*pc_c*tsp_dt)
     else
!-- analog energy density
        gas_eraddens(zr,zz,1)=gas_eraddens(zr,zz,1)+ep*elabfact* &
             d*thelp*cinv*dtinv
     endif
!-- depositing nonanalog absorbed energy
     gas_edep(zr,zz,1)=gas_edep(zr,zz,1)+ep* &
          (1d0-exp(-gas_fcoef(zr,zz,1)*gas_cap(g,zr,zz,1)* &
          elabfact*d*thelp))*elabfact
     if(gas_edep(zr,zz,1)/=gas_edep(zr,zz,1)) then
        stop 'transport2: invalid energy deposition'
     endif
!-- reducing particle energy
     ep = ep*exp(-gas_fcoef(zr,zz,1)*gas_cap(g,zr,zz,1) * &
          elabfact*d*thelp)
  endif

!
!-- checking which event occurs from min distance

!
!-- distance to census
  if(d==dcen) then
!-- censusing particle
     prt_done = .true.
     gas_numcensus(zr,zz,1) = gas_numcensus(zr,zz,1)+1

!
!-- distance to Thomson scatter
  elseif(d==dthm) then
!-- resampling direction
     r1 = rand()
     xi = 1d0 - 2d0*r1
     r1 = rand()
     om = pc_pi2*r1
!-- checking velocity dependence
     if(in_isvelocity) then
!-- calculating transformation factors
        dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
        azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
             sqrt(1d0-xi**2)*cos(om)+r/pc_c)
!-- transforming to lab:
!-- z-cosine
        xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
        if(xi>1d0) then
           xi = 1d0
        elseif(xi<-1d0) then
           xi = -1d0
        endif
!-- azimuthal direction angle
        if(azidotu<0d0) then
           om = azidotu+pc_pi2
        else
           om = azidotu
        endif
!-- recalculating dirdotu
        dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
!-- transforming to cmf, then to lab:
!-- wavelength
        wl = wl*(1d0-dirdotu*cinv)/elabfact
!-- energy weight
        ep = ep*elabfact/(1d0-dirdotu*cinv)
        ep0 = ep0*elabfact/(1d0-dirdotu*cinv)
     endif

!
!-- distance to effective collision
  elseif(d==dcol) then
!-- sampling
     r1 = rand()
!-- checking if analog
     if(prt_isimcanlog.and.r1<=gas_fcoef(zr,zz,1)) then
!-- effective absorption:
!-- ending particle
        isvacant=.true.
        prt_done=.true.
!-- adding comoving energy to deposition energy
        gas_edep(zr,zz,1) = gas_edep(zr,zz,1) + ep*elabfact
     else
!-- effectively scattered:
!-- resampling direction
        r1 = rand()
        xi = 1d0 - 2d0*r1
        r1 = rand()
        om = pc_pi2*r1
!-- checking velocity dependence
        if(in_isvelocity) then
!-- calculating transformation factors
           dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
           azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming to lab:
!-- z-cosine
           xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
           if(xi>1d0) then
              xi = 1d0
           elseif(xi<-1d0) then
              xi = -1d0
           endif
!-- azimuthal direction angle
           if(azidotu<0d0) then
              om = azidotu+pc_pi2
           else
              om = azidotu
           endif
!-- recalculating dirdotu
           dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
!-- transforming to cmf, then to lab:
!-- energy weight
           ep = ep*elabfact/(1d0-dirdotu*cinv)
           ep0 = ep0*elabfact/(1d0-dirdotu*cinv)
        endif
!-- redistributing wavelength
        denom2 = 0d0
        r1 = rand()
        do ig = 1, gas_ng
           iig = ig
           if ((r1>=denom2).and.(r1<denom2+gas_emitprob(ig,zr,zz,1))) exit
           denom2 = denom2+gas_emitprob(ig,zr,zz,1)
        enddo
        g = iig
!-- uniformly in new group
        r1 = rand()
        wl = 1d0/((1d0-r1)/gas_wl(g)+r1/gas_wl(g+1))
!-- transforming to lab
        if(in_isvelocity) then
           wl = wl*(1d0-dirdotu*cinv)
        endif
!-- checking if DDMC in new group
        if((gas_cap(g,zr,zz,1)+gas_sig(zr,zz,1)) * &
             min(dx(zr),dy(zz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
           ptcl%rtsrc = 2
           gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!-- transforming to cmf
           if(in_isvelocity) then
              ep = ep*(1d0-dirdotu*cinv)
              ep0 = ep0*(1d0-dirdotu*cinv)
              wl = wl/(1d0-dirdotu*cinv)
           endif
        endif
     endif

!
!-- distance to z-boundary
  elseif(d==dbz) then
     if(xi>=0d0) then
!-- checking if particle escapes top
        if(zz == gas_ny) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame group
           g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
           if(g>gas_ng.or.g<1) then
              if(g>gas_ng) then
                 g=gas_ng
                 wl=gas_wl(gas_ng+1)
              else
                 g=1
                 wl=gas_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           gas_luminos(g) = gas_luminos(g)+ep*dtinv
           gas_lumdev(g) = gas_lumdev(g)+(ep0*dtinv)**2
           gas_lumnum(g) = gas_lumnum(g)+1
!-- checking if above cell is DDMC
        elseif((gas_cap(g,zr,zz+1,1)+gas_sig(zr,zz+1,1)) * &
             min(dx(zr),dy(zz+1))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming z-cosine to lab
           if(in_isvelocity) then
              xi = (xi-z*cinv)/elabfact
              if(xi>1d0) then
                 xi = 1d0
              elseif(xi<-1d0) then
                 xi = -1d0
              endif
           endif
           help = (gas_cap(g,zr,zz+1,1)+gas_sig(zr,zz+1,1))*dy(zz+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppl*(1d0+1.5d0*abs(xi))) then
              ptcl%rtsrc = 2
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
              if(in_isvelocity) then
                 ep = ep*elabfact
                 ep0 = ep0*elabfact
                 wl = wl/elabfact
              endif
              zz = zz + 1
           else
!-- resampling z-cosine
              r1 = rand()
              r2 = rand()
              xi = -max(r1,r2)
!-- resampling azimuthal
              r1 = rand()
              om = pc_pi2*r1
              if(in_isvelocity) then
                 dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
                 azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                      sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming z-axis direction cosine to lab
                 xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
                 if(xi>1d0) then
                    xi = 1d0
                 elseif(xi<-1d0) then
                    xi = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(azidotu<0d0) then
                    om = azidotu+pc_pi2
                 else
                    om = azidotu
                 endif
              endif
           endif
        else
!-- IMC in upper cell
           zz = zz+1
        endif
!-- xi<0
     else
!-- checking if particle escapes bottom
        if(zz == 1) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame group
           g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
           if(g>gas_ng.or.g<1) then
              if(g>gas_ng) then
                 g=gas_ng
                 wl=gas_wl(gas_ng+1)
              else
                 g=1
                 wl=gas_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           gas_luminos(g) = gas_luminos(g)+ep*dtinv
           gas_lumdev(g) = gas_lumdev(g)+(ep0*dtinv)**2
           gas_lumnum(g) = gas_lumnum(g)+1
!-- checking if lower cell is DDMC
        elseif((gas_cap(g,zr,zz-1,1)+gas_sig(zr,zz-1,1)) * &
             min(dx(zr),dy(zz-1))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming z-cosine to lab
           if(in_isvelocity) then
              xi = (xi-z*cinv)/elabfact
              if(xi>1d0) then
                 xi = 1d0
              elseif(xi<-1d0) then
                 xi = -1d0
              endif
           endif
           help = (gas_cap(g,zr,zz-1,1)+gas_sig(zr,zz-1,1)) * &
                dy(zz-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppr*(1d0+1.5d0*abs(xi))) then
              ptcl%rtsrc = 2
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
              if(in_isvelocity) then
                 ep = ep*elabfact
                 ep0 = ep0*elabfact
                 wl = wl/elabfact
              endif
              zz = zz - 1
           else
!-- resampling z-cosine
              r1 = rand()
              r2 = rand()
              xi = max(r1,r2)
!-- resampling azimuthal
              r1 = rand()
              om = pc_pi2*r1
              if(in_isvelocity) then
                 dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
                 azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                      sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming z-axis direction cosine to lab
                 xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
                 if(xi>1d0) then
                    xi = 1d0
                 elseif(xi<-1d0) then
                    xi = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(azidotu<0d0) then
                    om = azidotu+pc_pi2
                 else
                    om = azidotu
                 endif
              endif
           endif
        else
!-- IMC in lower cell
           zz = zz-1
        endif
     endif

!
!-- distance to r-boundary
  elseif(d==dbr) then
     if(cos(om)>=0d0) then
!-- checking if particle escapes at outer radius
        if(zr == gas_nx) then
!-- ending particle
           isvacant = .true.
           prt_done = .true.
!-- retrieving lab frame group
           g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
           if(g>gas_ng.or.g<1) then
              if(g>gas_ng) then
                 g=gas_ng
                 wl=gas_wl(gas_ng+1)
              else
                 g=1
                 wl=gas_wl(1)
              endif
           endif
!-- tallying outbound luminosity
           gas_luminos(g) = gas_luminos(g)+ep*dtinv
           gas_lumdev(g) = gas_lumdev(g)+(ep0*dtinv)**2
           gas_lumnum(g) = gas_lumnum(g)+1
!-- checking if outer cell is DDMC
        elseif((gas_cap(g,zr+1,zz,1)+gas_sig(zr+1,zz,1)) * &
             min(dx(zr+1),dy(zz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming r-cosine to cmf
           if(in_isvelocity) then
              azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                   sqrt(1d0-xi**2)*cos(om)-r*cinv)
              if(azidotu<0d0) then
                 om = azidotu+pc_pi2
              else
                 om = azidotu
              endif
              xi = (xi-z*cinv)/elabfact
              if(xi>1d0) then
                 xi = 1d0
              elseif(xi<-1d0) then
                 xi = -1d0
              endif
           endif
!-- r-cosine
           mu0 = sqrt(1d0-xi**2)*cos(om)
           help = (gas_cap(g,zr+1,zz,1)+gas_sig(zr+1,zz,1))*dx(zr+1)*thelp
           ppl = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppl*(1d0+1.5d0*abs(mu0))) then
              ptcl%rtsrc = 2
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
              if(in_isvelocity) then
                 ep = ep*elabfact
                 ep0 = ep0*elabfact
                 wl = wl/elabfact
              endif
              zr = zr + 1
           else
!-- resampling direction
              r1 = rand()
              r2 = rand()
              mu0 = -max(r1,r2)
              r1 = rand()
              xi = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(max(1d0-xi**2-mu0**2,0d0)),mu0)
              if(om<0d0) then
                 om = om+pc_pi2
              endif
!-- extending azimuthal sample
              r1 = rand()
              if(r1 < 0.5d0) then
                 om = pc_pi2-om
              endif
              if(in_isvelocity) then
                 dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
                 azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                      sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming z-axis direction cosine to lab
                 xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
                 if(xi>1d0) then
                    xi = 1d0
                 elseif(xi<-1d0) then
                    xi = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(azidotu<0d0) then
                    om = azidotu+pc_pi2
                 else
                    om = azidotu
                 endif
              endif
           endif
        else
!-- IMC in outer cell
           zr = zr + 1
           if(abs(r-gas_xarr(zr))/gas_xarr(zr)<1d-10) then
              r = gas_xarr(zr)
           else
              write(*,*) db, rold, r, gas_xarr(zr)
              stop 'transport2: outer db'
           endif
        endif
!-- cos(om)<0
     else
        if(zr==1) then
           write(*,*) om, omold, r, rold, db
           stop 'transport2: cos(om)<0 and zr=1'
        endif
        if((gas_cap(g,zr-1,zz,1)+gas_sig(zr-1,zz,1)) * &
             min(dx(zr-1),dy(zz))*thelp >= prt_tauddmc &
             .and..not.in_puretran) then
!-- transforming r-cosine to cmf
           if(in_isvelocity) then
              azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                   sqrt(1d0-xi**2)*cos(om)-r*cinv)
              if(azidotu<0d0) then
                 om = azidotu+pc_pi2
              else
                 om = azidotu
              endif
              xi = (xi-z*cinv)/elabfact
              if(xi>1d0) then
                 xi = 1d0
              elseif(xi<-1d0) then
                 xi = -1d0
              endif
           endif
!-- r-cosine
           mu0 = sqrt(1d0-xi**2)*cos(om)
           help = (gas_cap(g,zr-1,zz,1)+gas_sig(zr-1,zz,1))*dx(zr-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
!-- sampling
           r1 = rand()
           if (r1 < ppr*(1d0+1.5d0*abs(mu0))) then
              ptcl%rtsrc = 2
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
              if(in_isvelocity) then
                 ep = ep*elabfact
                 ep0 = ep0*elabfact
                 wl = wl/elabfact
              endif
              zr = zr - 1
           else
!-- resampling direction
              r1 = rand()
              r2 = rand()
              mu0 = max(r1,r2)
              r1 = rand()
              xi = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(max(1d0-xi**2-mu0**2,0d0)),mu0)
              if(om<0d0) then
                 om = om+pc_pi2
              endif
!-- extending azimuthal sample
              r1 = rand()
              if(r1 < 0.5d0) then
                 om = pc_pi2-om
              endif
              if(in_isvelocity) then
                 dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
                 azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                      sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming z-axis direction cosine to lab
                 xi = (xi+z*cinv)/(1d0+dirdotu*cinv)
                 if(xi>1d0) then
                    xi = 1d0
                 elseif(xi<-1d0) then
                    xi = -1d0
                 endif
!-- transforming azimuthal angle to lab
                 if(azidotu<0d0) then
                    om = azidotu+pc_pi2
                 else
                    om = azidotu
                 endif
              endif
           endif
        else
!-- IMC in inner cell
           zr = zr - 1
           if(abs(r-gas_xarr(zr+1))/gas_xarr(zr+1)<1d-10) then
              r = gas_xarr(zr+1)
           else
              write(*,*) db, rold, r, gas_xarr(zr+1), sin(om)
              write(*,*) zr, gas_xarr(zr+1)/rold
              stop 'transport2: inner db'
           endif
        endif
     endif

!
!-- distance to doppler shift
  elseif(d==ddop) then
     if(.not.in_isvelocity) stop 'transport2: ddop and no velocity'
     if(g<gas_ng) then
!-- shifting group
        g = g+1
        wl = (gas_wl(g)+1d-6*(gas_wl(g+1)-gas_wl(g)))*elabfact
     else
!-- resampling wavelength in highest group
        r1 = rand()
        wl=1d0/(r1/gas_wl(g+1) + (1d0-r1)/gas_wl(gas_ng))
        wl = wl*elabfact
     endif
!-- check if ddmc region
     if ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1)) * &
          min(dx(zr),dy(zz))*thelp >= prt_tauddmc &
          .and..not.in_puretran) then
        ptcl%rtsrc = 2
        gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
        if(in_isvelocity) then
           ep = ep*elabfact
           ep0 = ep0*elabfact
           wl = wl/elabfact
        endif
     endif
  endif


end subroutine transport2
