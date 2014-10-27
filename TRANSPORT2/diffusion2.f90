subroutine diffusion2(ptcl,isvacant)

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
  !them through one DDMC diffusion event (Densmore, 2007).  If
  !the puretran boolean is set to false, this routine couples to the
  !analogous IMC transport routine through the advance. If puretran
  !is set to true, this routine is not used.
!##################################################
  real*8,parameter :: cinv = 1d0/pc_c
  integer, external :: binsrch
!
  integer :: ig, iig, g
  logical :: lhelp
  real*8 :: r1, r2, thelp, mu0
  real*8 :: denom, denom2, denom3
  real*8 :: ddmct, tau, tcensus
  real*8 :: elabfact, dirdotu, azidotu
  real*8 :: PU, PD, PR, PL, PA
!-- lumped quantities -----------------------------------------

  real*8 :: emitlump, speclump
  real*8 :: caplump
  real*8 :: specig
  real*8 :: opacleakllump, opacleakrlump
  real*8 :: opacleakdlump, opacleakulump
  real*8 :: mfphelp, ppl, ppr
  real*8 :: resopacleakl, resopacleakr
  real*8 :: resopacleakd, resopacleaku
  integer :: glump, gunlump
  integer :: glumps(gas_ng)
  real*8 :: glumpinv,dtinv,capinv(gas_ng)
  real*8 :: help

  integer,pointer :: zr,zz
  real*8,pointer :: r,z,xi,om,ep,ep0,wl
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dy
  dx(l) = gas_xarr(l+1) - gas_xarr(l)
  dx2(l)= gas_xarr(l+1)**2-gas_xarr(l)**2
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

!--------------------------------------------------------------
!
!-- shortcut
  dtinv = 1d0/tsp_dt
  capinv = 1d0/gas_cap(:,zr,zz,1)

!
!-- set expansion helper
  if(gas_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif
!
!-- looking up initial group
  g = binsrch(wl,gas_wl,gas_ng+1,in_ng)
!-- checking group bounds
  if(g>gas_ng.or.g<1) then
     if(g==gas_ng+1) then
        g = gas_ng
     elseif(g==0) then
        g = 1
     else
        stop 'diffusion2 (1): particle group invalid'
     endif
  endif

!-- lump testing ---------------------------------------------
!
  glump = 0
  gunlump = gas_ng
  glumps = 0
!
!-- find lumpable groups
  if(gas_cap(g,zr,zz,1)*min(dx(zr),dy(zz)) * &
       thelp>=prt_taulump) then
     do ig = 1, g-1
        if(gas_cap(ig,zr,zz,1)*min(dx(zr),dy(zz)) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=ig
        else
           glumps(gunlump)=ig
           gunlump=gunlump-1
        endif
     enddo
     do ig = g, gas_ng
        if(gas_cap(ig,zr,zz,1)*min(dx(zr),dy(zz)) &
             *thelp >= prt_taulump) then
           glump=glump+1
           glumps(glump)=ig
        else
           glumps(gunlump)=ig
           gunlump=gunlump-1
        endif
     enddo
  endif
!
  if(glump==0) then
     glump=1
     glumps(1)=g
!
     forall(ig=2:g) glumps(ig)=ig-1
     forall(ig=g+1:gas_ng) glumps(ig)=ig
!
  endif
  glumpinv = 1d0/glump
!
!-- lumping
  speclump = 0d0
  do ig = 1, glump
     iig = glumps(ig)
     specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1)*capinv(iig)
     speclump = speclump+specig
  enddo
  if(speclump>0d0.and.glump>1) then
     speclump = 1d0/speclump
  else
     speclump = 0d0
  endif

  emitlump = 0d0
  caplump = 0d0
  if(speclump>0d0) then
!
!-- calculating lumped values
     do ig = 1, glump
        iig = glumps(ig)
        specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1)*capinv(iig)
!-- emission lump
        emitlump = emitlump+gas_emitprob(iig,zr,zz,1)
!-- Planck x-section lump
        caplump = caplump+specig*gas_cap(iig,zr,zz,1)*speclump
     enddo
!-- leakage opacities
     opacleakllump = gas_opacleak(1,zr,zz,1)
     opacleakrlump = gas_opacleak(2,zr,zz,1)
     opacleakdlump = gas_opacleak(3,zr,zz,1)
     opacleakulump = gas_opacleak(4,zr,zz,1)
  else
!
!-- calculating unlumped values
     emitlump = gas_emitprob(g,zr,zz,1)
     caplump = gas_cap(g,zr,zz,1)

!-- inward
     if(zr==1) then
        opacleakllump = 0d0
     elseif((gas_cap(g,zr-1,zz,1)+ &
          gas_sig(zr-1,zz,1))*min(dx(zr-1),dy(zz))* &
          thelp<prt_tauddmc) then
!-- DDMC interface
        mfphelp = (gas_cap(g,zr,zz,1)+gas_sig(zr,zz,1))*dx(zr)*thelp
        ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakllump= ppl*(thelp*gas_xarr(zr))/ &
             (thelp**2*dx2(zr))
     else
!-- DDMC interior
        mfphelp = ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1))*dx(zr)+&
             (gas_sig(zr-1,zz,1)+gas_cap(g,zr-1,zz,1))*dx(zr-1))*thelp
        opacleakllump=(4d0/3d0)*(thelp*gas_xarr(zr))/ &
             (mfphelp*thelp**2*dx2(zr))
     endif

!-- outward
     if(zr==gas_nx) then
        lhelp = .true.
     else
        lhelp = (gas_cap(g,zr+1,zz,1)+ &
           gas_sig(zr+1,zz,1))*min(dx(zr+1),dy(zz))* &
           thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        mfphelp = (gas_cap(g,zr,zz,1)+gas_sig(zr,zz,1))*dx(zr)*thelp
        ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
        opacleakrlump= ppr*(thelp*gas_xarr(zr+1))/ &
             (thelp**2*dx2(zr))
     else
!-- DDMC interior
        mfphelp = ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1))*dx(zr)+&
             (gas_sig(zr+1,zz,1)+gas_cap(g,zr+1,zz,1))*dx(zr+1))*thelp
        opacleakrlump=(4d0/3d0)*(thelp*gas_xarr(zr+1))/ &
             (mfphelp*thelp**2*dx2(zr))
     endif

!-- downward
     if(zz==1) then
        lhelp = .true.
     else
        lhelp = (gas_cap(g,zr,zz-1,1)+ &
             gas_sig(zr,zz-1,1))*min(dx(zr),dy(zz-1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (gas_cap(g,zr,zz,1)+gas_sig(zr,zz,1))*dy(zz)*thelp
        ppl = 4d0/(3d0*help+6d0*pc_dext)
        opacleakdlump=0.5d0*ppl/(thelp*dy(zz))
     else
!-- DDMC interior
        help = ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1))*dy(zz)+&
             (gas_sig(zr,zz-1,1)+gas_cap(g,zr,zz-1,1))*dy(zz-1))*thelp
        opacleakdlump=(2d0/3d0)/(help*dy(zz)*thelp)
     endif

!-- upward
     if(zz==gas_ny) then
        lhelp = .true.
     else
        lhelp = (gas_cap(g,zr,zz+1,1)+ &
             gas_sig(zr,zz+1,1))*min(dx(zr),dy(zz+1)) * &
             thelp<prt_tauddmc
     endif
     if(lhelp) then
!-- DDMC interface
        help = (gas_cap(g,zr,zz,1)+gas_sig(zr,zz,1))*dy(zz)*thelp
        ppr = 4d0/(3d0*help+6d0*pc_dext)
        opacleakulump=0.5d0*ppr/(thelp*dy(zz))
     else
!-- DDMC interior
        help = ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1))*dy(zz)+&
             (gas_sig(zr,zz+1,1)+gas_cap(g,zr,zz+1,1))*dy(zz+1))*thelp
        opacleakulump=(2d0/3d0)/(help*dy(zz)*thelp)
     endif
  endif
!
!-------------------------------------------------------------
!

!-- calculating time to census or event
  denom = opacleakllump+opacleakrlump + &
       opacleakdlump+opacleakulump + &
       (1d0-emitlump)*(1d0-gas_fcoef(zr,zz,1))*caplump
  if(prt_isddmcanlog) then
     denom = denom+gas_fcoef(zr,zz,1)*caplump
  endif

  r1 = rand()
  tau = abs(log(r1)/(pc_c*denom))
  tcensus = tsp_t+tsp_dt-ptcl%tsrc
  ddmct = min(tau,tcensus)

!
!-- calculating energy depostion and density
  !
  if(prt_isddmcanlog) then
     gas_eraddens(zr,zz,1)= gas_eraddens(zr,zz,1)+ep*ddmct*dtinv
  else
     gas_edep(zr,zz,1) = gas_edep(zr,zz,1)+ep*(1d0-exp(-gas_fcoef(zr,zz,1) &!{{{
          *caplump*pc_c*ddmct))
     if(gas_fcoef(zr,zz,1)*caplump*min(dx(zr),dy(zz)) * &
          thelp>1d-6) then
        help = 1d0/(gas_fcoef(zr,zz,1)*caplump)
        gas_eraddens(zr,zz,1)= &
             gas_eraddens(zr,zz,1)+ep* &
             (1d0-exp(-gas_fcoef(zr,zz,1)*caplump*pc_c*ddmct))* &
             help*cinv*dtinv
     else
        gas_eraddens(zr,zz,1) = gas_eraddens(zr,zz,1)+ep*ddmct*dtinv
     endif
     ep=ep*exp(-gas_fcoef(zr,zz,1)*caplump*pc_c*ddmct)
!!}}}
  endif

!-- updating particle time
  ptcl%tsrc = ptcl%tsrc+ddmct

!-- stepping particle ------------------------------------
!
!
!-- check for census
  if (ddmct /= tau) then
     prt_done = .true.
     gas_numcensus(zr,zz,1)=gas_numcensus(zr,zz,1)+1

!-- otherwise, perform event
  else
     r1 = rand()!{{{
     help = 1d0/denom

!-- outward leak probability
     PR = opacleakrlump*help
!-- inward leak probability
     PL = opacleakllump*help
!-- upward leak probability
     PU = opacleakulump*help
!-- downward leak probability
     PD = opacleakdlump*help
!-- absorption probability
     if(prt_isddmcanlog) then
        PA = gas_fcoef(zr,zz,1)*caplump*help
     else
        PA = 0d0
     endif

!-- absorption sample
     if(r1>=0d0 .and. r1<PA) then
        isvacant = .true.
        prt_done = .true.
        gas_edep(zr,zz,1) = gas_edep(zr,zz,1)+ep

!-- inward leakage sample
     elseif (r1>=PA .and. r1<PA+PL) then
!{{{
!-- checking if at inner bound
        if (zr == 1) then
           stop 'diffusion1: non-physical inward leakage'

!-- sample adjacent group (assumes aligned g bounds)
        else
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakllump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,zr-1,zz,1) + &
                      gas_sig(zr-1,zz,1))*min(dx(zr-1),dy(zz)) * &
                      thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                         dx(zr)*thelp
                    ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakl = ppl*(thelp*gas_xarr(zr))/ &
                         (thelp**2*dx2(zr))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(zr,zz,1)+gas_cap(iig,zr,zz,1)) * &
                         dx(zr)+(gas_sig(zr-1,zz,1) + &
                         gas_cap(iig,zr-1,zz,1))*dx(zr-1))*thelp
                    resopacleakl = (4d0/3d0)*(thelp*gas_xarr(zr))/ &
                         (mfphelp*thelp**2*dx2(zr))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakl*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakl*speclump*help
              enddo
           else
              iig = g
           endif

!-- updating properties
           if((gas_sig(zr-1,zz,1)+gas_cap(iig,zr-1,zz,1)) * &
                min(dx(zr-1),dy(zz))*thelp>=prt_tauddmc) then
              zr = zr-1
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
              g = iig
           else
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
!
!-- method changed to IMC
              ptcl%rtsrc = 1
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              z = gas_yarr(zz)*(1d0-r1)+gas_yarr(zz+1)*r1
              r = gas_xarr(zr)

!-- current particle cell set to 1 left
              zr = zr-1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu0 = -max(r1,r2)
              r1 = rand()
              xi = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(abs(1d0-xi**2-mu0**2)),mu0)
              if(om<0d0) then
                 om = om+pc_pi2
              endif
!-- extending azimuthal sample
              r1 = rand()
              if(r1 < 0.5d0) then
                 om = pc_pi2-om
              endif

!-- doppler and aberration corrections
              if(gas_isvelocity) then
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
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 ep = ep*(1d0+dirdotu*cinv)
                 ep0 = ep0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif

        endif!}}}

!
!-- outward leakage sample
     elseif (r1>=PA+PL .and. r1<PA+PL+PR) then
!!{{{
!-- checking if at outer bound
        if (zr == gas_nx) then
           isvacant = .true.
           prt_done = .true.
           gas_eout = gas_eout+ep
!-- outbound luminosity tally
!-- sampling r, z
           r1 = rand()
           z = gas_yarr(zz)*(1d0-r1)+gas_yarr(zz+1)*r1
           r = gas_xarr(gas_nx+1)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           mu0 = max(r1,r2)
           r1 = rand()
           xi = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
           om = atan2(sqrt(abs(1d0-xi**2-mu0**2)),mu0)
           if(om<0d0) then
              om = om+pc_pi2
           endif
!-- extending azimuthal sample
           r1 = rand()
           if(r1 < 0.5d0) then
              om = pc_pi2-om
           endif
           if(gas_isvelocity) then
              dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
              elabfact = 1d0+dirdotu*cinv
           endif
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakrlump
              do ig = 1, glump
                 iig=glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                      dx(zr)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleakr = ppr*(thelp*gas_xarr(zr+1))/ &
                      (thelp**2*dx2(zr))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakr*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakr*speclump*help
              enddo

              r1 = rand()
              wl=1d0/(r1/gas_wl(iig+1) + (1d0-r1)/gas_wl(iig))
!-- changing from comoving frame to observer frame
              wl = wl/elabfact
!-- obtaining spectrum (lab) group
              iig = binsrch(wl,gas_wl,gas_ng+1,in_ng)
              if(iig>gas_ng.or.iig<1) then
                 if(iig>gas_ng) then
                    iig=gas_ng
                    wl=gas_wl(gas_ng+1)
                 else
                    iig=1
                    wl=gas_wl(1)
                 endif
              endif
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           else
              iig = g
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           endif

        else
!-- outward leaking in domain
!
!-- sample adjacent group (assumes aligned g bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakrlump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,zr+1,zz,1)+gas_sig(zr+1,zz,1)) * &
                      min(dx(zr+1),dy(zz))*thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                         dx(zr)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakr = ppr*(thelp*gas_xarr(zr+1))/ &
                         (thelp**2*dx2(zr))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(zr,zz,1)+gas_cap(iig,zr,zz,1)) * &
                         dx(zr)+&
                         (gas_sig(zr+1,zz,1)+gas_cap(iig,zr+1,zz,1)) * &
                         dx(zr+1))*thelp
                    resopacleakr = (4d0/3d0)*(thelp*gas_xarr(zr+1))/ &
                         (mfphelp*thelp**2*dx2(zr))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakr*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakr*speclump*help
              enddo
           else
              iig = g
           endif

!-- updating particle properties
           if((gas_sig(zr+1,zz,1)+gas_cap(iig,zr+1,zz,1)) * &
                min(dx(zr+1),dy(zz)) &
                *thelp >= prt_tauddmc) then
              zr = zr+1
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
              g = iig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
!
!-- method changed to IMC
              ptcl%rtsrc = 1
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              z = gas_yarr(zz)*(1d0-r1)+gas_yarr(zz+1)*r1
              r = gas_xarr(zr+1)

!-- current particle cell set to 1 left
              zr = zr+1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              mu0 = max(r1,r2)
              r1 = rand()
              xi = sqrt(1d0-mu0**2)*cos(pc_pi2*r1)
              om = atan2(sqrt(abs(1d0-xi**2-mu0**2)),mu0)
              if(om<0d0) then
                 om = om+pc_pi2
              endif
!-- extending azimuthal sample
              r1 = rand()
              if(r1 < 0.5d0) then
                 om = pc_pi2-om
              endif

!-- doppler and aberration corrections
              if(gas_isvelocity) then
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
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 ep = ep*(1d0+dirdotu*cinv)
                 ep0 = ep0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif
        endif

!
!-- downward leakage sample
     elseif(r1>=PA+PL+PR.and.r1<PA+PL+PR+PD) then

!-- checking if at outer bound
        if (zz == 1) then
           isvacant = .true.
           prt_done = .true.
           gas_eout = gas_eout+ep
!-- outbound luminosity tally
!-- sampling r, z
           r1 = rand()
           r = sqrt(gas_xarr(zr)**2*(1d0-r1)+gas_xarr(zr+1)**2*r1)
           z = gas_yarr(zz)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           xi = -max(r1,r2)
           r1 = rand()
           om = pc_pi2*r1
           if(gas_isvelocity) then
              dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
              elabfact = 1d0+dirdotu*cinv
           endif
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakdlump
              do ig = 1, glump
                 iig=glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                      dy(zz)*thelp
                 ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleakd = 0.5d0*ppl/(thelp*dy(zz))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakd*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakd*speclump*help
              enddo

              r1 = rand()
              wl=1d0/(r1/gas_wl(iig+1) + (1d0-r1)/gas_wl(iig))
!-- changing from comoving frame to observer frame
              wl = wl/elabfact
!-- obtaining spectrum (lab) group
              iig = binsrch(wl,gas_wl,gas_ng+1,in_ng)
              if(iig>gas_ng.or.iig<1) then
                 if(iig>gas_ng) then
                    iig=gas_ng
                    wl=gas_wl(gas_ng+1)
                 else
                    iig=1
                    wl=gas_wl(1)
                 endif
              endif
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           else
              iig = g
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           endif
        else

!-- sample adjacent group (assumes aligned g bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakdlump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,zr,zz-1,1)+gas_sig(zr,zz-1,1)) * &
                      min(dx(zr),dy(zz-1))*thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                         dy(zz)*thelp
                    ppl = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleakd = 0.5d0*ppl/(thelp*dy(zz))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(zr,zz,1)+gas_cap(iig,zr,zz,1)) * &
                         dy(zz)+&
                         (gas_sig(zr,zz-1,1)+gas_cap(iig,zr,zz-1,1)) * &
                         dy(zz-1))*thelp
                    resopacleakd = (2d0/3d0)/(mfphelp*thelp*dy(zz))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleakd*speclump*help)) exit
                 denom2 = denom2+specig*resopacleakd*speclump*help
              enddo
           else
              iig = g
           endif

!-- updating particle properties
           if((gas_sig(zr,zz-1,1)+gas_cap(iig,zr,zz-1,1)) * &
                min(dx(zr),dy(zz-1)) &
                *thelp >= prt_tauddmc) then
              zz = zz-1
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
              g = iig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
!
!-- method changed to IMC
              ptcl%rtsrc = 1
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              r = sqrt(gas_xarr(zr)**2*(1d0-r1)+gas_xarr(zr+1)**2*r1)
              z = gas_yarr(zz)

!-- current particle cell set to 1 left
              zz = zz-1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              xi = -max(r1,r2)
              r1 = rand()
              om = pc_pi2*r1

!-- doppler and aberration corrections
              if(gas_isvelocity) then
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
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 ep = ep*(1d0+dirdotu*cinv)
                 ep0 = ep0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif
        endif

!
!-- upward leakage sample
     elseif(r1>=PA+PL+PR+PD.and.r1<PA+PL+PR+PD+PU) then

!-- checking if at outer bound
        if (zz == gas_ny) then
           isvacant = .true.
           prt_done = .true.
           gas_eout = gas_eout+ep
!-- outbound luminosity tally
!-- sampling r, z
           r1 = rand()
           r = sqrt(gas_xarr(zr)**2*(1d0-r1)+gas_xarr(zr+1)**2*r1)
           z = gas_yarr(zz+1)
!-- sampling direction
           r1 = rand()
           r2 = rand()
           xi = max(r1,r2)
           r1 = rand()
           om = pc_pi2*r1
           if(gas_isvelocity) then
              dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
              elabfact = 1d0+dirdotu*cinv
           endif
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakulump
              do ig = 1, glump
                 iig=glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                      dy(zz)*thelp
                 ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                 resopacleaku = 0.5d0*ppr/(thelp*dy(zz))
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleaku*speclump*help)) exit
                 denom2 = denom2+specig*resopacleaku*speclump*help
              enddo

              r1 = rand()
              wl=1d0/(r1/gas_wl(iig+1) + (1d0-r1)/gas_wl(iig))
!-- changing from comoving frame to observer frame
              wl = wl/elabfact
!-- obtaining spectrum (lab) group
              iig = binsrch(wl,gas_wl,gas_ng+1,in_ng)
              if(iig>gas_ng.or.iig<1) then
                 if(iig>gas_ng) then
                    iig=gas_ng
                    wl=gas_wl(gas_ng+1)
                 else
                    iig=1
                    wl=gas_wl(1)
                 endif
              endif
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           else
              iig = g
              gas_luminos(iig) = gas_luminos(iig)+&
                   ep*dtinv*elabfact
              gas_lumdev(iig) = gas_lumdev(iig)+&
                   (ep*dtinv*elabfact)**2
              gas_lumnum(iig) = gas_lumnum(iig)+1
           endif
        else

!-- sample adjacent group (assumes aligned g bounds)
           if(speclump>0d0) then
              r1 = rand()
              denom2 = 0d0
              help = 1d0/opacleakulump
              do ig= 1, glump
                 iig = glumps(ig)
                 specig = gas_siggrey(zr,zz,1)*gas_emitprob(iig,zr,zz,1) * &
                      capinv(iig)
!-- calculating resolved leakage opacities
                 if((gas_cap(iig,zr,zz+1,1)+gas_sig(zr,zz+1,1)) * &
                      min(dx(zr),dy(zz+1))*thelp<prt_tauddmc) then
!-- DDMC interface
                    mfphelp = (gas_cap(iig,zr,zz,1)+gas_sig(zr,zz,1)) * &
                         dy(zz)*thelp
                    ppr = 4d0/(3d0*mfphelp+6d0*pc_dext)
                    resopacleaku = 0.5d0*ppr/(thelp*dy(zz))
                 else
!-- IMC interface
                    mfphelp = ((gas_sig(zr,zz,1)+gas_cap(iig,zr,zz,1)) * &
                         dy(zz)+&
                         (gas_sig(zr,zz+1,1)+gas_cap(iig,zr,zz+1,1)) * &
                         dy(zz+1))*thelp
                    resopacleaku = (2d0/3d0)/(mfphelp*thelp*dy(zz))
                 endif
                 if((r1>=denom2).and. &
                      (r1<denom2+specig*resopacleaku*speclump*help)) exit
                 denom2 = denom2+specig*resopacleaku*speclump*help
              enddo
           else
              iig = g
           endif

!-- updating particle properties
           if((gas_sig(zr,zz+1,1)+gas_cap(iig,zr,zz+1,1)) * &
                min(dx(zr),dy(zz+1)) &
                *thelp >= prt_tauddmc) then
              zz = zz+1
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
              g = iig
           else
!-- converting to IMC
              r1 = rand()
              wl = 1d0/(r1/gas_wl(iig+1)+(1d0-r1)/gas_wl(iig))
!
!-- method changed to IMC
              ptcl%rtsrc = 1
              gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!
!-- location set right bound of left cell
              r1 = rand()
              r = sqrt(gas_xarr(zr)**2*(1d0-r1)+gas_xarr(zr+1)**2*r1)
              z = gas_yarr(zz+1)

!-- current particle cell set to 1 left
              zz = zz+1
!-- particl angle sampled from isotropic b.c. inward
              r1 = rand()
              r2 = rand()
              xi = max(r1,r2)
              r1 = rand()
              om = pc_pi2*r1

!-- doppler and aberration corrections
              if(gas_isvelocity) then
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
!-- transforming wavelength to lab
                 wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
                 ep = ep*(1d0+dirdotu*cinv)
                 ep0 = ep0*(1d0+dirdotu*cinv)
              endif
!
!-- group reset
              g = iig
!
           endif
        endif

!-- effective scattering
     else
        denom2 = 1d0-emitlump
        help = 1d0/denom2
!
        denom3 = 0d0
        r1 = rand()

        do ig = gas_ng,glump+1,-1
           iig=glumps(ig)
           if((r1>=denom3).and.(r1<denom3+gas_emitprob(iig,zr,zz,1)*help)) exit
           denom3 = denom3+gas_emitprob(iig,zr,zz,1)*help
        enddo
!
        g = iig
        r1 = rand()
        wl = 1d0/((1d0-r1)/gas_wl(g) + r1/gas_wl(g+1))

        if ((gas_sig(zr,zz,1)+gas_cap(g,zr,zz,1)) * &
             min(dx(zr),dy(zz)) &
             *thelp >= prt_tauddmc) then
           ptcl%rtsrc = 2
        else
           ptcl%rtsrc = 1
           gas_methodswap(zr,zz,1)=gas_methodswap(zr,zz,1)+1
!-- direction sampled isotropically           
           r1 = rand()
           xi = 1d0 - 2d0*r1
           r1 = rand()
           om = pc_pi2*r1
!-- position sampled uniformly
           r1 = rand()
           r = sqrt(r1*gas_xarr(zr+1)**2+(1d0-r1)*gas_xarr(zr)**2)
           r1 = rand()
           z = r1*gas_yarr(zz+1)+(1d0-r1)*gas_yarr(zz)

!-- doppler and aberration corrections
           if(gas_isvelocity) then
!-- calculating transformation factors
              dirdotu = xi*z+sqrt(1d0-xi**2)*cos(om)*r
              azidotu = atan2(sqrt(1d0-xi**2)*sin(om), &
                   sqrt(1d0-xi**2)*cos(om)+r*cinv)
!-- transforming z-axis direction cosine to lab
              xi = (xi+z/pc_c)/(1d0+dirdotu/pc_c)
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
!-- transforming wavelength to lab
              wl = wl/(1d0+dirdotu*cinv)
!-- transforming energy weights to lab
              ep = ep*(1d0+dirdotu*cinv)
              ep0 = ep0*(1d0+dirdotu*cinv)
           endif
        endif
     endif
  endif!}}}


end subroutine diffusion2
