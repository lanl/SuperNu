subroutine transport1_gamgrey(ptcl,isvacant)

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
  integer :: binsrch
  real*8 :: r1, r2, thelp,thelpinv
  real*8 :: db, dcol, dcen, dthm, ddop, d
  real*8 :: siglabfact, dcollabfact, elabfact
  real*8 :: rold, P, told, muold
! real*8 :: x1, x2, xx0
  real*8 :: dtinv
  real*8 :: help
  real*8 :: ppl, ppr

  integer,pointer :: z
  real*8,pointer :: r, mu, E, E0
!-- statement function
  integer :: l
  real*8 :: dx
  dx(l) = grd_xarr(l+1) - grd_xarr(l)

  z => ptcl%zsrc
  r => ptcl%rsrc
  mu => ptcl%musrc
  E => ptcl%esrc
  E0 => ptcl%ebirth
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
!== DISTANCE CALCULATIONS
!
!-- distance to boundary = db
  if (z == 1) then
     db = abs(sqrt(grd_xarr(z+1)**2-(1.0-mu**2)*r**2)-mu*r)
  elseif (mu < -sqrt(1.0d0-(grd_xarr(z)/r)**2)) then
     db = abs(sqrt(grd_xarr(z)**2-(1.0d0-mu**2)*r**2)+mu*r)
  else
     db = abs(sqrt(grd_xarr(z+1)**2-(1.0d0-mu**2)*r**2)-mu*r)
  endif
!
!-- distance to fictitious collision = dcol
  if(prt_isimcanlog) then
     if(gas_capgam(z,1,1)>0d0) then
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        dcol = abs(log(r1)/(gas_capgam(z,1,1)*dcollabfact))
     else
        dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
     endif
  else
     dcol = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- distance to Thomson-type collision = dthm
  dthm = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
!
!-- distance to census = dcen
  dcen = abs(pc_c*(tsp_t+tsp_dt-ptcl%tsrc)*thelpinv)
!
!-- distance to Doppler shift = ddop
  if(grd_isvelocity) then
     ddop = abs(pc_c - r*mu)
  else
     ddop = 2d0*abs(pc_c*tsp_dt*thelpinv) !> dcen
  endif
!
!-- minimum distance = d
!  if(tsp_it==29) write(*,*) dcol,dthm,db,dcen,ddop
  d = min(dcol,dthm,db,dcen,ddop)
!
!== END OF DISTANCE CALCULATIONS
!
!-- position, angle, time update  
  rold = r
  r = sqrt((1.0d0-mu**2)*r**2+(d+r*mu)**2)
!  r = sqrt(r**2+d**2+2d0*d*r*mu)
  told = ptcl%tsrc
  ptcl%tsrc = ptcl%tsrc + thelp*d*cinv
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
        grd_edep(z,1,1)=grd_edep(z,1,1)+E*(1.0d0-exp( &
             -gas_capgam(z,1,1)*siglabfact*d*thelp))*elabfact
     !--
     if(gas_capgam(z,1,1)*dx(z)*thelp>1d-6) then     
        grd_eraddens(z,1,1) = grd_eraddens(z,1,1)+E* &
             (1.0d0-exp(-siglabfact*gas_capgam(z,1,1)*d*thelp))* &
             elabfact/(siglabfact*gas_capgam(z,1,1)*pc_c*tsp_dt)
     else
        grd_eraddens(z,1,1) = grd_eraddens(z,1,1)+E* &
             elabfact*d*dcollabfact*cinv*dtinv
     endif
     !--
     E = E*exp(-gas_capgam(z,1,1)*siglabfact*d*thelp)

  else
     !
     grd_eraddens(z,1,1) = grd_eraddens(z,1,1)+E* &
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
!
        E = E*elabfact*help
!        E0 = E0*elabfact/(1.0-mu*r*cinv)
     endif
     !
     !!}}}
  elseif (d == dcol) then  !fictitious scattering with implicit capture
     !!{{{
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     if(r1<=1d0.and.prt_isimcanlog) then
        isvacant = .true.
        prt_done = .true.
        grd_edep(z,1,1) = grd_edep(z,1,1) + E*elabfact
!-- velocity effects accounting
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
!
           E = E*elabfact*help
           
        endif
!
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
     endif
     !!}}}
  elseif (d == db) then   !------boundary crossing ----
     if (mu>=0.0d0) then!{{{
        if (z == grd_nx) then
           isvacant = .true.
           prt_done = .true.
!
!-- outbound luminosity tally
!-- velocity effects accounting
!
           tot_eright = tot_eright+E*elabfact
           flx_gamluminos(1,1) = flx_gamluminos(1,1)+E*dtinv
           flx_gamlumdev(1,1) = flx_gamlumdev(1,1)+(E*dtinv)**2
           flx_gamlumnum(1,1) = flx_gamlumnum(1,1)+1
        ! Checking if DDMC region right
           else
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = -max(r1,r2)
              if(grd_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
           endif
        ! End of check
        else
           z = z+1
           r = grd_xarr(z)
        endif
     else
        if (z==1) then
           if ((gas_capgam(z,1,1)*dx(z+1) &
                *thelp >= prt_tauddmc) &
                .and.(in_puretran.eqv..false.)) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(grd_isvelocity) then
                 mu = (mu-r*cinv)/(1.0-r*mu*cinv)
              endif
              help = gas_capgam(z,1,1)*dx(z+1)*thelp
              ppl = 4d0/(3d0*help+6d0*pc_dext)
              P = ppl*(1.0+1.5*abs(mu))

              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
              mu = -max(r1,r2)
              if(grd_isvelocity) then
                 mu = (mu+r*cinv)/(1.0+r*mu*cinv)
              endif
           else
              z = z+1
           endif
        elseif ((gas_capgam(z,1,1)*dx(z-1) &
             *thelp >= prt_tauddmc) &
             .and.(in_puretran.eqv..false.)) then
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           if(grd_isvelocity) then
              mu = (mu-r*cinv)/(1.0-r*mu*cinv)
!-- amplification
!
!              E0=E0*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
!              E = E*(1d0+2d0*min(0.055*prt_tauddmc,1d0)*r*cinv)
              if(mu<0d0) then
!-- velocity effects accounting
                 help = 1d0/abs(mu)
!
                 E0 = E0*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
                 E = E*(1d0+2d0*(0.55d0*help-1.25d0*abs(mu))*r*cinv)
              endif
!--
           endif
           help = gas_capgam(z,1,1)*dx(z-1)*thelp
           ppr = 4d0/(3d0*help+6d0*pc_dext)
           P = ppr*(1.0+1.5*abs(mu))
!--
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           r2 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu = max(r1,r2)
           if(grd_isvelocity) then
              mu = (mu+r*cinv)/(1.0+r*mu*cinv)
           endif
        else
           z = z-1
        endif
     endif!}}}
  elseif (d == dcen) then
     prt_done = .true.
     grd_numcensus(z,1,1) = grd_numcensus(z,1,1)+1
!     tot_erad = tot_erad + E*elabfact
!
  endif

end subroutine transport1_gamgrey
