subroutine particle_advance_gamgrey(nmpi)

  use particlemod
  use gridmod
  use physconstmod
  use inputparmod
  use timingmod
  use fluxmod
  implicit none
  integer,intent(in) :: nmpi
!##################################################
  !This subroutine propagates all existing particles that are not vacant
  !during a time step.  Particles may generally undergo a physical interaction
  !with the gas, cross a spatial cell boundary, or be censused for continued
  !propagation in the next time step.  Currently DDMC and IMC particle events
  !are being handled in separate subroutines but this may be changed to reduce
  !total subroutine calls in program.
!##################################################
  integer :: nhere, nemit, ndmy
  real*8 :: r1
  integer :: i, j, k, ii, iimpi
  integer,pointer :: ix, iy, iz
  real*8,pointer :: x,y,z,mu,om,e,e0
  real*8 :: t0,t1  !timing
  real*8 :: labfact, cmffact, mu1, mu2
  real*8 :: etot,pwr
  real*8 :: om0, mu0, x0, y0, z0
!
  real*8,parameter :: basefrac=.1d0
  real*8 :: base,edone,einv,invn,en
  integer :: n,ndone
  integer*8 :: nstot,nsavail,nsbase
!-- hardware
!
  type(packet),target :: ptcl

  grd_edep = 0d0
  flx_gamluminos = 0d0
  flx_gamlumdev = 0d0
  flx_gamlumnum = 0
!
  grd_eraddens = 0d0
  grd_numcensus = 0

!-- initializing volume numbers
  grd_nvol=0

!-- shortcut
  pwr = in_srcepwr

!-- total particle number
  nstot = nmpi*int(prt_ns,8)

!-- total energy (converted by pwr)
  etot = sum(grd_emitex**pwr)

!-- base (flat,constant) particle number per cell over ALL RANKS
  n = count(grd_emitex>0d0)  !number of cells that emit
  base = dble(nstot)/n  !uniform distribution
  base = basefrac*base

!-- number of particles available for proportional distribution
  nsbase = int(n*base,8)  !total number of base particles
  nsavail = nstot - nsbase

!-- total particle number per cell
  edone = 0d0
  ndone = 0
  invn = 1d0/nstot
  einv = 1d0/etot
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     en = grd_emitex(i,j,k)**pwr
     if(en==0d0) cycle
!-- continuously guide the rounding towards the correct cumulative value
     n = int(en*nsavail*einv + base)  !round down
     if(edone*einv>ndone*invn) n = n + 1  !round up
     grd_nvol(i,j,k) = n
     edone = edone + en
     ndone = ndone + n
  enddo
  enddo
  enddo

!--
  ix => ptcl%ix
  iy => ptcl%iy
  iz => ptcl%iz
  x => ptcl%x
  y => ptcl%y
  z => ptcl%z
  mu => ptcl%mu
  om => ptcl%om
  e => ptcl%e
  e0 => ptcl%e0
!-- unused
!    real*8 :: tsrc
!    real*8 :: wlsrc
!    integer :: rtsrc

  call time(t0)

  iimpi = 0
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     call sourcenumbers_roundrobin(iimpi,grd_emitex(i,j,k), &
        0d0,grd_nvol(i,j,k),nemit,nhere,ndmy)
  do ii=1,nhere
!-- adopt position
     ix = i
     iy = j
     iz = k

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- particle propagation
     select case(in_igeom)
     case(1)
!-- calculating position!{{{
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        x = (r1*grd_xarr(ix+1)**3 + &
             (1.0-r1)*grd_xarr(ix)**3)**(1.0/3.0)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!--
        if(grd_isvelocity) then
           x0 = x
           cmffact = 1d0+mu0*x0/pc_c !-- 1+dir*v/c
           mu = (mu0+x0/pc_c)/cmffact
        else
           mu = mu0
        endif!}}}
     case(2)
!-- calculating position!{{{
        r1 = rand()
        x = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
!-- must be inside cell
        x = min(x,grd_xarr(ix+1))
        x = max(x,grd_xarr(ix))
!
        r1 = rand()
        y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = x
           y0 = y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
!-- om
           om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
           if(om<0d0) om=om+pc_pi2
!-- mu
           mu = (mu0+y0/pc_c)/cmffact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
        else
           mu = mu0
           om = om0
        endif!}}}
     case(3)
!-- setting 2nd,3rd cell index!{{{
        iy = j
        iz = k
!-- calculating position
        r1 = rand()
        x = r1*grd_xarr(i+1) + (1d0-r1) * &
             grd_xarr(i)
        r1 = rand()
        y = r1*grd_yarr(j+1) + (1d0-r1) * &
             grd_yarr(j)
        r1 = rand()
        z = r1*grd_zarr(k+1) + (1d0-r1) * &
             grd_zarr(k)
!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1
!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = x
           y0 = y
           z0 = z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           mu = (mu0+z0/pc_c)/cmffact
           if(mu>1d0) then
              mu = 1d0
           elseif(mu<-1d0) then
              mu = -1d0
           endif
!-- om
           om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(om<0d0) om = om+pc_pi2
        else
           mu = mu0
           om = om0
        endif!}}}
     endselect
!
!-- emission energy per particle
     e = grd_emitex(ix,iy,iz)/nemit*cmffact
     e0 = e

!-----------------------------------------------------------------------
!-- Advancing particle until census, absorption, or escape from domain
     prt_done=.false.
!
     select case(in_igeom)
     case(1)
        do while (.not.prt_done)!{{{
           call transport1_gamgrey(ptcl)
!-- verify position
           if(.not.prt_done .and. (x>grd_xarr(ix+1) .or. x<grd_xarr(ix))) then
              write(0,*) 'prt_adv_ggrey: not in cell', &
                 ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
           endif
!-- transformation factor
           if(grd_isvelocity) then
              labfact = 1.0d0 - mu*x/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (e<1d-6*e0 .and. .not.prt_done) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
              else
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo!}}}
     case(2)
        do while (.not.prt_done)!{{{
           call transport2_gamgrey(ptcl)
!-- verify position
           if(.not.prt_done) then
              if(x>grd_xarr(ix+1) .or. x<grd_xarr(ix)) then
                write(0,*) 'prt_adv_ggrey: r not in cell', &
                   ix,x,grd_xarr(ix),grd_xarr(ix+1),mu
              endif
              if(y>grd_yarr(iy+1) .or. y<grd_yarr(iy)) then
                write(0,*) 'prt_adv_ggrey: r not in cell', &
                   iy,y,grd_yarr(iy),grd_yarr(iy+1),mu
              endif
           endif
!-- transformation factor
           if(grd_isvelocity) then
              labfact = 1d0-(mu*y+sqrt(1d0-mu**2) * &
                   cos(om)*x)/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (e<1d-6*e0 .and. .not.prt_done) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
              else
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo!}}}
     case(3)
        do while (.not.prt_done)!{{{
           call transport3_gamgrey(ptcl)
!-- transformation factor
           if(grd_isvelocity) then
              labfact = 1d0-(mu*z+sqrt(1d0-mu**2) * &
                   (cos(om)*x+sin(om)*y))/pc_c
           else
              labfact = 1d0
           endif
!-- Russian roulette for termination of exhausted particles
           if (e<1d-6*e0 .and. .not.prt_done) then
              r1 = rand()
              prt_tlyrand = prt_tlyrand+1
              if(r1<0.5d0) then
                 prt_done = .true.
                 grd_edep(ix,iy,iz) = grd_edep(ix,iy,iz) + e*labfact
              else
                 e = 2d0*e
                 e0 = 2d0*e0
              endif
           endif
        enddo!}}}
     endselect

  enddo !ii
  enddo !i
  enddo !j
  enddo !k

  call time(t1)
  call timereg(t_pcktgam, t1-t0)

end subroutine particle_advance_gamgrey
