subroutine initial_particles

  use randommod
  use gridmod
  use groupmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine instantiates the initial particles before
  !the first time step.
!##################################################
  type(packet) :: ptcl
!
  logical :: lhelp
  integer :: ig, i,j,k, iig, ipart
  integer, dimension(grd_nx,grd_ny,grd_nz) :: ijkused
  real*8 :: wl0, mu0, om0, ep0, x0, y0, z0
  real*8 :: denom2, mu1, mu2
  real*8 :: r1
  real*8 :: wl1,wl2,wl3,wl4
  real*8 :: cmffact,gm

!-- helper quantities
  wl1=grp_wlinv(grp_ng+1)
  wl2=grp_wlinv(1)

!-- instantiating initial particles
  i = 1
  j = 1
  k = 1
  ijkused = 0
  do ipart=1,prt_ninitnew

!-- incrementing to next vacant cell
     loopk: do k=k,grd_nz
        do j=j,grd_ny
           do i=i,grd_nx
              lhelp = ijkused(i,j,k)<grd_nvolinit(i,j,k)
              if(lhelp) exit loopk
           enddo
           i = 1
        enddo
        j = 1
     enddo loopk
!
!-- sanity check
     if(.not.lhelp) stop 'initial_particles: invalid particle'

!-- increasing cell occupancy
     ijkused(i,j,k) = ijkused(i,j,k)+1
!
!-- setting 1st cell index
     ptcl%ix = i

!-- setting particle index to not vacant
     prt_isvacant(ipart) = .false.
!
!-- calculating particle time
     ptcl%t = tsp_t
!
!-- setting type
     ptcl%itype = 1

!-- calculating wavelength
     denom2 = 0d0
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grp_ng
        wl3 = grp_wlinv(ig+1)
        wl4 = grp_wlinv(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2 + (wl4-wl3)/(wl2-wl1)) exit
        denom2 = denom2 + (wl4-wl3)/(wl2-wl1)
     enddo
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)*grp_wlinv(iig)+r1*grp_wlinv(iig+1))

!-- calculating direction cosine (comoving)
     r1 = rnd_r(rnd_state)
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = grd_evolinit(i,j,k)/dble(grd_nvolinit(i,j,k))

!
!-- selecting geometry
     select case(in_igeom)

!-- 3D spherical
     case(1)
!-- calculating position
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%y = r1*grd_yarr(j+1)+(1d0-r1)*grd_yarr(j)
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%z = r1*grd_zarr(k+1)+(1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))
!-- sampling azimuthal angle of direction
        r1 = rnd_r(rnd_state)
        ptcl%om = pc_pi2*r1
!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = ptcl%x
!-- 1+dir*v/c
           cmffact = 1d0+x0*mu0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
        else
           ptcl%mu = mu0
        endif

!-- 2D
     case(2)
!-- setting 2nd cell index
        ptcl%iy = j
!-- calculating position
        r1 = rnd_r(rnd_state)
        ptcl%x = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
        r1 = rnd_r(rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + &
             (1d0-r1)*grd_yarr(j)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
!-- sampling azimuthal angle of direction
        r1 = rnd_r(rnd_state)
        om0 = pc_pi2*r1

!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           gm = 1d0/sqrt(1d0-(x0**2+y0**2)/pc_c**2)
!-- om
           ptcl%om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+(gm*x0/pc_c) * &
                (1d0+gm*(cmffact-1d0)/(gm+1d0)))
           if(ptcl%om<0d0) ptcl%om=ptcl%om+pc_pi2
!-- mu
           ptcl%mu = (mu0+(gm*y0/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
                (gm*cmffact)
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif

!-- 3D
     case(3)
!-- setting 2nd, 3rd cell index
        ptcl%iy = j
        ptcl%iz = k
!-- calculating position
        r1 = rnd_r(rnd_state)
        ptcl%x = r1*grd_xarr(i+1) + &
             (1d0-r1)*grd_xarr(i)
        r1 = rnd_r(rnd_state)
        ptcl%y = r1*grd_yarr(j+1) + &
             (1d0-r1)*grd_yarr(j)
        r1 = rnd_r(rnd_state)
        ptcl%z = r1*grd_zarr(k+1) + &
             (1d0-r1)*grd_zarr(k)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
        ptcl%y = min(ptcl%y,grd_yarr(j+1))
        ptcl%y = max(ptcl%y,grd_yarr(j))
        ptcl%z = min(ptcl%z,grd_zarr(k+1))
        ptcl%z = max(ptcl%z,grd_zarr(k))
!-- sampling azimuthal angle of direction
        r1 = rnd_r(rnd_state)
        om0 = pc_pi2*r1

!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = ptcl%x
           y0 = ptcl%y
           z0 = ptcl%z
!-- 1+dir*v/c
           mu1 = sqrt(1d0-mu0**2)*cos(om0)
           mu2 = sqrt(1d0-mu0**2)*sin(om0)
           cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!-- mu
           ptcl%mu = (mu0+z0/pc_c)/cmffact
           if(ptcl%mu>1d0) then
              ptcl%mu = 1d0
           elseif(ptcl%mu<-1d0) then
              ptcl%mu = -1d0
           endif
!-- om
           ptcl%om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
           if(ptcl%om<0d0) ptcl%om = &
                ptcl%om+pc_pi2
        else
           ptcl%mu = mu0
           ptcl%om = om0
        endif

!-- 1D
     case(11)
!-- calculating position
        r1 = rnd_r(rnd_state)
        prt_tlyrand = prt_tlyrand+1
        ptcl%x = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)
!-- must be inside cell
        ptcl%x = min(ptcl%x,grd_xarr(i+1))
        ptcl%x = max(ptcl%x,grd_xarr(i))
!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = ptcl%x
!-- 1+dir*v/c
           cmffact = 1d0+x0*mu0/pc_c
!-- mu
           ptcl%mu = (mu0+x0/pc_c)/cmffact
        else
           ptcl%mu = mu0
        endif
     endselect

!-- if velocity-dependent, transforming energy, wavelength
     if(grd_isvelocity) then
        ptcl%e = ep0*cmffact
        ptcl%e0 = ep0*cmffact
        ptcl%wl = wl0/cmffact
     else
        ptcl%e = ep0
        ptcl%e0 = ep0
        ptcl%wl = wl0
     endif

!-- save particle result
!-----------------------
     prt_particles(ipart) = ptcl

  enddo !ipart


end subroutine initial_particles
