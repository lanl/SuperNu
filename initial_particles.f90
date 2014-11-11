subroutine initial_particles

  use gridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine instantiates the initial particles before
  !the first time step.
!##################################################
!
  logical :: lhelp
  integer :: ig, i,j,k, iig, ipart, ihelp, jhelp
  integer, dimension(grd_nx,grd_ny,grd_nz) :: ijkused
  real*8 :: wl0, mu0, om0, ep0, x0, y0
  real*8 :: denom2
  real*8 :: r1
  real*8 :: wl1,wl2,wl3,wl4
  real*8 :: cmffact, azitrfm

!-- helper quantities
  wl1=1d0/grd_wl(grd_ng+1)
  wl2=1d0/grd_wl(1)

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
     prt_particles(ipart)%zsrc = i

!-- setting particle index to not vacant
     prt_isvacant(ipart) = .false.
!
!-- calculating particle time
     prt_particles(ipart)%tsrc = tsp_t
!
!-- setting type
     prt_particles(ipart)%rtsrc = 1

!-- calculating wavelength
     denom2 = 0d0
     r1=rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, grd_ng
        wl3 = 1d0/grd_wl(ig+1)
        wl4 = 1d0/grd_wl(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2 + (wl4-wl3)/(wl2-wl1)) exit
        denom2 = denom2 + (wl4-wl3)/(wl2-wl1)
     enddo
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/grd_wl(iig)+r1/grd_wl(iig+1))

!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1

!-- calculating particle energy
     ep0 = grd_evolinit(i,j,k)/real(grd_nvolinit(i,j,k))

!
!-- selecting geometry
     select case(in_igeom)

!-- 1D
     case(1)
!-- calculating position
        r1 = rand()
        prt_tlyrand = prt_tlyrand+1
        prt_particles(ipart)%rsrc = (r1*grd_xarr(i+1)**3 + &
             (1d0-r1)*grd_xarr(i)**3)**(1d0/3d0)

!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = prt_particles(ipart)%rsrc
!-- 1+dir*v/c
           cmffact = 1d0+x0*mu0/pc_c
!-- mu
           prt_particles(ipart)%musrc = (mu0+x0/pc_c)/cmffact
        else
           prt_particles(ipart)%musrc = mu0
        endif

!-- 2D
     case(2)
!-- setting 2nd cell index
        prt_particles(ipart)%iy = j
!-- calculating position
        r1 = rand()
        prt_particles(ipart)%rsrc = sqrt(r1*grd_xarr(i+1)**2 + &
             (1d0-r1)*grd_xarr(i)**2)
        r1 = rand()
        prt_particles(ipart)%y = r1*grd_yarr(j+1) + &
             (1d0-r1)*grd_yarr(j)

!-- sampling azimuthal angle of direction
        r1 = rand()
        om0 = pc_pi2*r1

!-- if velocity-dependent, transforming direction
        if(grd_isvelocity) then
           x0 = prt_particles(ipart)%rsrc
           y0 = prt_particles(ipart)%y
!-- 1+dir*v/c
           cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
           azitrfm = atan2(sqrt(1d0-mu0**2)*sin(om0), &
                sqrt(1d0-mu0**2)*cos(om0)+x0/pc_c)
!-- mu
           prt_particles(ipart)%musrc = (mu0+y0/pc_c)/cmffact
           if(prt_particles(ipart)%musrc>1d0) then
              prt_particles(ipart)%musrc = 1d0
           elseif(prt_particles(ipart)%musrc<-1d0) then
              prt_particles(ipart)%musrc = -1d0
           endif
!-- om
           if(azitrfm >= 0d0) then
              prt_particles(ipart)%om = azitrfm
           else
              prt_particles(ipart)%om = azitrfm+pc_pi2
           endif
        else
           prt_particles(ipart)%musrc = mu0
           prt_particles(ipart)%om = om0
        endif

!-- 3D
     case(3)
        stop 'initial_particles: no 3D transport'
     endselect

!-- if velocity-dependent, transforming energy, wavelength
     if(grd_isvelocity) then
        prt_particles(ipart)%esrc = ep0*cmffact
        prt_particles(ipart)%ebirth = ep0*cmffact
        prt_particles(ipart)%wlsrc = wl0/cmffact
     else
        prt_particles(ipart)%esrc = ep0
        prt_particles(ipart)%ebirth = ep0
        prt_particles(ipart)%wlsrc = wl0
     endif

  enddo !ipart


end subroutine initial_particles
