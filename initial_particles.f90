subroutine initial_particles

  use gasgridmod
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
  integer :: ig, i,j,k iig, ipart
  integer, dimension(gas_nx,gas_ny,gas_nz) :: ijkused
  real*8 :: wl0, mu0, ep0, r0
  real*8 :: denom2
  real*8 :: r1,r3
  real*8 :: x1,x2,x3,x4

!-- helper quantities
  x1=1d0/gas_wl(gas_ng+1)
  x2=1d0/gas_wl(1)

!-- instantiating initial particles
  i = 1
  j = 1
  k = 1
  ijkused = 0
  do ipart=1,prt_ninitnew
     do k=k,gas_nz
     do j=j,gas_ny
     do i=i,gas_nx
        if(ijkused(i,j,k)<gas_nvolinit(i,j,k)) exit
     enddo
     enddo
     enddo
!--
     ijkused(i,j,k) = ijkused(i,j,k)+1
     denom2 = 0d0
     r1=rand()
     prt_tlyrand = prt_tlyrand+1
     do ig = 1, gas_ng
        x3 = 1d0/gas_wl(ig+1)
        x4 = 1d0/gas_wl(ig)
        iig = ig
        if(r1>=denom2.and.r1<denom2 + (x4-x3)/(x2-x1)) exit
        denom2 = denom2 + (x4-x3)/(x2-x1)
     enddo
!-- calculating wavelegth unformly
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))
!-- calculating radial position
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     
     prt_particles(ipart)%rsrc = (r3*gas_xarr(i+1)**3 + &
          (1.0-r3)*gas_xarr(i)**3)**(1.0/3.0)
     r0 = prt_particles(ipart)%rsrc
!-- calculating direction cosine (comoving)
     r1 = rand()
     prt_tlyrand = prt_tlyrand+1
     mu0 = 1d0-2d0*r1
     mu0 = max(mu0,1d-7)
!-- calculating particle time
     prt_particles(ipart)%tsrc = tsp_t
!-- calculating particle energy
     ep0 = gas_evolinit(i,1,1)/real(gas_nvolinit(i,1,1))
!             gas_eext=gas_eext+ep0
     if(gas_isvelocity) then
        prt_particles(ipart)%Esrc = ep0*(1.0+r0*mu0/pc_c)
        prt_particles(ipart)%Ebirth = ep0*(1.0+r0*mu0/pc_c)
!-- velocity effects accounting
!                 gas_evelo=gas_evelo-ep0*r0*mu0/pc_c
!
        prt_particles(ipart)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
        !
        prt_particles(ipart)%musrc = (mu0+r0/pc_c)/&
             (1.0+r0*mu0/pc_c)
     else
        prt_particles(ipart)%Esrc = ep0
        prt_particles(ipart)%Ebirth = ep0

        prt_particles(ipart)%wlsrc = wl0
        !
        prt_particles(ipart)%musrc = mu0
     endif
     prt_particles(ipart)%rtsrc = 1

!-- Setting i = zone of particle
     prt_particles(ipart)%zsrc = i
!-- Setting particle index to not vacant
     prt_isvacant(ipart) = .false.

  enddo !ipart


end subroutine initial_particles
