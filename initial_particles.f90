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
    integer :: ig, iir, iig, ipart
    integer, dimension(gas_nx) :: iirused
    real*8 :: wl0, mu0, Ep0, r0
    real*8 :: denom2
    real*8 :: r1,r3
    real*8 :: x1,x2,x3,x4
    logical :: isnotvacnt

!-- helper quantities
    x1=1d0/gas_wl(gas_ng+1)
    x2=1d0/gas_wl(1)

    !insantiating initial particles
    iir = 1
    iirused = 0
    do ipart = 1, prt_ninitnew
       isnotvacnt=.false.
       do while(.not.isnotvacnt)
          if(iirused(iir)<gas_nvolinit(iir,1,1)) then
             iirused(iir)=iirused(iir)+1
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
             !calculating wavelegth unformly
             r1 = rand()
             prt_tlyrand = prt_tlyrand+1
             wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))
             !calculating radial position
             r3 = rand()
             prt_tlyrand = prt_tlyrand+1
             prt_particles(ipart)%rsrc = (r3*gas_xarr(iir+1)**3 + &
                  (1.0-r3)*gas_xarr(iir)**3)**(1.0/3.0)
             r0 = prt_particles(ipart)%rsrc
             !calculating direction cosine (comoving)
             r1 = rand()
             prt_tlyrand = prt_tlyrand+1
             mu0 = 1d0-2d0*r1
             if(abs(mu0)<0.0000001d0) then
                   mu0=0.0000001d0
             endif
             !calculating particle time
             prt_particles(ipart)%tsrc = tsp_t
             !calculating particle energy
             Ep0 = gas_evolinit(iir,1,1)/real(gas_nvolinit(iir,1,1))
!             gas_eext=gas_eext+Ep0
             if(gas_isvelocity) then
                prt_particles(ipart)%Esrc = Ep0*(1.0+r0*mu0/pc_c)
                prt_particles(ipart)%Ebirth = Ep0*(1.0+r0*mu0/pc_c)
!-- velocity effects accounting
!                 gas_evelo=gas_evelo-Ep0*r0*mu0/pc_c
!
                prt_particles(ipart)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
                !
                prt_particles(ipart)%musrc = (mu0+r0/pc_c)/&
                     (1.0+r0*mu0/pc_c)
             else
                prt_particles(ipart)%Esrc = Ep0
                prt_particles(ipart)%Ebirth = Ep0

                prt_particles(ipart)%wlsrc = wl0
                !
                prt_particles(ipart)%musrc = mu0
             endif
             prt_particles(ipart)%rtsrc = 1

             !Setting iir = zone of particle
             prt_particles(ipart)%zsrc = iir
             !Setting particle index to not vacant
             prt_particles(ipart)%isvacant = .false.

             isnotvacnt = .true.
          else
             iir=iir+1
          endif
       enddo
    enddo


end subroutine initial_particles
