      subroutine fluxtally(it)
c     ------------------------
      use timestepmod
      use fluxmod
      use particlemod
      use miscmod
      implicit none
      integer,intent(in) :: it
************************************************************************
* Go through all particles and tally the ones that have left the domain
* and end up in this flux time step.
************************************************************************
      integer :: ipart
      integer :: ig,imu,iom
      real*8 :: help
      type(packet),pointer :: ptcl

      flx_lumtime = 0d0
      flx_luminos = 0d0
      flx_lumdev = 0d0
      flx_lumnum = 0
c
!!c$omp do schedule(static,1) !round-robin
      do ipart=1,prt_npartmax

c-- check vacancy
       if(prt_isvacant(ipart)) cycle
c
c-- active particle
       ptcl => prt_particles(ipart)
c
c-- check flux status
       if(ptcl%x/=huge(help)) cycle
c
c-- check flux time
       if(ptcl%t>tsp_tarr(it+1)) cycle
c
c-- retrieving lab frame flux group, polar, azimuthal bin
       ig = binsrch(ptcl%wl,flx_wl,flx_ng+1,.false.)
       imu = binsrch(ptcl%mu,flx_mu,flx_nmu+1,.false.)
       iom = binsrch(ptcl%om,flx_om,flx_nom+1,.false.)
c
c-- tallying outbound luminosity
       flx_lumtime(ig,imu,iom) = flx_lumtime(ig,imu,iom)+ptcl%t
       flx_luminos(ig,imu,iom) = flx_luminos(ig,imu,iom)+ptcl%e
       flx_lumdev(ig,imu,iom) = flx_lumdev(ig,imu,iom)+ptcl%e**2
       flx_lumnum(ig,imu,iom) = flx_lumnum(ig,imu,iom)+1

c-- mark particle slot occupied or vacant
       prt_isvacant(ipart) = .true.

      enddo !ipart
!!c$omp end do nowait

!-- convert to flux per second
      help = 1d0/tsp_dt
      flx_luminos = flx_luminos*help
      flx_lumdev = flx_lumdev*help**2
c
      end subroutine fluxtally
