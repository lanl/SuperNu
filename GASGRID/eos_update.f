      subroutine eos_update(do_output)
c     --------------------------------
      use gasgridmod
      use ionsmod
      use timestepmod, only:tsp_it
      use timingmod
      implicit none
      logical,intent(in) :: do_output
************************************************************************
* Solve the eos for given temperatures.
************************************************************************
      integer :: i,niter,iion,nion,istat
      integer :: j,k,iz,ii
      real*8 :: t0,t1
      real*8 :: ndens
      real*8 :: pdens(ion_nion,gas_nx,gas_ny,gas_nz)
c
c-- loop over all gas_vals cells
      call time(t0)
      do k=1,gas_nz
      do j=1,gas_ny
      do i=1,gas_nx
       ndens = gas_vals2(i,j,k)%natom/gas_vals2(i,j,k)%vol !atom number density
       call ion_solve_eos(gas_vals2(i,j,k)%natom1fr(1:),
     &   gas_temp(i,j,k),ndens,gas_vals2(i,j,k)%nelec,niter)
c
c-- debug output
!      write(6,*) ir,niter !DEBUG
c-- mark gcell to need new opacities
       if(niter>1) gas_vals2(i,j,k)%opdirty = .true.
c
c-- store occupation numbers of each ion's ground states
       do iz=1,gas_nelem
        do ii=1,ion_grndlev(iz,i,j,k)%ni
         ion_grndlev(iz,i,j,k)%g(ii) = ion_el(iz)%i(ii)%glev(1)
         ion_grndlev(iz,i,j,k)%oc(ii) =
     &     gas_vals2(i,j,k)%natom*gas_vals2(i,j,k)%natom1fr(iz)*
     &     ion_el(iz)%i(ii)%glev(1) * ion_el(iz)%i(ii)%n /
     &     (ion_el(iz)%i(ii)%q * gas_vals2(i,j,k)%volcrp) !number density, not number
         !write(6,*) iz,ii,ion_grndlev(iz,i,j,k)%oc(ii) !ion_el(iz)%i(ii)%nlev,ion_el(iz)%i(ii)%glev(1) !DEBUG
        enddo !ii
       enddo !iz
c
c-- store partial densities
       if(do_output) then
        iion = 0!{{{
        do iz=1,gas_nelem
         do ii=1,ion_grndlev(iz,i,j,k)%ni
          iion = iion + 1
          pdens(iion,i,j,k) = ion_el(iz)%i(ii)%n
         enddo
        enddo!}}}
       endif
c
      enddo !ix
      enddo !iy
      enddo !iz
      call time(t1)
      call timereg(t_eos, t1-t0)
c
c
c-- print partial densities
      if(do_output) then
       open(4,file='output.pdens',status='unknown',position='append',
     &   iostat=istat)
       if(istat/=0) stop 'eos_update: error opening output file'
c
c-- write header
       write(4,*) '# partial densities',tsp_it
       write(4,*) '# nr, nelem'
       write(4,*) gas_nx,gas_ny,gas_nz,gas_nelem
       write(4,*) '# nion'
       write(4,'(40i12)') (ion_el(iz)%ni,iz=1,gas_nelem)
       write(4,*) '# ions'
       do iz=1,gas_nelem
        write(4,'(40i12)') (iz*100+i,i=1,ion_grndlev(iz,1)%ni)
       enddo
c
c-- electron density
       write(4,'(2a12)') '# nelec','elec_dens' ![atom^-1],[cm^-3]
       do k=1,gas_nz
       do j=1,gas_ny
       do i=1,gas_nx
        write(4,'(1p,2e12.4)') gas_vals2(i,j,k)%nelec,
     &    gas_vals2(i,j,k)%nelec*gas_vals2(i,j,k)%natom/gas_vals2(i,j,k)%vol
       enddo !i
       enddo !j
       enddo !k
c
c-- partial densities
       nion = 0
       do iz=1,gas_nelem
        write(4,'("#",40i12)') (iz*100+i,i=1,ion_grndlev(iz,1)%ni)
        do k=1,gas_nz
        do j=1,gas_ny
        do i=1,gas_nx
         write(4,'(1p,40e12.4)') (pdens(nion+i,i,j,k)*
     &     gas_vals2(i,j,k)%natom1fr(iz),
     &     i=1,ion_grndlev(iz,1)%ni)
        enddo !i
        enddo !j
        enddo !k
        nion = nion + ion_grndlev(iz,1)%ni
       enddo !iz
c
       close(4)
      endif !do_output
c
      end subroutine eos_update
