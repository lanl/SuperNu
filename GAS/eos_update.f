      subroutine eos_update(do_output)
c     --------------------------------
      use gasmod
      use ionsmod
      use timestepmod, only:tsp_it
      use timingmod
      implicit none
      logical,intent(in) :: do_output
************************************************************************
* Solve the eos for given temperatures.
************************************************************************
      integer :: i,ll,niter,iion,nion,istat
      integer :: iz,ii
      real*8 :: t0,t1
      real*8 :: ndens
      real*8 :: pdens(ion_nion,gas_ncell)
c
c-- loop over all gas_vals cells
      t0 = t_time()
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       ndens = gas_natom(i)/gas_vol(i) !atom number density
       call ions_solve_eos(gas_natom1fr(1,i),
     &   gas_temp(i),ndens,gas_nelec(i),niter)
c
c-- store occupation numbers of each ion's ground states
       do iz=1,gas_nelem
        do ii=1,ion_grndlev(iz,i)%ni
         ion_grndlev(iz,i)%ginv(ii) = 1d0/ion_el(iz)%i(ii)%glev(1)
         ion_grndlev(iz,i)%oc(ii) =
     &     gas_natom(i)*gas_natom1fr(iz,i)*
     &     ion_el(iz)%i(ii)%glev(1) * ion_el(iz)%i(ii)%n /
     &     (ion_el(iz)%i(ii)%q * gas_vol(i)) !number density, not number
         !write(6,*) iz,ii,ion_grndlev(iz,i)%oc(ii) !ion_el(iz)%i(ii)%nlev,ion_el(iz)%i(ii)%glev(1) !DEBUG
        enddo !ii
       enddo !iz
c
c-- store partial densities
       if(do_output) then
        iion = 0!{{{
        do iz=1,gas_nelem
         do ii=1,ion_grndlev(iz,i)%ni
          iion = iion + 1
          pdens(iion,i) = ion_el(iz)%i(ii)%n
         enddo
        enddo!}}}
       endif
c
      enddo !ix
      t1 = t_time()
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
       write(4,*) gas_ncell,gas_nelem
       write(4,*) '# nion'
       write(4,'(40i12)') (ion_el(iz)%ni,iz=1,gas_nelem)
       write(4,*) '# ions'
       do iz=1,gas_nelem
        write(4,'(40i12)') (iz*100+ll,ll=1,ion_grndlev(iz,1)%ni)
       enddo
c
c-- electron density
       write(4,'(2a12)') '# nelec','elec_dens' ![atom^-1],[cm^-3]
       do i=1,gas_ncell
        write(4,'(1p,2e12.4)') gas_nelec(i),
     &    gas_nelec(i)*gas_natom(i)/gas_vol(i)
       enddo !i
c
c-- partial densities
       nion = 0
       do iz=1,gas_nelem
        write(4,'("#",40i12)') (iz*100+ll,ll=1,ion_grndlev(iz,1)%ni)
        do i=1,gas_ncell
         write(4,'(1p,40e12.4)') (pdens(nion+ll,i)*
     &     gas_natom1fr(iz,i),ll=1,ion_grndlev(iz,1)%ni)
        enddo !i
        nion = nion + ion_grndlev(iz,1)%ni
       enddo !iz
c
       close(4)
      endif !do_output
c
      end subroutine eos_update
