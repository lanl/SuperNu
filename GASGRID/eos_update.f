      subroutine eos_update(do_output)
c     --------------------------------
      use gasgridmod
      use ionsmod
      use timestepmod, only:tim_itim,tim_itc
      use timingmod
      implicit none
      logical,intent(in) :: do_output
************************************************************************
* Solve the eos for given temperatures.
************************************************************************
      integer :: i,niter,iion,nion
      integer :: icg,iz,ii
      real :: t0,t1
      real*8 :: ndens
      real*8 :: pdens(ion_nion,gg_ncg)
c
c-- loop over all ggrid cells
      call time(t0)
      do icg=1,gg_ncg
       ndens = ggrid2(icg)%natom/ggrid2(icg)%vol !atom number density
       call ion_solve_eos(ggrid2(icg)%natom1fr(1:),
     &   ggrid2(icg)%temp,ndens,ggrid2(icg)%nelec,niter)
c
c-- debug output
!      write(6,*) icg,niter !DEBUG
c-- mark gcell to need new opacities
       if(niter>1) ggrid2(icg)%opdirty = .true.
c
c-- store occupation numbers of each ion's ground states
       do iz=1,gg_nelem
        do ii=1,ion_grndlev(iz,icg)%ni
         ion_grndlev(iz,icg)%g(ii) = ion_el(iz)%i(ii)%glev(1)
         ion_grndlev(iz,icg)%oc(ii) =
     &     ggrid2(icg)%natom*ggrid2(icg)%natom1fr(iz)*
     &     ion_el(iz)%i(ii)%glev(1) * ion_el(iz)%i(ii)%n /
     &     (ion_el(iz)%i(ii)%q * ggrid2(icg)%volcrp) !number density, not number
         !write(*,*) iz,ii,ion_grndlev(iz,icg)%oc(ii) !ion_el(iz)%i(ii)%nlev,ion_el(iz)%i(ii)%glev(1) !DEBUG
        enddo !ii
       enddo !iz
c
c-- store partial densities
       if(do_output) then
        iion = 0!{{{
        do iz=1,gg_nelem
         do ii=1,ion_grndlev(iz,icg)%ni
          iion = iion + 1
          pdens(iion,icg) = ion_el(iz)%i(ii)%n
         enddo
        enddo!}}}
       endif
c
      enddo !icg
      call time(t1)
      call timereg(t_eos, t1-t0)
c
c
c-- print partial densities
      if(do_output) then
       write(8,*)!{{{
       write(8,*) 'partial densities:',tim_itim,tim_itc
c-- electron density
       write(8,'(2a12)') 'nelec','elec_dens' ![atom^-1],[cm^-3]
       do icg=1,gg_ncg
        write(8,'(1p,2e12.4)') ggrid2(icg)%nelec,
     &    ggrid2(icg)%nelec*ggrid2(icg)%natom/ggrid2(icg)%vol
       enddo
c-- partial densities
       nion = 0
       write(8,*)
       do iz=1,gg_nelem
        write(8,'(40i12)') (iz*100+i,i=1,ion_grndlev(iz,1)%ni)
        do icg=1,gg_ncg
         write(8,'(1p,40e12.4)') (pdens(nion+i,icg),
     &     i=1,ion_grndlev(iz,1)%ni)
        enddo !icg
        write(8,*)
        nion = nion + ion_grndlev(iz,1)%ni
       enddo !iz!}}}
      endif !do_output
c
      end subroutine eos_update
