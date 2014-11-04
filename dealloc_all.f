      subroutine dealloc_all
c     ----------------------
      use mpimod, only:impi,impi0
      use ionsmod, only:ion_dealloc
      use bbxsmod
      use gasgridmod
      use particlemod
      use inputstrmod
      use inputparmod
      implicit none
************************************************************************
* deallocate all that was used till the end of the program. Any
* allocatable arrays remaining allocated after this had to be dealt
* with earlier.  This helps to catch memory leaks! (drr)
************************************************************************
c-- ionsmod
      !write(*,*) 'here 1.5 ...'
      if(impi==impi0) call ion_dealloc
c-- bbxsmod
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f
!c-- gasgridmod
      deallocate(gas_wl)
c-- gasgridmod
      deallocate(gas_numcensus)
      deallocate(gas_xarr,gas_yarr,gas_zarr)
      deallocate(gas_edep,gas_temp,gas_methodswap)
      deallocate(gas_cap)
      deallocate(gas_emitprob,gas_opacleak)
      deallocate(gas_eraddens)
      if(impi==impi0) deallocate(gas_siggrey)
      deallocate(gas_fcoef)
      deallocate(gas_sig)
      if(impi==impi0) deallocate(gas_capgam)
      deallocate(gas_emit,gas_emitex,gas_nvol,gas_nvolex)
      deallocate(gas_evolinit,gas_nvolinit)
c-- particlemod
      if(impi==impi0.and..not.in_norestart) deallocate(prt_tlyrandarr)
      deallocate(prt_particles)
c-- inputstrmod
      if(impi==impi0) deallocate(str_xleft)
      if(impi==impi0) deallocate(str_mass)
      if(allocated(str_massfr)) deallocate(str_massfr)
      if(allocated(str_iabund)) deallocate(str_iabund)
c-- fluxmod
      deallocate(flx_wl,flx_mu,flx_om)
      deallocate(flx_luminos,flx_lumdev,flx_lumnum)

      end subroutine dealloc_all
