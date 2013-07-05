      subroutine dealloc_all
c     ----------------------
      use mpimod, only:impi,impi0
      use ionsmod, only:ion_dealloc
      use bbxsmod
      use gasgridmod
      use particlemod
      use inputstrmod
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
      if(impi==impi0) deallocate(gas_vals2,gas_temphist,gas_dwl)
      if(impi==impi0) deallocate(gas_caprosl,gas_caprosr)
      !deallocate(gas_vals2,gas_temphist,gas_dwl)
c-- gasgridmod
      deallocate(gas_numcensus,gas_rarr,gas_drarr)
      deallocate(gas_edep,gas_temp,gas_tempb)
      if(impi==impi0) deallocate(gas_rhob)
      deallocate(gas_cap)
      deallocate(gas_emitprob,gas_opacleakl,gas_opacleakr)
      deallocate(gas_ppl,gas_ppr)
      deallocate(gas_eraddens)
      if(impi==impi0) deallocate(gas_siggrey)
      deallocate(gas_fcoef)
      deallocate(gas_sig)
      if(impi==impi0) deallocate(gas_sigbl,gas_sigbr)
      if(impi==impi0) deallocate(gas_exsource,gas_capgam)
      deallocate(gas_curvcent)
      deallocate(gas_emit,gas_emitex,gas_nvol,gas_nvolex)
c-- particlemod
      deallocate(prt_particles)
c-- inputstrmod
      deallocate(str_velright)
      deallocate(str_velleft)
      deallocate(str_mass)
      if(allocated(str_massfr)) deallocate(str_massfr)
      if(allocated(str_abundlabl)) deallocate(str_abundlabl)
      if(allocated(str_iabund)) deallocate(str_iabund)

      end subroutine dealloc_all
