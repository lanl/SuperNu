      subroutine dealloc_all
c     ----------------------
      use mpimod, only:impi,impi0
      use ionsmod, only:ion_dealloc
      use bbxsmod
      use gasgridmod
      use particlemod
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
      deallocate(gas_wl,gas_cap)
*     if(impi==impi0) deallocate(gas_vals2,gas_temphist,gas_dwl)
      deallocate(gas_vals2,gas_temphist,gas_dwl)
c-- gasgridmod
      deallocate(gas_numcensus,gas_rarr,gas_drarr)
      deallocate(gas_edep,gas_tempb)
      deallocate(gas_sigmapg,gas_sigmargleft)
      deallocate(gas_sigmargright,gas_emitprobg,gas_sigmal,gas_sigmar)
      deallocate(gas_ppl,gas_ppr)
      deallocate(gas_eraddensg)
      deallocate(gas_siggrey,gas_fcoef)
      deallocate(gas_sig,gas_sigbl,gas_sigbr)
      deallocate(gas_capgam)
      deallocate(gas_exsource)
      deallocate(gas_curvcent)
c-- particlemod
      deallocate(prt_particles)

      end subroutine dealloc_all
