      subroutine dealloc_all
c     ----------------------
      use mpimod, only:impi,impi0
      use ionsmod, only:ion_dealloc
      use bbxsmod
      use gasgridmod
      use particlemod
      IMPLICIT NONE
************************************************************************
* deallocate all that was used till the end of the program. Any
* allocatable arrays remaining allocated after this had to be dealt
* with earlier.  This helps to catch memory leaks! (drr)
************************************************************************
c-- ionsmod
      if(impi==impi0) call ion_dealloc
c-- bbxsmod
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f
!c-- ggridmod
!      deallocate(ggrid,ggrid_wl,ggrid_cap,ggrid_icapbb)
!      if(impi==impi0) deallocate(ggrid2,ggrid2_temp,ggrid2_dwl)
c-- gasgridmod
      deallocate(gas_numcensus,gas_rarr,gas_drarr)
      deallocate(gas_edep,gas_tempb,gas_vals2)
      deallocate(gas_sigmapg,gas_sigmargleft)
      deallocate(gas_sigmargright,gas_emitprobg,gas_sigmal,gas_sigmar)
      deallocate(gas_ppl,gas_ppr,gas_ppick)
      deallocate(gas_sigmap,gas_fcoef)
c-- particlemod
      deallocate(prt_particles)
c
      end subroutine dealloc_all
