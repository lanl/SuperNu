      subroutine dealloc_all
c     ----------------------
      use mpimod, only:impi,impi0
      use ionsmod, only:ion_dealloc
      use bbxsmod
      use gridmod
      use gasmod
      use particlemod
      use inputparmod
      use fluxmod
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
c-- gasmod
      deallocate(grd_numcensus)
      deallocate(grd_xarr,grd_yarr,grd_zarr)
      deallocate(grd_edep,grd_temp,grd_methodswap)
      deallocate(grd_cap)
      deallocate(grd_opacleak)
      deallocate(grd_eraddens)
      if(impi==impi0) deallocate(grd_capgrey)
      deallocate(grd_fcoef)
      deallocate(grd_sig)
      deallocate(grd_emit,grd_emitex,grd_nvol)
      deallocate(grd_evolinit,grd_nvolinit)
c-- particlemod
      if(impi==impi0.and..not.in_norestart) deallocate(prt_tlyrandarr)
      deallocate(prt_particles)
c-- fluxmod
      deallocate(flx_wl,flx_mu,flx_om)
      deallocate(flx_luminos,flx_lumdev,flx_lumnum)

      end subroutine dealloc_all
