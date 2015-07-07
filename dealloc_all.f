      subroutine dealloc_all
c     ----------------------
      use mpimod
      use ionsmod
      use bbxsmod
      use gridmod
      use groupmod
      use gasmod
      use particlemod
      use fluxmod
      implicit none
************************************************************************
* deallocate all that was used till the end of the program. Any
* allocatable arrays remaining allocated after this had to be dealt
* with earlier.  This helps to catch memory leaks! (drr)
************************************************************************
c-- ionsmod
      call ions_dealloc
      call gas_dealloc
      call grid_dealloc
      call flux_dealloc
      deallocate(prt_particles,prt_isvacant)
      call mpimod_dealloc
      deallocate(grp_wl,grp_wlinv)
      if(allocated(bb_xs)) deallocate(bb_xs) !only impi==impi0, but only if nobbopac==f

      end subroutine dealloc_all
