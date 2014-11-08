      subroutine opacity_planckmean
c     -----------------------------
      use physconstmod
      use mpimod
      use gasmod
      use miscmod
      implicit none
************************************************************************
* Convert opacity into grey opacities
************************************************************************
      integer :: i
!     integer :: ig
!     real*8 :: x1,x2
      real*8 :: help,hlparr(gas_ng+1)
c-- Planck opacity
      gas_siggrey = 0d0
      do i=1,gas_ncell
       help = pc_h*pc_c/(pc_kb*gas_temp(i))
       hlparr = help/gas_wl
       gas_siggrey(i) = 15d0/pc_pi**4 *
     &   sum(gas_cap(:,i)*specint(hlparr(2:),hlparr(:gas_ng),3))
!      do ig=1,gas_ng
!       x1 = pc_h*pc_c/(gas_wl(ig + 1)*pc_kb*gas_temp(i))
!       x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(i))
!       gas_siggrey(i) = gas_siggrey(i) + 15d0/pc_pi**4*
!    &    gas_cap(ig,i)*specint(x1,x2,3)
!      enddo
      enddo !i
c
      end subroutine opacity_planckmean
