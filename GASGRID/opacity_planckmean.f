      subroutine opacity_planckmean
c     -----------------------------
      use physconstmod
      use mpimod
      use gasgridmod
      use miscmod
      implicit none
************************************************************************
* Convert opacity into grey opacities
************************************************************************
      integer :: i,j,k,l
!     integer :: ig
!     real*8 :: x1,x2
      real*8 :: help,hlparr(gas_ng+1)
c-- Planck opacity
      gas_siggrey = 0d0
      do k=1,gas_nz
      do j=1,gas_ny
      do i=1,gas_nx
       help = pc_h*pc_c/(pc_kb*dd_temp(i,j,k))
       hlparr = help/gas_wl
       gas_siggrey(i,j,k) = 15d0/pc_pi**4 *
     &   sum(gas_cap(:,i,j,k)*specint(hlparr(2:),hlparr(:gas_ng),3))
!      do ig=1,gas_ng
!       x1 = pc_h*pc_c/(gas_wl(ig + 1)*pc_kb*dd_temp(i,j,k))
!       x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*dd_temp(i,j,k))
!       gas_siggrey(i,j,k) = gas_siggrey(i,j,k) + 15d0/pc_pi**4*
!    &    gas_cap(ig,i,j,k)*specint(x1,x2,3)
!      enddo
      enddo !i
      enddo !j
      enddo !k
c
      end subroutine opacity_planckmean
