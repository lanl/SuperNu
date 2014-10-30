subroutine emission_probability

  use gasgridmod
  use physconstmod
  use miscmod, only:specint
  implicit none

!-----------------------
  !multigroup volume emission probabilities
!-----------------------

  integer :: i,j,k, ig
  real*8 :: x1, x2

!-- init
  gas_emitprob = 0d0

!-- Calculating grouped volume emission probabilities:
  if(gas_opacanaltype=='pick') then
     do k=1,gas_nz
     do j=1,gas_ny
     do i=1,gas_nx
        gas_emitprob(1,i,j,k) = gas_ppick(1)*gas_cap(1,i,j,k)/gas_siggrey(i,j,k)
        gas_emitprob(2,i,j,k) = gas_ppick(2)*gas_cap(2,i,j,k)/gas_siggrey(i,j,k)
!       gas_emitprob(3:gas_ng,i,j,k) = 0d0  !-- not necessary
     enddo !i
     enddo !j
     enddo !k
  else
     if(gas_ng==1) then
        gas_emitprob = 1d0
     else
        do k=1,gas_nz
        do j=1,gas_ny
        do i=1,gas_nx
           do ig=1,gas_ng
              x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_temp(i,j,k))
              x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(i,j,k))
              if(gas_siggrey(i,j,k)<=0d0) then
!                gas_emitprob(ig,i,j,k) = 0d0  !-- not necessary
              else
                 gas_emitprob(ig,i,j,k) = 15d0*specint(x1,x2,3)*gas_cap(ig,i,j,k)/ &
                      (gas_siggrey(i,j,k)*pc_pi**4)
              endif
           enddo !ig
        enddo !i
        enddo !j
        enddo !k
     endif
  endif

end subroutine emission_probability
