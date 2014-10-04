subroutine temperature_update

  use gasgridmod
  use timestepmod
  use physconstmod
  use manufacmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine updates the material state.  At the moment
  !it updates material temperature and the approximate amount
  !of gamma ray energy introduced in the time step.
!##################################################
  integer :: i,j,k, ig
  real*8 :: dtemp, dtemp2
  real*8,parameter :: tauni = 8.8d0*86400.0d0
  real*8,parameter :: tauco = 111.3d0*86400.0d0
  real*8 :: temphelp

  !calculating radiation energy density
  gas_vals2%eraddens = 0d0
  gas_eraddens=gas_eraddens/gas_vals2%vol
  gas_vals2%eraddens = gas_eraddens

!-- calculating temperature
  if(allocated(gas_temppreset)) then
!-- apply read-in temperature profile
   gas_temp = gas_temppreset(:,:,:,tsp_it)
  else
!-- calculate temp correction
   do k=1,gas_nz
   do j=1,gas_ny
   do i=1,gas_nx
     dtemp = gas_edep(i,j,k)/gas_vals2(i,j,k)%vol !new
     !write(6,*) gas_edep(i,j,k), gas_vals2(i,j,k)%vol
     dtemp = (dtemp - tsp_dt*gas_fcoef(i,j,k)*gas_siggrey(i,j,k)*pc_c*gas_vals2(i,j,k)%ur)/gas_vals2(i,j,k)%bcoef
     
     dtemp2 = (gas_fcoef(i,j,k)/gas_vals2(i,j,k)%bcoef)*tsp_dt* &
          gas_vals2(i,j,k)%matsrc

     gas_temp(i,j,k)=gas_temp(i,j,k)+dtemp+dtemp2

   enddo !i
   enddo !j
   enddo !k
  endif

  gas_vals2%ur = pc_acoef*gas_temp**4

!-- reset physical gas_siggreyold (used in fleck_factor, only approximate)
  if(gas_isvelocity.and.in_opacanaltype=='none') then
     gas_siggreyold=gas_siggrey*(tsp_t/(tsp_t+tsp_dt))**3
  endif
!
!-- summing comoving material energy
  gas_emat = 0d0
  do k=1,gas_nz
  do j=1,gas_ny
  do i=1,gas_nx
     gas_emat=gas_emat+ &
          gas_vals2(i,j,k)%bcoef*gas_temp(i,j,k)* &
          gas_vals2(i,j,k)%vol
  enddo !i
  enddo !j
  enddo !k

end subroutine temperature_update
