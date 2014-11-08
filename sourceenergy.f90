subroutine sourceenergy(nmpi)

  use gasgridmod
  use totalsmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none
  integer,intent(in) :: nmpi

!##################################################
!This subroutine computes the distribution of source particles each
!time step.  A fraction of the source particle number prt_ns is given
!to each cell based on the amount of energy emitted by the cell.
!##################################################

  integer :: i
  integer :: ihelp
! tot_esurf for any new prt_particles from a surface source
! prt_nsurf = number of surface prt_particles
! prt_nnew = total number of new prt_particles~=prt_ns

  tot_esurf = 0d0
  dd_emit = 0d0
  dd_emitex = 0d0
  
  ! Calculating dd_emitex from analytic distribution
  call analytic_source

  ! Calculating fictitious emission energy per cell: loop
  do i = 1, dd_ncell
     
     dd_emit(i) =  tsp_dt*dd_fcoef(i)*dd_siggrey(i)*pc_c* &
          dd_ur(i)*dd_vol(i)
!
     dd_emit(i) = dd_emit(i) + &
          tsp_dt*dd_vol(i)*(1d0-dd_fcoef(i))*&
          dd_matsrc(i)
!-- accounting for total material source energy:
!-- dd_fcoef of matsrc goes directly in temperature equation
!-- and remaining 1-dd_fcoef is thermal radiation.
     if(tsp_it==1) then
        tot_eext = tot_eext+tsp_dt*dd_vol(i)*&
             dd_matsrc(i)
     else
!-- rtw: tot_eext is only broadcast for tsp_it==1
        tot_eext = tot_eext+tsp_dt*dd_vol(i)*&
             dd_matsrc(i)*dble(nmpi)
     endif
!
!-- accounting for thermal decay radiation source energy
     if(.not.in_novolsrc .and. in_srctype=='none') then
        dd_emit(i) = dd_emit(i) + dd_nisource(i)
        if(tsp_it==1) then
           tot_eext = tot_eext + dd_nisource(i)
        else
!-- rtw: tot_eext is only broadcast for tsp_it==1
           tot_eext = tot_eext + dble(nmpi)*dd_nisource(i)
        endif
     endif
  enddo

end subroutine sourceenergy
