subroutine tabular_source(it,nt,tcenter,srctype)
  use tbsrcmod
  use gasmod
  use miscmod, only:binsrch
  implicit none
  integer,intent(in) :: it, nt
  real*8,intent(in) :: tcenter
  character(4),intent(in) :: srctype
!--------------------------------------------------
!- Calculate heating rate from source tables.
!- This calculation uses Barnes & Kasen (2016) for
!- the thermalization fractions.
!--------------------------------------------------
  integer,parameter :: icols(4) = [1,2,4,5]
  real*8,allocatable,save :: tdyn(:),twnd(:)
  integer :: itdyn, itwnd, icol, i
  real*8 :: helparr(gas_ncell)
  real*8 :: therm_fracs(gas_ncell,tbs_ncol-3)

!-- sanity checks
  if(tbs_ntdyn<=0 .and. tbs_ntwnd<=0) stop &
       'tabular_source: ntdyn<=0 && ntwnd<=0'
  if(srctype/='tabl') stop &
       'tabular_source: invalid srctype'

!-- reset material source
  gas_matsrc = 0d0

!-- initialize helper and thermalization fraction array
  helparr = 0d0
  therm_fracs = 0d0

!-- allocate and set time arrays once
  if(.not.allocated(tdyn) .and. tbs_ntdyn>0) then
     allocate(tdyn(tbs_ntdyn))
     tdyn = tbs_dynhr(1,:)
  endif
  if(.not.allocated(twnd) .and. tbs_ntwnd>0) then
     allocate(twnd(tbs_ntwnd))
     twnd = tbs_wndhr(1,:)
  endif

!-- set the helper array
  if(any(tbs_dynopt==1)) then
     where(gas_rho>0d0) helparr = 2d0/(tcenter*gas_rho)
  endif

!-- calculate the dynamical ejecta thermalization fractions
  do icol=1,tbs_ncol-3
     if(tbs_dynopt(icol)==1) then
        where(helparr < huge(helparr) .and. helparr > 0d0)
           therm_fracs(:,icol) = log(1d0+(helparr*tbs_dynpars(icol))) / &
                (helparr*tbs_dynpars(icol))
        endwhere
     else
        therm_fracs(:,icol) = tbs_dynpars(icol)
     endif
  enddo

!-- sanity check
  if(any(therm_fracs>1d0)) stop 'tabular_source: therm_fracs>1'

!-- calculate dynamical ejecta rate
  if(tbs_ntdyn>0) then
!-- find rates at time closest to tcenter
     itdyn = binsrch(tcenter,tdyn,tbs_ntdyn,.false.)
!-- exclude column for energy from gamma (hard-coded as 6 for now)
     helparr = 0d0
     do i=1,4
        icol = icols(i)
        helparr = helparr+therm_fracs(:,icol)*tbs_dynhr(icol+3,itdyn)
     enddo
!-- now multiply by dynamical ejecta fraction
     gas_matsrc = helparr * gas_dynfr
  endif

!-- calculate the wind thermalization fractions
  do icol=1,tbs_ncol-3
     if(tbs_wndopt(icol)==1) then
        where(helparr < huge(helparr) .and. helparr > 0d0)
           therm_fracs(:,icol) = log(1d0+(helparr*tbs_wndpars(icol))) / &
                (helparr*tbs_wndpars(icol))
        endwhere
     else
        therm_fracs(:,icol) = tbs_wndpars(icol)
     endif
  enddo

!-- sanity check
  if(any(therm_fracs>1d0)) stop 'tabular_source: therm_fracs>1'

!-- calculate wind rate
  if(tbs_ntwnd>0) then
!-- find rates at time closest to tcenter
     itwnd = binsrch(tcenter,twnd,tbs_ntwnd,.false.)
!-- exclude column for energy from gamma (hard-coded as 6 for now)
     helparr = 0d0
     do i=1,4
        icol = icols(i)
        helparr = helparr+therm_fracs(:,icol)*tbs_wndhr(icol+3,itwnd)
     enddo
!-- add wind contribution to material source
     gas_matsrc = gas_matsrc + helparr * (1d0-gas_dynfr)
  endif

!-- convert source to per-volume units
  gas_matsrc = gas_matsrc * gas_rho

!-- remove sources from zero-mass cells (and catch nans)
  where(gas_mass<=0d0) gas_matsrc = 0d0

!-- create external gamma source for grey mc transport
  gas_decaygamma = 0d0
  if(tbs_ntdyn>0) &
       gas_decaygamma = gas_dynfr*tbs_dynhr(6,itdyn)
  if(tbs_ntwnd>0) gas_decaygamma = gas_decaygamma + &
       (1d0-gas_dynfr)*tbs_wndhr(6,itwnd)
!-- convert source to per-volume units
  gas_decaygamma = gas_decaygamma * gas_rho
!-- remove sources from zero-mass cells (and catch nans)
  where(gas_mass<=0d0) gas_decaygamma = 0d0

!-- clean up
  if(it==nt) then
     if(allocated(tdyn)) deallocate(tdyn)
     if(allocated(twnd)) deallocate(twnd)
  endif

end subroutine tabular_source
