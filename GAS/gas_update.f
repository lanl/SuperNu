      subroutine gas_update(impi,it)
c     ------------------------------
      use gridmod
      use physconstmod
      use miscmod
      use ionsmod
      use timestepmod
      use totalsmod
      use gasmod
      use inputparmod
      use timingmod
      implicit none
      integer,intent(in) :: impi,it
************************************************************************
* Update the part of the gas grid that depends on time and temperature.
* The non-changing part is computed in gas_setup.
* The work done here is:
* - nuclear decay (energy release and chemical composition update)
* - temperature and volume
* - LTE EOS: ionization balance and electron density
* - opacities
************************************************************************
      logical,save :: lfirst=.true.
      logical :: do_output,lexist
      integer :: i,j,l,istat
      real*8 :: help
      real*8 :: dtempfrac = 0.99d0
      real*8 :: natom1fr(gas_ncell,-2:-1) !todo: memory storage order?
      real*8 :: natom2fr(gas_ncell,-2:-1)
      real*8 :: nisource0(gas_ncell)
c-- previous values
      real*8,allocatable,save :: tempalt(:),capgreyalt(:)
!     real*8 :: hlparr(grd_nx),hlparrdd(gas_ncell)
c-- timing
      real*8 :: t0,t1
c
c-- begin
      t0 = t_time()
c
c-- nuclear decay
c================
c-- Get ni56 and co56 abundances on begin and end of the time step.!{{{
c-- The difference between these two has decayed.
c
c-- initial decay, prior to first time step
      if(tsp_it==1 .and. grd_isvelocity.and.in_srctype=='none') then
       call update_natomfr(0d0)!{{{
       forall(i=-2:-1) natom1fr(:,i) = gas_natom1fr(i,:)
c-- end of time step
       call update_natomfr(tsp_t)
       forall(i=-2:-1) natom2fr(:,i) = gas_natom1fr(i,:)
c-- energy deposition
       nisource0 =  !per average atom (mix of stable and unstable)
     &   (natom1fr(:,gas_ini56) - natom2fr(:,gas_ini56)) *
     &     (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(:,gas_ico56) - natom2fr(:,gas_ico56)) *
     &     pc_qhl_co56
c-- total, units=ergs
       nisource0 = nisource0 * gas_natom
       tot_eext0 = sum(nisource0)!}}}
      endif
c
c-- current time step
      if(grd_isvelocity.and.in_srctype=='none') then
c-- beginning of time step
       help = tsp_t
       call update_natomfr(help)
       forall(i=-2:-1) natom1fr(:,i) = gas_natom1fr(i,:)
c-- end of time step
       call update_natomfr(tsp_t + tsp_dt)
       forall(i=-2:-1) natom2fr(:,i) = gas_natom1fr(i,:)
c
c-- update the abundances for the center time
       !call update_natomfr(tsp_tcenter)
       call update_natomfr(tsp_t)
       call update_ye
c
c-- energy deposition
       gas_nisource =  !per average atom (mix of stable and unstable)
     &   (natom1fr(:,gas_ini56) - natom2fr(:,gas_ini56)) *
     &     (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(:,gas_ico56) - natom2fr(:,gas_ico56)) *
     &     pc_qhl_co56
c-- total, units=ergs
       gas_nisource = gas_nisource * gas_natom
c
c-- We assume that half of the source energy goes into doppler losses
c-- the other half into heating the radiation field.
c-- We keep adding nisource/2 until the radiation field converges.
c-- The doppler losses are then equal to nisource/2.  We give two full
c-- nisource shots to bring the radiation field up to speed and help
c-- converge more quickly.
       if(it<1 .and. it>in_ntres+1) gas_nisource = gas_nisource*.5
      endif
!}}}
c
c
c-- update volume
c========================================
      i = 0
      do l=gas_icell1,gas_icell1+gas_ncell-1
       i = i+1
       gas_vol(i) = grd_vol(l)
      enddo !l
c
c
c-- update density, heat capacity
c================================
      gas_rho = gas_mass/gas_vol
c-- Calculating power law heat capacity
      gas_bcoef = in_gas_cvcoef * gas_temp**in_gas_cvtpwr *
     &  gas_rho**in_gas_cvrpwr
c-- temperature
      gas_ur = pc_acoef*gas_temp**4
c
c-- sanity check temperatures
      if(any(gas_temp/=gas_temp)) stop 'gas_temp NaN'
      if(any(gas_temp<=0d0)) stop 'gas_temp<=0'
c
c
c-- totals
c=========
c-- add initial thermal input to dd_eext
      if(tsp_it==1) then
!-- total comoving material energy
       tot_emat = sum(gas_bcoef*gas_temp*gas_vol)
       tot_eext = tot_eext + tot_emat  !was initialized either in totalsmod or in totals_startup
      endif
c
c
c
c-- compute the starting tempurature derivative in the fleck factor
      if(lfirst .or. in_opacanaltype/='none') then
c-- temporarily change!{{{
       gas_temp = dtempfrac*gas_temp
       if(in_opacanaltype=='none') then
        if(.not.in_noeos) call eos_update(.false.)
       endif
c
       if(in_opacanaltype/='none') then
        call analytic_opacity
       else
        if(in_ngs==0) then
         call physical_opacity
        else
         call physical_opacity_subgrid
        endif
       endif
       call opacity_planckmean
c
c-- save
       if(.not.allocated(tempalt)) then
        allocate(tempalt(gas_ncell))
        allocate(capgreyalt(gas_ncell))
       endif
       tempalt = gas_temp
       capgreyalt = gas_capgrey/gas_rho !per gram
c
c-- change back
       gas_temp = gas_temp/dtempfrac
!}}}
      endif
c
c
c
c-- solve LTE EOS
c================
      do_output = (in_pdensdump=='each' .or.
     &  (in_pdensdump=='one' .and. tsp_it==1))
      if(.not.in_noeos) call eos_update(do_output)
c
c
c
c-- calculate opacities
c======================
c-- gamma opacity
      gas_capgam = in_opcapgam*gas_ye*gas_rho
c
c
c-- simple analytical group/grey opacities: Planck and Rosseland 
      if(in_opacanaltype/='none') then
       call analytic_opacity
      else
c-- calculate physical opacities
c-- test existence of input.opac file
       inquire(file='input.opac',exist=lexist)
       if(.not.lexist) then
c-- calculate opacities
        if(in_ngs==0) then
         call physical_opacity
        else
         call physical_opacity_subgrid
        endif
       else
c-- read in opacities
        open(4,file='input.opac',status='old',iostat=istat)!{{{
        if(istat/=0) stop 'read_opac: no file: input.opac'
c-- read header
        read(4,*,iostat=istat)
        if(istat/=0) stop 'read_opac: file empty: input.opac'
c-- read each cell individually
        do j=1,tsp_it
c-- skip delimiter
         read(4,*,iostat=istat)
         if(istat/=0) stop 'read_opac: delimiter error: input.opac'
c-- read data
         do i=1,gas_ncell
          read(4,*,iostat=istat) help,gas_sig(i),gas_cap(:,i)
          if(istat/=0) stop 'read_opac: body error: input.opac'
         enddo !i
        enddo !j
        close(4)
        write(6,*) 'read_opac: read successfully'
!}}}
       endif
      endif
      call opacity_planckmean
c
c
c-- write out opacities
c----------------------
      if(trim(in_opacdump)=='off') then !{{{
c-- donothing
      else
       open(4,file='output.opac',status='unknown',position='append')
      endif !off
c
c-- write opacity grid
      inquire(4,opened=do_output)
      if(do_output) then
c-- header
       if(tsp_it==1) write(4,'("#",3i8)') gas_ncell,tsp_nt
       write(4,'("#",3i8)') tsp_it
c-- body
       do i=1,gas_ncell
        write(4,'(1p,9999e12.4)') gas_temp(i),gas_sig(i),gas_cap(:,i)
       enddo
c-- close file
       close(4)
      endif !do_output !}}}
c
c
c-- Calculating Fleck factor, leakage opacities
      call fleck_factor(tempalt,capgreyalt)
c
c
c-- save previous values for gentile-fleck factor calculation in next iter
      tempalt = gas_temp
      capgreyalt = gas_capgrey/gas_rho
c
      lfirst = .false.
c
      t1 = t_time()
      call timereg(t_gasupd,t1-t0)
c
      end subroutine gas_update
c
c
c
      subroutine update_natomfr(t)
c     ----------------------------!{{{
      use physconstmod
      use gasmod
      use inputparmod
      implicit none
      real*8,intent(in) :: t
************************************************************************
* update natom fractions for nuclear decay
************************************************************************
      real*8 :: expni,expco,help
c
      expni = exp(-t/pc_thl_ni56)
      expco = exp(-t/pc_thl_co56)
c
c-- update Fe
      help = 1d0 + (pc_thl_co56*expco - pc_thl_ni56*expni)/
     &  (pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Fe < 0'
      gas_natom1fr(26,:) = gas_natom0fr(gas_ini56,:)*help+!initial Ni56
     &  gas_natom0fr(gas_ico56,:)*(1d0-expco) +          !initial Co56
     &  gas_natom0fr(0,:)                                !initial Fe (stable)
c
c-- update Co56 and Co
      help = pc_thl_co56*(expni - expco)/(pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Co < 0'
c-- Co56
      gas_natom1fr(gas_ico56,:) =
     &  gas_natom0fr(gas_ini56,:)*help +  !initial Ni56
     &  gas_natom0fr(gas_ico56,:)*expco   !initial Co56
c-- Co
      gas_natom1fr(27,:) = gas_natom1fr(gas_ico56,:) +  !unstable
     &  gas_natom0fr(1,:)                                      !initial Co (stable)
c
c-- update Ni56 and Ni
c-- Ni56
      gas_natom1fr(gas_ini56,:) =
     &  gas_natom0fr(gas_ini56,:)*expni  !initial Ni56
c-- Ni
      gas_natom1fr(28,:) = gas_natom1fr(gas_ini56,:) + !unstable
     &  gas_natom0fr(2,:)                              !initial Ni (stable)
c!}}}
      end subroutine update_natomfr



      subroutine update_ye
c     --------------------!{{{
      use gasmod
      use elemdatamod
      implicit none
************************************************************************
* update nuclear electron fractions
************************************************************************
      integer :: l
c
      gas_ye = 0d0
      do l=1,gas_nelem
c-- wrong memory order, but this is a small array
       gas_ye = gas_ye + gas_natom1fr(l,:)*l/elem_data(l)%m
      enddo
c!}}}
      end subroutine update_ye
