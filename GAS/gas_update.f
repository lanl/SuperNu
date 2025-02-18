* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine gas_update(it)
c     -------------------------
      use gridmod
      use physconstmod
      use nucdatamod
      use miscmod
      use ionsmod
      use timestepmod
      use totalsmod
      use gasmod
      use inputparmod
      use timingmod
      implicit none
      integer,intent(in) :: it
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
      real*8 :: natom1fr(-2*gas_nchain:gas_nelem,gas_ncell)
      real*8 :: natom2fr(-2*gas_nchain:gas_nelem,gas_ncell)
      real*8 :: decay0(gas_ncell)
      real*8 :: temp(gas_ncell)
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
       natom1fr = gas_natom1fr
c-- end of time step
       call update_natomfr(tsp_t)
       natom2fr = gas_natom1fr
c-- energy deposition
       decay0 =  !per average atom (mix of stable and unstable)
     &   (natom1fr(gas_ini56,:) - natom2fr(gas_ini56,:)) *
     &     (nuc_qhl_ni56 + nuc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(gas_ico56,:) - natom2fr(gas_ico56,:)) *
     &     nuc_qhl_co56 !+
c-- off for backwards compatibility
c    &   (natom2fr(26,:) - natom1fr(26,:))*nuc_q_poskin !beta decay
       tot_eext0 = sum(decay0*gas_natom) !total, units=ergs !}}}
      endif
c
c-- current time step
      if(grd_isvelocity.and.in_srctype=='none') then
c-- beginning of time step
       call update_natomfr(tsp_t)
       natom1fr = gas_natom1fr
c-- end of time step
       call update_natomfr(tsp_t + tsp_dt)
       natom2fr = gas_natom1fr
c
c-- update the abundances for the center time
       !call update_natomfr(tsp_tcenter)
       call update_natomfr(tsp_t)
c-- sanity check
       if(any(gas_natom1fr<0d0)) stop 'gas_update: natom1fr<0'
!c-- print change in electron fraction
!       do i=1,gas_ncell
!        write(0,*) i,gas_ye(i),gas_ye0(i),
!     &   2*(gas_ye(i) - gas_ye0(i))/(gas_ye(i) + gas_ye0(i))
!       enddo
c
c-- energy deposition
c-- gamma decay
       gas_decaygamma =  !per average atom (mix of stable and unstable)
     &   (natom1fr(gas_ini56,:) - natom2fr(gas_ini56,:)) *
     &     (nuc_qhl_ni56 + nuc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(gas_ico56,:) - natom2fr(gas_ico56,:)) *
     &     nuc_qhl_co56
       gas_decaygamma = gas_decaygamma * gas_natom !total, units=ergs
c-- beta decay (off for backwards compatibility)
c      gas_decaybeta = (natom2fr(26,:) - natom1fr(26,:))*nuc_q_poskin
c      gas_decaybeta = gas_decaybeta * gas_natom !total, units=ergs
c
c-- We assume that half of the source energy goes into doppler losses
c-- the other half into heating the radiation field.
c-- We keep adding nisource/2 until the radiation field converges.
c-- The doppler losses are then equal to nisource/2.  We give two full
c-- nisource shots to bring the radiation field up to speed and help
c-- converge more quickly.
       if(it<1 .and. it>tsp_itrestart+1) then
        gas_decaygamma = gas_decaygamma*.5
        gas_decaybeta = gas_decaybeta*.5
       endif
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
c-- update density, start temperature derivative
c===============================================
      gas_rho = gas_mass/gas_vol
c-- temperature
      gas_ur = pc_acoef*gas_temp**4
c-- temperature truncated by opacity table bounds
      help = pc_ev/pc_kb
      temp = max(0.01d0*help,min(5d0*help,gas_temp))
c
c-- sanity check temperatures
      if(any(gas_temp/=gas_temp)) stop 'gas_temp NaN'
      if(any(gas_temp<=0d0)) stop 'gas_temp<=0'
c
c
c
c-- compute the starting tempurature derivative in the fleck factor
      if(lfirst .or. in_opacanaltype/='none') then
c-- temporarily change!{{{
       gas_temp = dtempfrac*gas_temp
       if(in_opacanaltype=='none') then
        if(.not.in_noeos) call eos_update(.false.,dtempfrac*temp)
       endif
c
       if(in_opacanaltype/='none') then
        call analytic_opacity
       else
c-- initialize opacity
        gas_sig = 0d0
        gas_cap = 0.0
        if(in_doemiss) then
          if(.not.allocated(gas_em_cap)) then
            stop 'in_doemiss and gas_em_cap not allocated'
          endif
          if(.not.allocated(gas_em_capgrey)) then
            stop 'in_doemiss and gas_em_capgrey not allocated'
          endif
          gas_em_cap = 0.0
        endif
        if(.not.in_notbopac) then
          call tabular_opacity(in_doemiss)
        else
          call physical_opacity(dtempfrac*temp)
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
       if (in_doemiss) then
         capgreyalt = gas_em_capgrey/gas_rho !per gram
       else
         capgreyalt = gas_capgrey/gas_rho !per gram
       endif
c
c-- change back
       gas_temp = gas_temp/dtempfrac
!}}}
      endif
c
c
c
c-- solve LTE EOS, heat capacity
c===============================
      do_output = (in_io_pdensdump=='each' .or.
     &  (in_io_pdensdump=='one' .and. tsp_it==1))
      if(.not.in_noeos) call eos_update(do_output,temp)
c
      if(in_gas_cvcoef>0d0) then
c-- calculate power law heat capacity
       gas_bcoef = in_gas_cvcoef * gas_temp**in_gas_cvtpwr *
     &   gas_rho**in_gas_cvrpwr
      else
c-- calculate physical heat capacity
       if(in_noeos) then
         gas_bcoef = 1.5d0*pc_kb*gas_natom/gas_vol
       else
c-- add free electron contribution if physical EOS was used
         gas_bcoef = 1.5d0*pc_kb*(1d0+gas_nelec)*gas_natom /
     &      gas_vol
       endif
      endif
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
c-- calculate opacities
c======================
c-- gamma opacity
      if(in_srctype=='tabl'.or.in_srctype=='rpro') then
        gas_capgam = in_opcapgam*gas_rho
      else
        gas_capgam = in_opcapgam*gas_ye*gas_rho
      endif
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
c-- initialize opacity
        gas_sig = 0d0
        gas_cap = 0.0
        if(in_doemiss) then
          if(.not.allocated(gas_em_cap)) then
            stop 'in_doemiss and gas_em_cap not allocated'
          endif
          if(.not.allocated(gas_em_capgrey)) then
            stop 'in_doemiss and gas_em_capgrey not allocated'
          endif
          gas_em_cap = 0.0
        endif
c-- calculate opacities
        if(.not.in_notbopac) then
          call tabular_opacity(in_doemiss)
        else
          call physical_opacity(temp)
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
      if(trim(in_io_opacdump)=='off') then !{{{
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
c-- Calculating Fleck factor
      call fleck_factor(tempalt,capgreyalt)
c
c
c-- save previous values for gentile-fleck factor calculation in next iter
      tempalt = gas_temp
      if (in_doemiss) then
        capgreyalt = gas_em_capgrey/gas_rho !per gram
      else
        capgreyalt = gas_capgrey/gas_rho !per gram
      endif
c
      lfirst = .false.
c
c-- clean up
      if(tsp_it==tsp_nt .and. tsp_it>1 .and. allocated(tempalt)) then
       deallocate(tempalt,capgreyalt)
      endif
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
      use nucdatamod
      use gasmod
      implicit none
      real*8,intent(in) :: t
************************************************************************
* recalculate natom for nuclear decay from scratch
************************************************************************
      integer :: i
      real*8 :: help
      real*8 :: x(gas_ncell,0:2)
      real*8 :: natom(gas_ncell)
      real*8 :: dye(gas_ncell) !delta natom*ye
c
c-- save norm for conservation check
      natom = sum(gas_natom1fr(22:28,:),dim=1)
c
c-- zero
      gas_natom1fr(22:28,:) = 0d0
      dye = 0d0
c
c-- ni56->co56->fe56
c-- initial radioactive natom
      x(:,2) = gas_natom0fr(-2,:,1)
      x(:,1) = gas_natom0fr(-1,:,1)
      x(:,0) = 0d0
c-- subtract original ye*natom
      dye = dye - x(:,2)*(28d0/56)
      dye = dye - x(:,1)*(27d0/56)
c-- update
      call nucdecay3(gas_ncell,t,nuc_thl_ni56,nuc_thl_co56,x)
c-- current radioactive natom
      gas_natom1fr(gas_ini56,:) = x(:,2)
      gas_natom1fr(gas_ico56,:) = x(:,1)
c-- add decayed fraction to total natom
      forall(i=0:2) gas_natom1fr(26+i,:) = gas_natom1fr(26+i,:) + x(:,i)
c-- add current ye*natom
      dye = dye + x(:,2)*(28d0/56)
      dye = dye + x(:,1)*(27d0/56)
      dye = dye + x(:,0)*(26d0/56)
c
c-- fe52->mn52->cr52
c-- initial radioactive natom
      x(:,2) = gas_natom0fr(-2,:,2)
      x(:,1) = gas_natom0fr(-1,:,2)
      x(:,0) = 0d0
c-- subtract original ye*natom
      dye = dye - x(:,2)*(26d0/52)
      dye = dye - x(:,1)*(25d0/52)
c-- update
      call nucdecay3(gas_ncell,t,nuc_thl_fe52,nuc_thl_mn52,x)
c-- current radioactive natom
      gas_natom1fr(gas_ife52,:) = x(:,2)
      gas_natom1fr(gas_imn52,:) = x(:,1)
c-- add decayed fraction to total natom
      forall(i=0:2) gas_natom1fr(24+i,:) = gas_natom1fr(24+i,:) + x(:,i)
c-- add current ye*natom
      dye = dye + x(:,2)*(26d0/52)
      dye = dye + x(:,1)*(25d0/52)
      dye = dye + x(:,0)*(24d0/52)
c
c-- cr48->v48->ti48
c-- initial radioactive natom
      x(:,2) = gas_natom0fr(-2,:,3)
      x(:,1) = gas_natom0fr(-1,:,3)
      x(:,0) = 0d0
c-- subtract original ye*natom
      dye = dye - x(:,2)*(24d0/48)
      dye = dye - x(:,1)*(23d0/48)
c-- update
      call nucdecay3(gas_ncell,t,nuc_thl_cr48,nuc_thl_v48,x)
c-- current radioactive natom
      gas_natom1fr(gas_icr48,:) = x(:,2)
      gas_natom1fr(gas_iv48,:) = x(:,1)
c-- add decayed fraction to total natom
      forall(i=0:2) gas_natom1fr(22+i,:) = gas_natom1fr(22+i,:) + x(:,i)
c-- add current ye*natom
      dye = dye + x(:,2)*(24d0/48)
      dye = dye + x(:,1)*(23d0/48)
      dye = dye + x(:,0)*(22d0/48)
c
c-- add stable fraction to total natom
      forall(i=1:2) gas_natom1fr(26+i,:) = gas_natom1fr(26+i,:) +
     &  gas_natom0fr(i,:,1)
      forall(i=1:2) gas_natom1fr(24+i,:) = gas_natom1fr(24+i,:) +
     &  gas_natom0fr(i,:,2)
      forall(i=0:2) gas_natom1fr(22+i,:) = gas_natom1fr(22+i,:) +
     &  gas_natom0fr(i,:,3)
c
c-- natom conservation check
      do i=1,gas_ncell
       help = sum(gas_natom1fr(22:28,i))
       if(abs(help-natom(i))>1d-14*natom(i)) stop
     &   'update_natomfr: natom not conserved'
      enddo
c
c-- calculate Ye
      gas_ye = gas_ye0 + dye
c!}}}
      end subroutine update_natomfr
c vim: fdm=marker
