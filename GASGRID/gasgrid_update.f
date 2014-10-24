      subroutine gasgrid_update
c     -----------------------
      use physconstmod
      use miscmod, only:warn
      use ionsmod
      use timestepmod
      use gasgridmod
      use inputparmod
      use timingmod
      use profiledatamod
      implicit none
************************************************************************
* Update the part of the gas grid that depends on time and temperature.
* The non-changing part is computed in gasgrid_setup.
* The work done here is:
* - nuclear decay (energy release and chemical composition update)
* - temperature and volume
* - LTE EOS: ionization balance and electron density
* - opacities
************************************************************************
      logical :: do_output,lexist,planckcheck
      integer :: i,j,k,l,ig,it,istat
      real*8 :: help,x1,x2
      real*8,external :: specint
      real*8 :: dtempfrac = 0.99d0
      real*8 :: natom1fr(gas_nx,gas_ny,gas_nz,-2:-1) !todo: memory storage order?
      real*8 :: natom2fr(gas_nx,gas_ny,gas_nz,-2:-1)
c-- gamma opacity
      real*8,parameter :: ye=.5d0 !todo: compute this value
c-- timing
      real*8 :: t0,t1
c
c-- begin
      call time(t0)
c
c-- nuclear decay
c================
c-- Get ni56 and co56 abundances on begin and end of the time step.!{{{
c-- The difference between these two has decayed.
      if(gas_isvelocity.and.gas_srctype=='none') then
c-- beginning of time step
       help = tsp_t
       call update_natomfr(help)
       forall(l=-2:-1) natom1fr(:,:,:,l) = gas_vals2%natom1fr(l)
c-- end of time step
       call update_natomfr(tsp_t + tsp_dt)
       forall(l=-2:-1) natom2fr(:,:,:,l) = gas_vals2%natom1fr(l)
c
c-- update the abundances for the center time
       !call update_natomfr(tsp_tcenter)
       call update_natomfr(tsp_t)
c
c-- energy deposition
       gas_vals2%nisource =  !per average atom (mix of stable and unstable)
     &   (natom1fr(:,:,:,gas_ini56) - natom2fr(:,:,:,gas_ini56)) *
     &     (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(:,:,:,gas_ico56) - natom2fr(:,:,:,gas_ico56)) *
     &     pc_qhl_co56
c-- total, units=ergs
       gas_vals2%nisource = gas_vals2%nisource *
     &   gas_vals2%natom
c-- use gamma deposition profiles if data available
       if(prof_ntgam>0) then
        help = sum(gas_vals2%nisource)
!       write(6,*) 'ni56 source:',help
        if(gas_ny>1 .or. gas_nz>1) stop 'gg_update: gam_prof: no 2D/3D'
        gas_vals2(:,1,1)%nisource = help * gamma_profile(tsp_t)
       endif
      endif
!}}}
c
c
c
c-- update volume and density 
c============================
      if(gas_isvelocity) then!{{{
       help = gas_velout*tsp_t
      else
       if(gas_ny>1) stop 'gg_update: help: no 2D'
       help = gas_lx
      endif
      !gas_vals2%vol = gas_vals2%volr*(gas_velout*tsp_tcenter)**3 !volume in cm^3
      gas_vals2%vol = gas_vals2%volr*help**3 !volume in cm^3
      gas_vals2%volcrp = gas_vals2%vol !effective volume in cm^3
c
c-- density
      gas_vals2%rho = gas_vals2%mass/gas_vals2%vol
c!}}}
c
c
c-- update interpolated density and temperatures at cell edges
c=============================================================
!Calculating power law heat capacity
      gas_vals2%bcoef = in_cvcoef * gas_temp**in_cvtpwr *
     &  gas_vals2%rho**in_cvrpwr

c-- add initial thermal input to gas_eext
      if(tsp_it==1) then
       gas_eext = sum(gas_vals2%bcoef*gas_temp*gas_vals2%vol)
      endif
c
c
!     return !DEBUG
c
c
c
c-- opacity
c==========
c!{{{
c
c-- compute the starting tempurature derivative in the fleck factor
      if(tsp_it==1.or.in_opacanaltype/='none') then
       gas_temp=dtempfrac*gas_temp
       if(gas_isvelocity .and. in_opacanaltype=='none') then
        call eos_update(.false.)
       endif
c
       call analytic_opacity
       if(in_opacanaltype=='none') then
        if(in_ngs==0) then
         call physical_opacity
        else
         call physical_opacity_subgrid
        endif
       endif
c
c-- gamma opacity
       gas_capgam = in_opcapgam*ye*
     &   gas_vals2%mass/gas_vals2%volcrp
c
       gas_siggreyprevit = gas_siggrey
       gas_temp = gas_temp/dtempfrac
      endif
c
c-- solve LTE EOS
c================
      if(gas_isvelocity) then
       do_output = (in_pdensdump=='each' .or. !{{{
     &   (in_pdensdump=='one' .and. tsp_it==1))
c
       call eos_update(do_output)
      endif
c
c
c-- simple physical group/grey opacities: Planck and Rosseland 
      call analytic_opacity
c-- add physical opacities
c-- rtw: must avoid reset in group_opacity routine
      if(in_opacanaltype=='none') then
c-- test existence of input.opac file
       inquire(file='input.opac',exist=lexist)
       if(lexist) then
c-- read in opacities
        open(4,file='input.opac',status='old',iostat=istat)!{{{
        if(istat/=0) stop 'read_opac: no file: input.opac'
c-- read header
        read(4,*,iostat=istat)
        if(istat/=0) stop 'read_opac: file empty: input.opac'
c-- read each cell individually
        do it=1,tsp_it
c-- skip delimiter
         read(4,*,iostat=istat)
         if(istat/=0) stop 'read_opac: delimiter error: input.opac'
c-- read data
         if(gas_ny>1 .or. gas_nz>1) stop 'read_opac: no 2D/3D'
         do i=1,gas_nx
          read(4,*,iostat=istat) help,gas_sig(i,1,1),gas_cap(:,i,1,1)
          if(istat/=0) stop 'read_opac: body error: input.opac'
         enddo !i
        enddo !it
        close(4)
        write(6,*) 'read_opac: read successfully'
!}}}
       elseif(in_ngs==0) then
c-- calculate opacities
        call physical_opacity
       else
        call physical_opacity_subgrid
       endif
c
      endif
c
c-- Planck opacity
      planckcheck = (.not.in_nobbopac .or. .not.in_nobfopac .or.
     &  .not.in_noffopac)
!Ryan,why is this conditional (drr 14/05/31)?
      if(planckcheck) then
       gas_siggrey = 0d0
       do k=1,gas_nz
       do j=1,gas_ny
       do i=1,gas_nx
        do ig=1,gas_ng
         x1 = pc_h*pc_c/(gas_wl(ig + 1)*pc_kb*gas_temp(i,j,k))
         x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(i,j,k))
         gas_siggrey(i,j,k) = gas_siggrey(i,j,k)+
     &     15d0*gas_cap(ig,i,j,k)*specint(x1,x2,3)/pc_pi**4
        enddo
       enddo !i
       enddo !j
       enddo !k
      endif
      !write(*,*) gas_siggrey(1)
      !write(*,*) gas_cap(:,1)
      !gas_siggrey(:)=0.5*gas_cap(2,:)
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
       if(gas_ny>1 .or. gas_nz>1) stop 'output.opac: no 2D/3D'
       if(tsp_it==1) write(4,'("#",3i8)') gas_ng,gas_nx,tsp_nt
       write(4,'("#",3i8)') tsp_it
c-- body
       do i=1,gas_nx
        write(4,'(1p,9999e12.4)') gas_temp(i,1,1),gas_sig(i,1,1),
     &    (gas_cap(j,i,1,1),j=1,gas_ng)
       enddo
c-- close file
       close(4)
      endif !do_output !}}}
c
c-- Calculating Fleck factor, leakage opacities
      call fleck_factor(dtempfrac)
c-- Calculating emission probabilities for each group in each cell
      call emission_probability
c-- Calculating IMC-DDMC albedo coefficients and DDMC leakage opacities
      call leakage_opacity
c
      call time(t1)
      call timereg(t_gasupd,t1-t0)
c
      end subroutine gasgrid_update
c
c
c
      subroutine update_natomfr(tsince)
c     -------------------------------!{{{
      use physconstmod
      use gasgridmod
      use inputparmod
      implicit none
      real*8,intent(in) :: tsince
************************************************************************
* update natom fractions for nuclear decay
************************************************************************
      real*8 :: expni,expco,help
c
      expni = exp(-tsince/pc_thl_ni56)
      expco = exp(-tsince/pc_thl_co56)
c
c-- update Fe
      help = 1d0 + (pc_thl_co56*expco - pc_thl_ni56*expni)/
     &  (pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Fe < 0'
      gas_vals2%natom1fr(26) = gas_vals2%natom0fr(gas_ini56)*help+!initial Ni56
     &  gas_vals2%natom0fr(gas_ico56)*(1d0-expco) +                  !initial Co56
     &  gas_vals2%natom0fr(0)                                        !initial Fe (stable)
c
c-- update Co56 and Co
      help = pc_thl_co56*(expni - expco)/(pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Co < 0'
c-- Co56
      gas_vals2%natom1fr(gas_ico56) =
     &  gas_vals2%natom0fr(gas_ini56)*help +  !initial Ni56
     &  gas_vals2%natom0fr(gas_ico56)*expco   !initial Co56
c-- Co
      gas_vals2%natom1fr(27) = gas_vals2%natom1fr(gas_ico56) +  !unstable
     &  gas_vals2%natom0fr(1)                                      !initial Co (stable)
c
c-- update Ni56 and Ni
c-- Ni56
      gas_vals2%natom1fr(gas_ini56) =
     &  gas_vals2%natom0fr(gas_ini56)*expni  !initial Ni56
c-- Ni
      gas_vals2%natom1fr(28) = gas_vals2%natom1fr(gas_ini56) + !unstable
     &  gas_vals2%natom0fr(2)                              !initial Ni (stable)
c!}}}
      end subroutine update_natomfr
