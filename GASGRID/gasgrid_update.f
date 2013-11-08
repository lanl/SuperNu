      subroutine gasgrid_update
c     -----------------------
      use mpimod, only:nmpi
      use physconstmod
      use miscmod, only:warn
      use ionsmod
      use timestepmod
      use gasgridmod
      use inputparmod
      use timingmod
      use gammaprofmod
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
      logical :: do_output
      integer :: i,j,ir
      real*8 :: help
      real*8 :: dwl(gas_ng)
      real*8 :: natom1fr(gas_nr,-2:-1) !todo: memory storage order?
      real*8 :: natom2fr(gas_nr,-2:-1)
c-- gamma opacity
      real*8,parameter :: ye=.5d0 !todo: compute this value
c-- distribute packets
      real*8 :: chiross(gas_nr),capplanck(gas_nr)
c-- timing
      real :: t0,t1
c
c-- begin
c
      write(7,*)
      write(7,*) 'update gas grid:'
      write(7,*) '---------------------------'
      if(tsp_it==1) then
       write(6,*)
       write(6,*) 'update gas grid:'
       write(6,*) '---------------------------'
      endif
c
c
      call time(t0)
c
c
c-- nuclear decay
c================
c-- Get ni56 and co56 abundances on begin and end of the time step.!{{{
c-- The difference between these two has decayed.
      if(gas_isvelocity.and.gas_srctype=='none') then
c-- beginning of time step
       help = tsp_texp
       if(tsp_it==1) help = 0d0 !start accounting from day 0
       call update_natomfr(help)
       forall(i=-2:-1) natom1fr(:,i) = gas_vals2(:)%natom1fr(i)
c-- end of time step
       call update_natomfr(tsp_texp + tsp_dt)
       forall(i=-2:-1) natom2fr(:,i) = gas_vals2(:)%natom1fr(i)
c
c-- update the abundances for the center time
       !call update_natomfr(tsp_tcenter)
       call update_natomfr(tsp_texp)
c
c-- energy deposition
       gas_vals2(:)%nisource =  !per average atom (mix of stable and unstable)
     &   (natom1fr(:,gas_ini56) - natom2fr(:,gas_ini56)) *
     &    (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(:,gas_ico56) - natom2fr(:,gas_ico56))*pc_qhl_co56
c-- total, units=ergs
       gas_vals2(:)%nisource = gas_vals2(:)%nisource *gas_vals2(:)%natom
c-- use gamma deposition profiles if data available
       if(gamprf_nt>0) then
        help = sum(gas_vals2%nisource)
        i = min(tsp_it,gamprf_nt)
        gas_vals2(:)%nisource = help * gamprf_prof(:,i)
        write(6,'(a,2f5.2)') 'using gamma prof: t,int=',gamprf_time(i),
     &    sum(gamprf_prof(:,i))
       endif
      endif
!}}}
c
c
c
c-- update volume and density 
c============================
      if(gas_isvelocity) then!{{{
       help = gas_velout*tsp_texp
      else
       help = gas_l0+gas_lr
      endif
      !gas_vals2%vol = gas_vals2%volr*(gas_velout*tsp_tcenter)**3 !volume in cm^3
      gas_vals2%vol = gas_vals2%volr*help**3 !volume in cm^3
      gas_vals2%volcrp = gas_vals2%vol !effective volume in cm^3
c
c-- density
      gas_vals2%rho = gas_vals2%mass/gas_vals2%vol
c
c-- keep track of temperature evolution
      gas_temphist(:,tsp_it) = gas_temp!}}}
c
c
c-- update interpolated density and temperatures at cell edges
c=============================================================
!Interpolating cell boundary temperatures (in keV currently): loop!{{{
      if(gas_isshell) then
       gas_tempb(1) = gas_templ0
      else
       gas_tempb(1) = gas_temp(1)
      endif
!gas_tempb(1) = 1.0
      do ir=2,gas_nr
       gas_tempb(ir) = (gas_temp(ir)**4 +
     &   gas_temp(ir-1)**4)/2.0
       gas_tempb(ir) = gas_tempb(ir)**0.25
      enddo
      gas_tempb(gas_nr+1) = gas_temp(gas_nr)
!Interpolating cell boundary densities (in g/cm^3): loop
      gas_rhob(1) = gas_vals2(1)%rho
      do ir = 2, gas_nr
       !gas_rhob(ir)=(gas_vals2(ir)%rho*gas_vals2(ir)%vol+ &
       !     gas_vals2(ir-1)%rho*gas_vals2(ir-1)%vol)/ &
       !     (gas_vals2(ir)%vol+gas_vals2(ir-1)%vol)
       !gas_rhob(ir) = (gas_vals2(ir)%rho*gas_vals2(ir-1)%rho)**0.5d0
         gas_rhob(ir)=1d0/(.5d0/gas_vals2(ir-1)%rho+
     & .5d0/gas_vals2(ir)%rho)
      enddo
      gas_rhob(gas_nr+1) = gas_vals2(gas_nr)%rho

!Calculating power law heat capacity
      do ir=1,gas_nr
       gas_vals2(ir)%bcoef = gas_cvcoef *
     &   gas_temp(ir)**gas_cvtpwr *
     &   gas_vals2(ir)%rho**gas_cvrpwr
      enddo!}}}
c
c
c
c-- reset counters
c=================
      gas_erad = 0.0   !Total radiation energy

      gas_eraddens =0d0 !radiation density field
      gas_vals2%eraddens =0d0
c
!     return !DEBUG
c
c
c-- solve LTE EOS
c================
      if(gas_isvelocity) then
         do_output = (in_pdensdump=='each' .or. !{{{
     &        (in_pdensdump=='one' .and. tsp_it==1))
c
            call eos_update(do_output)
            if(tsp_it==1) write(6,'(1x,a27,2(f8.2,"s"))')
     &           'eos timing                :',t_eos !}}}
      endif
c     
c
c
c-- opacity per rcell unit
c========================
      calc_opac: if(tsp_it==0) then
c!{{{
c-- gamma opacity
       gas_capgam = in_opcapgam*ye*
     &   gas_vals2(:)%mass/gas_vals2(:)%volcrp
c!}}}
      else calc_opac !tsp_it
c!{{{
c-- simple physical group/grey opacities: Planck and Rosseland 
       call analytic_opacity
c-- add physical opacities
       call physical_opacity
       !write(*,*) gas_siggrey(1)
       !write(*,*) gas_cap(:,1)
       !gas_siggrey(:)=0.5*gas_cap(2,:)
c
c-- write out opacities (additional gray opacity not included!)
c--------------------------------------------------------------
       if(trim(in_opacdump)=='off') then !{{{
c-- do nothing
       else
c-- open file descriptor
        if(tsp_it==1) then
         open(4,file='opacdump',action='write',status='replace')
c-- write wl-grid
         write(4,'(a,3i8)') '#',gas_ng,gas_nr
         write(4,'(a,3i8,1p,2e12.4)') '#',tsp_it,tsp_it,0, 0., 0.
         write(4,'(1p,10e12.4)') gas_wl
        elseif(trim(in_opacdump)=='each') then
         open(4,file='opacdump',action='write',status='replace')
        elseif(trim(in_opacdump)=='all') then
         open(4,file='opacdump',action='write',status='old',
     &     position='append')
        endif
       endif !off
c
c-- write opacity grid
       inquire(4,opened=do_output)
       if(do_output) then
        do ir=1,gas_nr
c-- convert from opacity in redona's rcell units to opacity per cm
         write(4,'(a,3i8,1p,2e12.4)') '#',tsp_it,tsp_it,ir,
     &     gas_temp(ir),gas_sig
         write(4,'(1p,10e12.4)') (gas_cap(j,ir),j=1,gas_ng)
        enddo
c-- close file
        close(4)
       endif !do_output !}}}
c
c-- Calculating Fleck factor, leakage opacities
       call fleck_factor
c-- Calculating emission probabilities for each group in each cell
       call emission_probability
c-- Calculating IMC-DDMC albedo coefficients and DDMC leakage opacities
       call leakage_opacity
c
c
c-- Rosseland opacity
c-- normalization integral first
       dwl = gas_wl(2:) - gas_wl(:gas_ng)
       forall(ir=1:gas_nr)
     &  chiross(ir) = sum(dplanckdtemp(gas_wl,gas_temp(ir))*dwl)
       forall(ir=1:gas_nr)
     &  capplanck(ir) = sum(planck(gas_wl,gas_temp(ir))*dwl)
c-- check against analytic solution
c      write(7,'(1p,10e12.4)') (chiross(ir),ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (4/pi*sb*gas_temp(ir)**3,ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (capplanck(ir),ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (sb/pi*gas_temp(ir)**4,ir=1,gas_nr)
c-- now the opacity weighting integral
       forall(ir=1:gas_nr)
     &  chiross(ir) = chiross(ir) /
     &    sum(dplanckdtemp(gas_wl,gas_temp(ir))*dwl/
     &    (gas_caprosl(:,ir) + gas_sig(ir)))
       forall(ir=1:gas_nr)
     &  capplanck(ir) = sum(planck(gas_wl,gas_temp(ir))*
     &    gas_cap(:,ir)*dwl) / capplanck(ir)
c-- Rosseland output
       write(7,*) 'mean opacities:'
       write(7,'(a8,3(3a11,1x))') 'ir','chi_Ross','cap_B','sig',
     &   'cap_min','cap_max','cap_mean',
     &   'capr_min','capr_max','capr_mean'
       do ir=1,gas_nr
        write(7,'(i8,1p,3(3e11.4,1x))') ir,
     &    chiross(ir),capplanck(ir),gas_sig(ir),
     &    minval(gas_cap(:,ir)),maxval(gas_cap(:,ir)),
     &     sum(gas_cap(:,ir))/gas_ng,
     &    minval(gas_caprosl(:,ir)),maxval(gas_caprosl(:,ir)),
     &     sum(gas_caprosl(:,ir))/gas_ng
       enddo
c
c-- timing output
       if(tsp_it==1 .and. tsp_it==1)
     &   write(6,'(1x,a27,3(f8.2,"s"))') 'opacity timing: bb|bf|ff  :',
     &   t_bb(1),t_bf(1),t_ff(1) !}}}
      endif calc_opac !tsp_it
c
c
c
c-- output
c=========
c-- to stdout!{{{
c
c-- energy depots
      if(tsp_it==1) then
       write(6,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &   sum(gas_vals2(:)%nisource)
      endif !tsp_it
c-- totals
      write(7,*)
      write(7,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &  sum(gas_vals2(:)%nisource)
c-- arrays
*     write(7,'(a6,4a12)')'ir','edep/vol','enostor/vol','rho',
      write(7,'(a6,4a12)')'ir','edep/dt','rho',
     &  'nelec','volcrp/vol'
      do i=1,gas_nr,10
       write(7,'(i6,1p,4e12.4)') (j,
     &  gas_vals2(j)%nisource/tsp_dt,
     &  gas_vals2(j)%mass/gas_vals2(j)%vol,
     &  gas_vals2(j)%nelec,gas_vals2(j)%volcrp/gas_vals2(j)%vol,
     &  j=i,min(i+9,gas_nr))
      enddo
!c
!c-- scattering coefficients
!      if(tsp_it>0) then
!       write(7,*)
!       write(7,*) 'sig'
!       write(7,'(1p,10e12.4)') gas_sig
!      endif!}}}
c
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
      gas_vals2(:)%natom1fr(26) = gas_vals2(:)%natom0fr(gas_ini56)*help+!initial Ni56
     &  gas_vals2(:)%natom0fr(gas_ico56)*(1d0-expco) +                  !initial Co56
     &  gas_vals2(:)%natom0fr(0)                                        !initial Fe (stable)
c
c-- update Co56 and Co
      help = pc_thl_co56*(expni - expco)/(pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Co < 0'
c-- Co56
      gas_vals2(:)%natom1fr(gas_ico56) =
     &  gas_vals2(:)%natom0fr(gas_ini56)*help +  !initial Ni56
     &  gas_vals2(:)%natom0fr(gas_ico56)*expco   !initial Co56
c-- Co
      gas_vals2(:)%natom1fr(27) = gas_vals2(:)%natom1fr(gas_ico56) +  !unstable
     &  gas_vals2(:)%natom0fr(1)                                      !initial Co (stable)
c
c-- update Ni56 and Ni
c-- Ni56
      gas_vals2(:)%natom1fr(gas_ini56) =
     &  gas_vals2(:)%natom0fr(gas_ini56)*expni  !initial Ni56
c-- Ni
      gas_vals2(:)%natom1fr(28) = gas_vals2(:)%natom1fr(gas_ini56) + !unstable
     &  gas_vals2(:)%natom0fr(2)                              !initial Ni (stable)
c!}}}
      end subroutine update_natomfr
