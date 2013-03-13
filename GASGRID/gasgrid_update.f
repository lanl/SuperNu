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
      IMPLICIT NONE
************************************************************************
* Update the part of the gas grid that depends on time and temperature.
* The non-changing part is computed in gasgrid_setup.
* The work done here is:
* - nuclear decay (energy release and chemical composition update)
* - temperature and volume
* - LTE EOS: ionization balance and electron density
* - opacities
************************************************************************
      LOGICAL :: do_output
      integer :: i,iw,j,ir,k
      REAL*8 :: help
      REAL*8 :: natom1fr(gas_nr,-2:-1) !todo: memory storage order?
      REAL*8 :: natom2fr(gas_nr,-2:-1)
      REAL*8 :: capbcum(gas_ng)
c-- gamma opacity
      REAL*8,parameter :: ye=.5d0 !todo: compute this value
c-- thomson scattering
      REAL*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- distribute packets
      integer :: mpacket !# packets to generate on each mpi rank
      integer :: nlower  !# ranks with 1 packet less
      REAL*8 :: enemit(gas_nr)
      REAL*8 :: chiross(gas_nr),capplanck(gas_nr)
c-- timing
      real :: t0,t1
c
c-- begin
c
      write(7,*)
      write(7,*) 'update gas grid:'
      write(7,*) '---------------------------'
      if(tsp_tn==1) then
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
      call update_natomfr(tsp_texp)
      forall(i=-2:-1) natom1fr(:,i) = gas_vals2(:)%natom1fr(i)
      call update_natomfr(tsp_texp + tsp_dt)
      forall(i=-2:-1) natom2fr(:,i) = gas_vals2(:)%natom1fr(i)
c
c-- update the abundances for the center time
      call update_natomfr(tsp_tcenter)
c
c-- energy deposition
!     gas_vals2(:)%nisource = gas_vals2(:)%natom*(
!    &  (natom1fr(:,gas_ini56) - natom2fr(:,gas_ini56)) *
!    &   (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
!    &  (natom1fr(:,gas_ico56) - natom2fr(:,gas_ico56))*pc_qhl_co56)
!}}}
c
c
c
c-- update temperature and volume
c================================
      gas_temphist(:,tsp_tn) = gas_vals2%temp
      gas_vals2%vol = gas_vals2%volr*(gas_velout*tsp_tcenter)**3 !volume in cm^3
      gas_vals2%volcrp = gas_vals2%vol !effective volume in cm^3
c
c
c
c-- solve LTE EOS
c================
      do_output = (in_pdensdump=='each' .or. !{{{
     &  (in_pdensdump=='one' .and. tsp_tn==1))
c
      call eos_update(do_output)
      if(tsp_tn==1) write(6,'(1x,a27,2(f8.2,"s"))')
     &  'eos timing                :',t_eos!}}}
c
c
c
c-- opacity per rcell unit
c========================
      calc_opac: if(tsp_tn==0) then
c!{{{
c-- gamma opacity
       gas_capgam = in_opcapgam*ye*
     &   gas_vals2(:)%mass/gas_vals2(:)%volcrp
c!}}}
      else calc_opac !tsp_tn
c!{{{
c-- thomson scattering
       gas_sig = cthomson*gas_vals2(:)%nelec*
     &   gas_vals2(:)%natom/gas_vals2(:)%volcrp
c
c-- 
       call opacity_calculation
c
c-- write out opacities (additional gray opacity not included!)
c--------------------------------------------------------------
       if(trim(in_opacdump)=='off') then !{{{
c-- do nothing
       else
c-- open file descriptor
        if(tsp_tn==1) then
         open(4,file='opacdump',action='write',status='replace')
c-- write wl-grid
         write(4,'(a,3i8)') '#',gas_ng,gas_nr
         write(4,'(a,3i8,1p,2e12.4)') '#',tsp_tn,tsp_tn,0, 0., 0.
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
         write(4,'(a,3i8,1p,2e12.4)') '#',tsp_tn,tsp_tn,ir,
     &     gas_vals2(ir)%temp,gas_sig
         write(4,'(1p,10e12.4)') (gas_cap(ir,j),j=1,gas_ng)
        enddo
c-- close file
        close(4)!}}}
       endif !do_output
c
c-- add some gray opacity, for testing
       if(in_opcap>0.) then
        forall(i=1:gas_ng) gas_cap(:,i) = gas_cap(:,i) + in_opcap*
     &    gas_vals2(:)%mass/gas_vals2(:)%volcrp
       endif
c
c
c-- Rosseland opacity
c-- normalization integral first
       forall(ir=1:gas_nr)
     &  chiross(ir) = sum(dplanckdtemp(gas_wl,gas_vals2(ir)%temp)*
     &    gas_dwl)
       forall(ir=1:gas_nr)
     &  capplanck(ir) = sum(planck(gas_wl,gas_vals2(ir)%temp)*
     &    gas_dwl)
c-- check against analytic solution
c      write(7,'(1p,10e12.4)') (chiross(ir),ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (4/pi*sb*gas_vals2(ir)%temp**3,ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (capplanck(ir),ir=1,gas_nr)
c      write(7,'(1p,10e12.4)') (sb/pi*gas_vals2(ir)%temp**4,ir=1,gas_nr)
c-- now the opacity weighting integral
       forall(ir=1:gas_nr)
     &  chiross(ir) = chiross(ir) /
     &    sum(dplanckdtemp(gas_wl,gas_vals2(ir)%temp)*gas_dwl/
     &    (gas_cap(ir,:) + gas_sig(ir)))
       forall(ir=1:gas_nr)
     &  capplanck(ir) = sum(planck(gas_wl,gas_vals2(ir)%temp)*
     &    gas_cap(ir,:)*gas_dwl) / capplanck(ir)
c-- Rosseland output
       write(7,*) 'mean opacities:'
       write(7,'(a8,7a12)') 'ir',
     &   'tau_Ross','mfp_Ross','chi Ross','cap_B','  cap_B',
     &   ' sig','    sig'
       write(7,'(a8,7a12)') 'units:',
     &   '     [1]','    [cu]',' [cm^-1]',' [cu]','[cm^-1]',
     &   '[cu]','[cm^-1]'
       do ir=1,gas_nr
        write(7,'(i8,1p,7e12.4)') ir,
     &    sum(chiross(ir:)),1d0/chiross(ir),chiross(ir),
     &    capplanck(ir),capplanck(ir),
     &    gas_sig(ir),gas_sig(ir)
       enddo
c
c-- timing output
       if(tsp_tn==1 .and. tsp_tn==1)
     &   write(6,'(1x,a27,3(f8.2,"s"))') 'opacity timing: bb|bf|ff  :',
     &   t_bb(1),t_bf(1),t_ff(1) !}}}
c
      endif calc_opac !tsp_tn
c
c
c
c-- output
c=========
c-- to stdout!{{{
c
c-- energy depots
      if(tsp_tn==1) then
       write(6,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &   sum(gas_vals2(:)%nisource)
      endif !tsp_tn
c-- totals
      write(7,*)
      write(7,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &  sum(gas_vals2(:)%nisource)
c-- arrays
*     write(7,'(a6,5a12)')'ir','edep/vol','enostor/vol','rho',
      write(7,'(a6,5a12)')'ir','edep/dt','rho',
     &  'nelec','volcrp/vol'
      do i=1,gas_nr,10
       write(7,'(i6,1p,5e12.4)') (j,
     &  gas_vals2(j)%nisource/tsp_dt,
     &  gas_vals2(j)%mass/gas_vals2(j)%vol,
     &  gas_vals2(j)%nelec,gas_vals2(j)%volcrp/gas_vals2(j)%vol,
     &  j=i,min(i+9,gas_nr))
      enddo
!c
!c-- scattering coefficients
!      if(tsp_tn>0) then
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
c     -------------------------------
      use physconstmod
      use gasgridmod
      use inputparmod
      IMPLICIT NONE
      REAL*8,intent(in) :: tsince
************************************************************************
* update natom fractions for nuclear decay
************************************************************************
      REAL*8 :: expni,expco,help
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
c
      end subroutine update_natomfr
