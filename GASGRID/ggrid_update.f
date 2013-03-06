      subroutine ggrid_update
c     -----------------------
      use mpimod, only:nmpi
      use physconstmod
      use miscmod, only:warn
      use ionsmod
      use timestepmod
      use ggridmod
      use inputparmod
      use rgridmod, only:rg_ncr
      use timingmod
      implicit none
************************************************************************
* Update the part of the gas grid that depends on time and temperature.
* The non-changing part is computed in ggrid_setup.
* The work done here is:
* - nuclear decay (energy release and chemical composition update)
* - temperature and volume
* - LTE EOS: ionization balance and electron density
* - opacities
************************************************************************
      logical :: do_output
      integer :: i,iw,j,icg,k
      real*8 :: help
      real*8 :: natom1fr(gg_ncg,-2:-1) !todo: memory storage order?
      real*8 :: natom2fr(gg_ncg,-2:-1)
      real*8 :: capbcum(in_nwlg)
c-- gamma opacity
      real*8,parameter :: ye=.5d0 !todo: compute this value
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- distribute packets
      integer :: mpacket !# packets to generate on each mpi rank
      integer :: nlower  !# ranks with 1 packet less
      real*8 :: enemit(gg_ncg)
      real*8 :: chiross(gg_ncg),capplanck(gg_ncg)
c-- timing
      real :: t0,t1
c
c-- begin
c
      write(7,*)
      write(7,*) 'update gas grid:'
      write(7,*) '---------------------------'
      if(tim_itc==1) then
       write(6,*)
       write(6,*) 'update gas grid:'
       write(6,*) '---------------------------'
      endif
c
c
      call time(t0)
c
c-- constants
c============
!bug  gg_cellength  = gg_vout*tim_cen/in_nr !found 04/apr/2011
      gg_cellength  = gg_vout*tim_cen/(in_nr + .5d0) !converts rcell length units to cm [cm/cell_unit]
c     rcellvol = pc_pi4/3.*(gg_vout*tim_cen)**3/rg_ncr
c
c
c
c-- nuclear decay
c================
c-- Get ni56 and co56 abundances on begin and end of the time step.!{{{
c-- The difference between these two has decayed.
      call update_natomfr(tim_bot)
      forall(i=-2:-1) natom1fr(:,i) = ggrid2(:)%natom1fr(i)
      call update_natomfr(tim_top)
      forall(i=-2:-1) natom2fr(:,i) = ggrid2(:)%natom1fr(i)
c
c-- update the abundances for the center time
      call update_natomfr(tim_cen)
c
c-- fill ggrid structure
c-- energy deposition (initialize local, replace with transport result later, in compute_local_engdep)
      if(tim_itc==0) then !gamma transport
c      ggrid(:)%engdep = dtim*(ggrid2(i)%natom1fr(gg_ini56)*pc_qhl_ni56/thl_ni56 + !rate
c    &   ggrid2(:)%natom1fr(gg_ico56)*pc_qhl_co56/thl_co56)*ggrid2(i)%natom
       ggrid(:)%engdep = ggrid2(:)%natom*(
     &   (natom1fr(:,gg_ini56) - natom2fr(:,gg_ini56)) *
     &    (pc_qhl_ni56 + pc_qhl_co56) +!ni56 that decays adds to co56
     &   (natom1fr(:,gg_ico56) - natom2fr(:,gg_ico56))*pc_qhl_co56)
      endif !tim_itc !}}}
c
c
c
c-- update temperature and volume
c================================
      if(any(ggrid2_temp(:,tim_itim)<=0d0)) then!{{{
       if(tim_itc>1) stop 'ggrid_update: temp==0 invalid'
       if(tim_itim==1) stop 'ggrid_update: temp==0 bad initialization'
c-- copy results from previous time-step
       ggrid2_temp(:,tim_itim) = ggrid2_temp(:,tim_itim-1)
      endif
c
      ggrid2%temp = ggrid2_temp(:,tim_itim)
      ggrid2%vol = ggrid2%volr*(gg_vout*tim_cen)**3 !volume in cm^3
      ggrid2%volcrp = pc_pi4/3d0*(gg_vout*tim_cen)**3 !effective volume in cm^3
     &  *ggrid%ncrp/rg_ncr!}}}
c
c
c
c-- solve LTE EOS
c================
      do_output = tim_itc>=in_ntc .and. !{{{
     &  (in_pdensdump=='each' .or.
     &  (in_pdensdump=='one' .and. tim_itim==1))
c
      call eos_update(do_output)
      if(tim_itc==1) write(6,'(1x,a27,2(f8.2,"s"))')
     &  'eos timing                :',t_eos!}}}
c
c
c
c-- opacity per rcell unit
c========================
      calc_opac: if(tim_itc==0) then
c!{{{
c-- gamma opacity
       ggrid(:)%capgam = in_opcapgam*ye*
     &   ggrid2(:)%mass/ggrid2(:)%volcrp*gg_cellength !gg_cellength converts cm^-1 to 1/rcell
c!}}}
      else calc_opac !tim_itc
c!{{{
c-- thomson scattering
       ggrid(:)%sig = cthomson*ggrid2(:)%nelec*
     &   ggrid2(:)%natom/ggrid2(:)%volcrp*gg_cellength 
c
c-- 
       call opacity_calculation
c
c-- write out opacities (additional gray opacity not included!)
c--------------------------------------------------------------
       if(trim(in_opacdump)=='off') then !{{{
c-- do nothing
       elseif(tim_itc>=in_ntc) then
c-- open file descriptor
        if(tim_itim==1) then
         open(4,file='opacdump',action='write',status='replace')
c-- write wl-grid
         write(4,'(a,3i8)') '#',in_nwlg,gg_ncg
         write(4,'(a,3i8,1p,2e12.4)') '#',tim_itim,tim_itc,0, 0., 0.
         write(4,'(1p,10e12.4)') ggrid_wl
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
        do icg=1,gg_ncg
c-- convert from opacity in redona's rcell units to opacity per gram
*        help = ggrid2(icg)%volcrp/(ggrid2(icg)%mass*gg_cellength)
c-- convert from opacity in redona's rcell units to opacity per cm
         help = 1d0/gg_cellength
         write(4,'(a,3i8,1p,2e12.4)') '#',tim_itim,tim_itc,icg,
     &     ggrid2(icg)%temp,ggrid(icg)%sig*help
         write(4,'(1p,10e12.4)') (ggrid_cap(icg,j)*help,j=1,in_nwlg)
        enddo
c-- close file
        close(4)!}}}
       endif !do_output
c
c-- add some gray opacity, for testing
       if(in_opcap>0.) then
        forall(i=1:in_nwlg) ggrid_cap(:,i) = ggrid_cap(:,i) + in_opcap*
     &    ggrid2(:)%mass/ggrid2(:)%volcrp*gg_cellength !convert from opacity in redona's rcell units to opacity per gram
       endif
c
c
c-- Rosseland opacity
c-- normalization integral first
       forall(icg=1:gg_ncg)
     &  chiross(icg) = sum(dplanckdtemp(ggrid_wl,ggrid2(icg)%temp)*
     &    ggrid2_dwl)
       forall(icg=1:gg_ncg)
     &  capplanck(icg) = sum(planck(ggrid_wl,ggrid2(icg)%temp)*
     &    ggrid2_dwl)
c-- check against analytic solution
c      write(7,'(1p,10e12.4)') (chiross(icg),icg=1,gg_ncg)
c      write(7,'(1p,10e12.4)') (4/pi*sb*ggrid2(icg)%temp**3,icg=1,gg_ncg)
c      write(7,'(1p,10e12.4)') (capplanck(icg),icg=1,gg_ncg)
c      write(7,'(1p,10e12.4)') (sb/pi*ggrid2(icg)%temp**4,icg=1,gg_ncg)
c-- now the opacity weighting integral
       forall(icg=1:gg_ncg)
     &  chiross(icg) = chiross(icg) /
     &    sum(dplanckdtemp(ggrid_wl,ggrid2(icg)%temp)*ggrid2_dwl/
     &    (ggrid_cap(icg,:) + ggrid(icg)%sig))
       forall(icg=1:gg_ncg)
     &  capplanck(icg) = sum(planck(ggrid_wl,ggrid2(icg)%temp)*
     &    ggrid_cap(icg,:)*ggrid2_dwl) / capplanck(icg)
c-- Rosseland output
       write(7,*) 'mean opacities:'
       write(7,'(a8,7a12)') 'icg',
     &   'tau_Ross','mfp_Ross','chi Ross','cap_B','  cap_B',
     &   ' sig','    sig'
       write(7,'(a8,7a12)') 'units:',
     &   '     [1]','    [cu]',' [cm^-1]',' [cu]','[cm^-1]',
     &   '[cu]','[cm^-1]'
       do icg=1,gg_ncg
        write(7,'(i8,1p,7e12.4)') icg,
     &    sum(chiross(icg:)),1d0/chiross(icg),chiross(icg)/gg_cellength,
     &    capplanck(icg),capplanck(icg)/gg_cellength,
     &    ggrid(icg)%sig,ggrid(icg)%sig/gg_cellength
       enddo
c
c-- timing output
       if(tim_itc==1 .and. tim_itim==1)
     &   write(6,'(1x,a27,3(f8.2,"s"))') 'opacity timing: bb|bf|ff  :',
     &   t_bb(1),t_bf(1),t_ff(1) !}}}
c
      endif calc_opac !tim_itc
c
c
c
c-- output
c=========
c-- to stdout!{{{
c
c-- energy depots
      if(tim_itc==1) then
       write(6,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &   sum(ggrid%engdep)
      endif !tim_itc
c-- totals
      write(7,*)
      write(7,'(1x,a,1p,e12.4)') 'energy deposition (Lagr)  :',
     &  sum(ggrid%engdep)
c-- arrays
*     write(7,'(a6,5a12)')'icg','engdep/vol','enostor/vol','rho',
      write(7,'(a6,5a12)')'icg','engdep/dt','rho',
     &  'nelec','volcrp/vol'
      do i=1,gg_ncg,10
       write(7,'(i6,1p,5e12.4)') (j,
     &  ggrid(j)%engdep/tim_dt,
     &  ggrid2(j)%mass/ggrid2(j)%vol,
     &  ggrid2(j)%nelec,ggrid2(j)%volcrp/ggrid2(j)%vol,
     &  j=i,min(i+9,gg_ncg))
      enddo
!c
!c-- scattering coefficients
!      if(tim_itc>0) then
!       write(7,*)
!       write(7,*) 'sig'
!       write(7,'(1p,10e12.4)') ggrid%sig
!      endif!}}}
c
c
      call time(t1)
      call timereg(t_ggupd,t1-t0)
c
      end subroutine ggrid_update
c
c
c
      subroutine update_natomfr(tsince)
c     -------------------------------
      use physconstmod
      use ggridmod
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
      ggrid2(:)%natom1fr(26) = ggrid2(:)%natom0fr(gg_ini56)*help + !initial Ni56
     &  ggrid2(:)%natom0fr(gg_ico56)*(1d0-expco) +                 !initial Co56
     &  ggrid2(:)%natom0fr(0)                                      !initial Fe (stable)
c
c-- update Co56 and Co
      help = pc_thl_co56*(expni - expco)/(pc_thl_ni56 - pc_thl_co56)
      if(help.lt.0) stop 'update_natomfr: Ni->Co < 0'
c-- Co56
      ggrid2(:)%natom1fr(gg_ico56) = ggrid2(:)%natom0fr(gg_ini56)*help +!initial Ni56
     &  ggrid2(:)%natom0fr(gg_ico56)*expco                              !initial Co56
c-- Co
      ggrid2(:)%natom1fr(27) = ggrid2(:)%natom1fr(gg_ico56) + !unstable
     &  ggrid2(:)%natom0fr(1)                              !initial Co (stable)
c
c-- update Ni56 and Ni
c-- Ni56
      ggrid2(:)%natom1fr(gg_ini56) = ggrid2(:)%natom0fr(gg_ini56)*expni !initial Ni56
c-- Ni
      ggrid2(:)%natom1fr(28) = ggrid2(:)%natom1fr(gg_ini56) + !unstable
     &  ggrid2(:)%natom0fr(2)                              !initial Ni (stable)
c
      end subroutine update_natomfr
