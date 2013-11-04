      subroutine physical_opacity
c     ---------------------------
c$    use omp_lib
      use physconstmod
      use inputparmod
      use ffxsmod
      use bfxsmod, only:bfxs
      use bbxsmod, only:bb_xs,bb_nline
      use ionsmod
      use gasgridmod
      use miscmod
      use timingmod
      implicit none
************************************************************************
* compute bound-free and bound-bound opacity.
************************************************************************
      integer :: icg
      real*8 :: wlinv
c-- timing
      real :: t0,t1
c-- helper arrays
      real*8 :: grndlev(gas_nr,ion_iionmax-1,gas_nelem)
      real*8 :: grndlev2(gas_nr,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(gas_nr)
      real*8 :: hlparr(gas_nr)
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: yend,dydx,dy !extrapolation
      integer :: iu,igg
      real*8 :: cap8
c-- bfxs
      integer :: ig,iz,ii,ie
      real*8 :: en,xs,wl,wll,wlr
      integer :: ilines,ilinee
c-- bbxs
      integer :: i,iwl
      real*8 :: phi,ocggrnd,expfac,wl0,dwl
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8 :: cap(gas_nr,in_ngs)
c-- special functions
      integer :: binsrch
      real*8 :: specint, x1, x2
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- planck opacity addition condition
      logical :: planckcheck
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*gas_temp)
c
c-- thomson scattering
      if(.not.in_nothmson) then
         gas_sig = cthomson*gas_vals2(:)%nelec*
     &        gas_vals2(:)%natom/gas_vals2(:)%volcrp
      endif
c
c-- ground level occupation number
      do iz=1,gas_nelem
       forall(icg=1:gas_nr,ii=1:min(iz,ion_el(iz)%ni - 1))
     &   grndlev(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)/
     &   ion_grndlev(iz,icg)%g(ii)
      enddo !iz
      do iz=1,gas_nelem
       forall(icg=1:gas_nr,ii=1:min(iz,ion_el(iz)%ni - 1))
     &   grndlev2(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)
      enddo !iz
c
c-- find the start point: set end before first line that falls into a group
      wlr = gas_wl(1)  !in cm
      ilines = 0
      do ilinee=ilines,bb_nline-1
       wl0 = bb_xs(ilinee+1)%wl0*pc_ang  !in cm
       if(wl0 > wlr) exit
      enddo
c
c-- bb,bf,ff opacities - group by group
      do ig=1,gas_ng
c-- right edge of the group
       wlr = gas_wl(ig+1)  !in cm
c-- bb loop start end end points
       ilines = ilinee + 1  !-- prevous end point is new starting point
       do ilinee=ilines,bb_nline-1
        wl0 = bb_xs(ilinee+1)%wl0*pc_ang  !in cm
        if(wl0 > wlr) exit
       enddo
c
       call group_opacity(ig)
c
       if(any(cap==0d0)) call warn('opacity_calc','some cap==0')
c
c-- assume evenly spaced subgroup bins
       gas_cap(ig,:) = sum(cap,dim=2)/in_ngs !assume evenly spaced subgroup bins
c-- todo: calculate gas_caprosl and gas_caprosr with cell-boundary temperature values
       gas_caprosl(ig,:) = in_ngs/sum(1/cap,dim=2) !assume evenly spaced subgroup bins
       gas_caprosr(ig,:) = gas_caprosl(ig,:)
      enddo !ig
c
c-- sanity check
      if(count(gas_caprosl-gas_cap > 1e-7*gas_cap)>0) then
       do icg=1,gas_nr
        write(6,*) icg
        write(6,'(1p,20e12.4)') gas_cap(:,icg)
        write(6,'(1p,20e12.4)') gas_caprosl(:,icg)
       enddo
c
       stop 'opacity_calc:'//
     &  'gas_capros > cas_cap'
      endif
c
c-- computing Planck opacity (rev 216)
      planckcheck = (.not.in_nobbopac.or..not.in_nobfopac.or.
     & .not.in_noffopac)
      if(planckcheck) then
         gas_siggrey = 0d0
         do icg=1,gas_nr
            do ig=1,gas_ng
               x1 = pc_h*pc_c/(gas_wl(ig+1)*pc_kb*gas_temp(icg))
               x2 = pc_h*pc_c/(gas_wl(ig)*pc_kb*gas_temp(icg))
               gas_siggrey(icg)=gas_siggrey(icg)+
     &              15d0*gas_cap(ig,icg)*specint(x1,x2,3)/pc_pi**4
            enddo
         enddo
      endif
c
      contains
c
      subroutine group_opacity(ig)
c     ----------------------------!{{{
      implicit none
      integer,intent(in) :: ig
************************************************************************
* Calculate bb,bf,ff opacity for one wl group using a refined wl subgrid
************************************************************************
      integer :: igs
c
c-- left group-boundary wavelength
      wll = gas_wl(ig)  !in cm
c-- subgroup width
      dwl = (gas_wl(ig+1) - wll)/in_ngs
c
c-- reset
      cap = 0d0
c
c-- bound-bound
      if(.not. in_nobbopac) then
       call time(t0)!{{{
c
      igs = 1
c$omp parallel do
c$omp& schedule(static)
c$omp& private(iz,ii,wl0,wlinv,phi,caphelp,expfac,ocggrnd)
c$omp& firstprivate(grndlev,hckt,igs)
c$omp& shared(cap)
       do i=ilines,ilinee
        iz = bb_xs(i)%iz
        ii = bb_xs(i)%ii
        wl0 = bb_xs(i)%wl0*pc_ang  !in cm
        wlinv = 1d0/wl0  !in cm
c-- igs pointer
        do igs=igs,in_ngs-1
         if(wl0 <= wll+igs*dwl) exit
        enddo
c-- profile function
        phi = 1d0/dwl
!       write(6,*) 'phi',phi
c-- evaluate caphelp
        do icg=1,gas_nr
         if(.not.gas_vals2(icg)%opdirty) cycle !opacities are still valid
         ocggrnd = grndlev(icg,ii,iz)
c-- oc high enough to be significant?
*        if(ocggrnd<=1d-30) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         if(ocggrnd<=0d0) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         expfac = 1d0 - exp(-hckt(icg)*wlinv)
         caphelp = phi*bb_xs(i)%gxs*ocggrnd * wl0**2/pc_c *
     &     exp(-bb_xs(i)%chilw*hckt(icg))*expfac
!        if(caphelp==0.) write(6,*) 'cap0',cap(icg,igs),phi,
!    &     bb_xs(i)%gxs,ocggrnd,exp(-bb_xs(i)%chilw*hckt(icg)),expfac
         if(caphelp==0.) cycle
         cap(icg,igs) = cap(icg,igs) + caphelp
        enddo !icg
c-- vectorized alternative is slower
cslow   where(gas_vals2(:)%opdirty .and. grndlev(:,ii,iz)>1d-30)
cslow    cap(:,igs) = cap(:,igs) +
cslow&     phi*bb_xs(i)%gxs*grndlev(:,ii,iz)*
cslow&     exp(-bb_xs(i)%chilw*hckt(:))*(1d0 - exp(-wlinv*hckt(:)))
cslow   endwhere
       enddo !i
c$omp end parallel do
c
       call time(t1)
       call timereg(t_bb, t1-t0)!}}}
      endif !in_nobbopac
c
c
c-- bound-free
      if(.not. in_nobfopac) then
       call time(t0)!{{{
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& firstprivate(grndlev2)
c$omp& shared(cap)
       do igs=1,in_ngs
        wl = wll + (igs-.5d0)*dwl !-- subgroup bin center value
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(icg=1:gas_nr)
*         forall(icg=1:gas_nr,gas_vals2(icg)%opdirty)
     &      cap(icg,igs) = cap(icg,igs) +
     &      xs*pc_mbarn*grndlev2(icg,ii,iz)
         enddo !ie
        enddo !iz
!       write(6,*) 'wl done:',igs !DEBUG
!       write(6,*) cap(:,igs) !DEBUG
       enddo !igs
c$omp end parallel do
c
       call time(t1)
       call timereg(t_bf, t1-t0)!}}}
      endif !in_nobfopac
c
c
c-- free-free
      if(.not. in_noffopac) then
       call time(t0)!{{{
c
c-- simple variant: nearest data grid point
       hlparr = (gas_vals2%natom/gas_vals2%vol)**2*gas_vals2%nelec
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,iu,help,cap8,gg,igg,gff,yend,dydx,dy)
c$omp& firstprivate(hckt,hlparr)
c$omp& shared(cap)
       do igs=1,in_ngs
        wl = wll + (igs-.5d0)*dwl !-- subgroup bin center value
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        do icg=1,gas_nr
         u = hckt(icg)*wlinv
         iu = nint(10d0*(log10(u) + 4d0)) + 1
c
         help = c1*sqrt(hckt(icg))*(1d0 - exp(-u))*wl**3*hlparr(icg)
         if(iu<1 .or. iu>ff_nu) then
          call warn('opacity_calc','ff: iu out of data limit')
          iu = min(iu,ff_nu)
          iu = max(iu,1)
         endif
c-- element loop
         cap8 = 0d0
         do iz=1,gas_nelem
          gg = iz**2*pc_rydberg*hckt(icg)
          igg = nint(5d0*(log10(gg) + 4d0)) + 1
c-- gff is approximately constant in the low igg data-limit, do trivial extrapolation:
          igg = max(igg,1)
          if(igg<=ff_ngg) then
           gff = ff_gff(iu,igg)
          else
c-- extrapolate
           yend = ff_gff(iu,ff_ngg)
           dydx = .5d0*(yend - ff_gff(iu,ff_ngg-2))
           dy = dydx*(igg - ff_ngg)
           if(abs(dy)>abs(yend - 1d0) .or. !don't cross asymptotic value
     &       sign(1d0,dy)==sign(1d0,yend - 1d0)) then !wrong slope
c-- asymptotic value
            gff = 1d0
           else
            gff = yend + dydx*(igg - ff_ngg)
           endif
          endif
c-- cross section
          cap8 = cap8 + help*gff*iz**2*gas_vals2(icg)%natom1fr(iz)
         enddo !iz
         cap(icg,igs) = cap(icg,igs) + cap8
        enddo !icg
       enddo !igs
c$omp end parallel do
c
       call time(t1)
       call timereg(t_ff, t1-t0)!}}}
      endif !in_noffopac!}}}
      end subroutine group_opacity
c
      end subroutine physical_opacity
