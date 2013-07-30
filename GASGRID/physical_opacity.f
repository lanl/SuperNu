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
      integer :: iw,iz,ii,ie
      real*8 :: en,xs,wl
c-- bbxs
      integer :: i,iwl
      real*8 :: phi,ocggrnd,expfac,wl0,dwl
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8 :: cap(gas_nr,gas_ng)
c-- special functions
      integer :: binsrch
      real*8 :: specint, x1, x2
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- planck opacity addition condition
      logical :: planckcheck
!c
!c-- constants
!old  wlhelp = 1d0/log(in_wlmax/dble(in_wlmin))
!old  wlminlg = log(dble(in_wlmin))
c
c-- reset
      cap = 0d0
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
c-- bound-bound
      if(.not. in_nobbopac) then
       call time(t0)!{{{

       do iz=1,gas_nelem
        forall(icg=1:gas_nr,ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)/
     &    ion_grndlev(iz,icg)%g(ii)
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(iz,ii,wl0,dwl,wlinv,iwl,phi,caphelp,expfac,ocggrnd)
c$omp& firstprivate(grndlev,hckt)
c$omp& shared(cap)
       do i=1,bb_nline
        iz = bb_xs(i)%iz
        ii = bb_xs(i)%ii
        wl0 = bb_xs(i)%wl0*pc_ang  !in cm
        wlinv = 1d0/wl0  !in cm
c-- iwl pointer
        iwl = binsrch(wl0,gas_wl,gas_ng+1)  !todo: thread safe?
c--
        if(iwl<1) cycle
        if(iwl>gas_ng) cycle
        dwl = gas_wl(iwl+1) - gas_wl(iwl)  !in cm
c-- profile function
!old    phi = gas_ng*wlhelp*wl0/pc_c !line profile
        phi = wl0**2/(dwl*pc_c)
!       write(6,*) 'phi',phi
c-- evaluate caphelp
        do icg=1,gas_nr
         if(.not.gas_vals2(icg)%opdirty) cycle !opacities are still valid
         ocggrnd = grndlev(icg,ii,iz)
c-- oc high enough to be significant?
*        if(ocggrnd<=1d-30) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         if(ocggrnd<=0d0) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         expfac = 1d0 - exp(-hckt(icg)*wlinv)
         caphelp = phi*bb_xs(i)%gxs*ocggrnd*
     &     exp(-bb_xs(i)%chilw*hckt(icg))*expfac
!        if(caphelp==0.) write(6,*) 'cap0',cap(icg,iwl),phi,
!    &     bb_xs(i)%gxs,ocggrnd,exp(-bb_xs(i)%chilw*hckt(icg)),expfac
         if(caphelp==0.) cycle
         cap(icg,iwl) = cap(icg,iwl) + caphelp
        enddo !icg
c-- vectorized alternative is slower
cslow   where(gas_vals2(:)%opdirty .and. grndlev(:,ii,iz)>1d-30)
cslow    cap(:,iwl) = cap(:,iwl) +
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
       do iz=1,gas_nelem
        forall(icg=1:gas_nr,ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& firstprivate(grndlev)
c$omp& shared(cap)
       do iw=1,gas_ng
        wl = gas_wl(iw)  !in cm
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(icg=1:gas_nr)
*         forall(icg=1:gas_nr,gas_vals2(icg)%opdirty)
     &      cap(icg,iw) = cap(icg,iw) +
     &      xs*pc_mbarn*grndlev(icg,ii,iz)
         enddo !ie
        enddo !iz
!       write(6,*) 'wl done:',iw !DEBUG
!       write(6,*) cap(:,iw) !DEBUG
       enddo !iw
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
       do iw=1,gas_ng
        wl = gas_wl(iw)  !in cm
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
         cap(icg,iw) = cap(icg,iw) + cap8
        enddo !icg
       enddo !iw
c$omp end parallel do
c
       call time(t1)
       call timereg(t_ff, t1-t0)!}}}
      endif !in_noffopac
c
      if(any(cap==0d0))
     & call warn('opacity_calc','some cap==0')
c
      gas_cap = gas_cap + transpose(cap)
c
c-- computing Planck opacity (rev 216)
      planckcheck = (.not.in_nobbopac.or..not.in_nobfopac.or.
     & .not.in_noffopac)
      if(planckcheck) then
         gas_siggrey = 0d0
         do icg=1,gas_nr
            do iw=1,gas_ng
               x1 = pc_h*pc_c/(gas_wl(iw+1)*pc_kb*gas_temp(icg))
               x2 = pc_h*pc_c/(gas_wl(iw)*pc_kb*gas_temp(icg))
               gas_siggrey(icg)=gas_siggrey(icg)+
     &              15d0*gas_cap(iw,icg)*specint(x1,x2,3)/pc_pi**4
            enddo
         enddo
      endif
c      hckt = pc_h*pc_c/(pc_kb*gas_temp)
      end subroutine physical_opacity
