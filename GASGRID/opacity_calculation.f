      subroutine opacity_calculation
c     ------------------------------
c$    use omp_lib
      use physconstmod
      use fluxmod
      use inputparmod
      use ffxsmod
      use bfxsmod, only:bfxs
      use bbxsmod, only:bb_xs,bb_nline
      use ionsmod
      use ggridmod
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
      real*8 :: grndlev(gg_ncg,ion_iionmax-1,gg_nelem)
      real*8 :: hckt(gg_ncg)
      real*8 :: hlparr(gg_ncg)
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
      real*8 :: phi,ocggrnd,expfac,wl0
      real :: cap
c
c-- reset
      ggrid_cap = 0.
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*ggrid2%temp)
c
c
c-- bound-bound
      if(.not. in_nobbopac) then
       call time(t0)!{{{

       do iz=1,gg_nelem
        forall(icg=1:gg_ncg,ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)/
     &    ion_grndlev(iz,icg)%g(ii)
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(iz,ii,wl0,wlinv,iwl,phi,cap,expfac,ocggrnd)
c$omp& firstprivate(grndlev,hckt)
c$omp& shared(ggrid_cap)
       do i=1,bb_nline
        iz = bb_xs(i)%iz
        ii = bb_xs(i)%ii
        wl0 = bb_xs(i)%wl0 !in ang
        wlinv = 1d0/(wl0*pc_ang)
c-- iwl pointer
        iwl = int((flx_wlhelp*(in_nwlg - 1d0))*(log(dble(wl0)) - !sensitive to multiplication order!
     &    flx_wlminlg)) + 1
        if(iwl<1) cycle
        if(iwl>in_nwlg) cycle
c-- profile function
        phi = (in_nwlg-1d0)*flx_wlhelp*(wl0*pc_ang)/pc_c !line profile
!       write(*,*) 'phi',phi
c-- evaluate cap
        do icg=1,gg_ncg
         if(.not.ggrid2(icg)%opdirty) cycle !opacities are still valid
         ocggrnd = grndlev(icg,ii,iz)
c-- oc high enough to be significant?
*        if(ocggrnd<=1d-30) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         if(ocggrnd<=0d0) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         expfac = 1d0 - exp(-hckt(icg)*wlinv)
         cap = sngl(phi*bb_xs(i)%gxs*ocggrnd*
     &     exp(-bb_xs(i)%chilw*hckt(icg))*expfac)
!        if(cap==0.) write(6,*) 'cap0',ggrid_cap(icg,iwl),phi,
!    &     bb_xs(i)%gxs,ocggrnd,exp(-bb_xs(i)%chilw*hckt(icg)),expfac
         if(cap==0.) cycle
         ggrid_cap(icg,iwl) = ggrid_cap(icg,iwl) + cap
        enddo !icg
c-- vectorized alternative is slower
cslow   where(ggrid2(:)%opdirty .and. grndlev(:,ii,iz)>1d-30)
cslow    ggrid_cap(:,iwl) = ggrid_cap(:,iwl) +
cslow&     sngl(phi*bb_xs(i)%gxs*grndlev(:,ii,iz)*
cslow&     exp(-bb_xs(i)%chilw*hckt(:))*(1d0 - exp(-wlinv*hckt(:))))
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
       do iz=1,gg_nelem
        forall(icg=1:gg_ncg,ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(icg,ii,iz) = ion_grndlev(iz,icg)%oc(ii)
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& firstprivate(grndlev)
c$omp& shared(ggrid_cap)
       do iw=1,in_nwlg
        wl = ggrid_wl(iw)
        en = pc_h*pc_c/(pc_ev*pc_ang*wl) !photon energy in eV
        do iz=1,gg_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(icg=1:gg_ncg)
*         forall(icg=1:gg_ncg,ggrid2(icg)%opdirty)
     &      ggrid_cap(icg,iw) = ggrid_cap(icg,iw) +
     &      sngl(xs*pc_mbarn*grndlev(icg,ii,iz))
         enddo !ie
        enddo !iz
!       write(*,*) 'wl done:',iw !DEBUG
!       write(*,*) ggrid_cap(:,iw) !DEBUG
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
       hlparr = (ggrid2%natom/ggrid2%vol)**2*ggrid2%nelec
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,iu,help,cap8,gg,igg,gff,yend,dydx,dy)
c$omp& firstprivate(hckt,hlparr)
c$omp& shared(ggrid_cap)
       do iw=1,in_nwlg
        wl = pc_ang*ggrid_wl(iw)
        wlinv = 1d0/wl
c-- gcell loop
        do icg=1,gg_ncg
         u = hckt(icg)*wlinv
         iu = nint(10d0*(log10(u) + 4d0)) + 1
c
         help = c1*sqrt(hckt(icg))*(1d0 - exp(-u))*wl**3*hlparr(icg)
         if(iu<1 .or. iu>ff_nu) stop 'opacity_calc: ff: iu wrong'
c-- element loop
         cap8 = 0d0
         do iz=1,gg_nelem
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
          cap8 = cap8 + help*gff*iz**2*ggrid2(icg)%natom1fr(iz)
         enddo !iz
         ggrid_cap(icg,iw) = ggrid_cap(icg,iw) + sngl(cap8)
        enddo !icg
       enddo !iw
c$omp end parallel do
c
       call time(t1)
       call timereg(t_ff, t1-t0)!}}}
      endif !in_noffopac
c
      if(any(ggrid_cap==0.))
     & call warn('opacity_calc','some ggrid_cap==0')
c
c-- convert to opacity per rcell
      ggrid_cap = ggrid_cap*gg_cellength !gg_cellength converts cm^-1 to 1/rcell
c
      end subroutine opacity_calculation
