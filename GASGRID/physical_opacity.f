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
      integer :: i,m
      real*8 :: wlinv
c-- timing
      real*8 :: t0,t1,t2,t3,t4
c-- helper arrays
      real*8 :: grndlev(dd_ncell,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(dd_ncell)
      real*8 :: hlparr(dd_ncell)
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: yend,dydx,dy !extrapolation
      integer :: iu,igg
      real*8 :: cap1(dd_ncell)
c-- bfxs
      integer :: il,ig,iz,ii,ie
      real*8 :: en,xs,wl
c-- bbxs
      real*8 :: phi,ocggrnd,expfac,wl0,dwl
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8 :: cap(dd_ncell,gas_ng)
!-- special functions
!     integer :: binsrch
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- warn once
      logical :: lwarn
c
c-- initialize
      dd_sig = 0d0
      cap = 0d0
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*dd_temp)
c
c-- thomson scattering
      if(.not. in_nothmson) then
       dd_sig = cthomson*dd_nelec*dd_natom/dd_vol
      endif
c
      call time(t0)
c
c-- bound-bound
      if(.not. in_nobbopac) then
       do iz=1,gas_nelem!{{{
        do i=1,dd_ncell
        forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)/
     &    ion_grndlev(iz,i)%g(ii)
        enddo !i
       enddo !iz
c
       ig = 0
c$omp parallel do
c$omp& schedule(static)
c$omp& private(iz,ii,wl0,dwl,wlinv,phi,caphelp,expfac,ocggrnd,cap1)
c$omp& firstprivate(ig)
c$omp& shared(grndlev,hckt,cap)
       do il=1,bb_nline
        wl0 = bb_xs(il)%wl0*pc_ang  !in cm
c-- ig pointer
!       ig = binsrch(wl0,gas_wl,gas_ng+1,in_ng)  !todo: thread safe?
        do ig=ig,gas_ng
         if(gas_wl(ig+1)>wl0) exit
        enddo !ig
c-- line in group
        if(ig<1) cycle
        if(ig>gas_ng) cycle !can't exit in omp
        dwl = gas_wl(ig+1) - gas_wl(ig)  !in cm
c
        iz = bb_xs(il)%iz
        ii = bb_xs(il)%ii
        wlinv = 1d0/wl0  !in cm
c-- profile function
!old    phi = gas_ng*wlhelp*wl0/pc_c !line profile
        phi = wl0**2/(dwl*pc_c)
c-- evaluate caphelp
        do i=1,dd_ncell
         ocggrnd = grndlev(i,ii,iz)
c-- oc high enough to be significant?
*        if(ocggrnd<=1d-30) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         if(ocggrnd<=0d0) then
          caphelp = 0d0
         else
          expfac = 1d0 - exp(-hckt(i)*wlinv)
          caphelp = phi*bb_xs(il)%gxs*ocggrnd*
     &      exp(-bb_xs(il)%chilw*hckt(i))*expfac
!         if(caphelp==0.) write(6,*) 'cap0',cap(i,ig),phi,
!    &      bb_xs(il)%gxs,ocggrnd,exp(-bb_xs(il)%chilw*hckt(i)),expfac
         endif
         cap1(i) = caphelp
        enddo !i
        cap(:,ig) = cap(:,ig) + cap1
       enddo !il
c$omp end parallel do!}}}
      endif !in_nobbopac
c
      call time(t1)
c
c-- bound-free
      if(.not. in_nobfopac) then
c!{{{
       do iz=1,gas_nelem
        do i=1,dd_ncell
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)
        enddo !i
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& firstprivate(grndlev)
c$omp& shared(cap)
       do ig=1,gas_ng
        wl = gas_wl(ig)  !in cm
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(i=1:dd_ncell)
     &      cap(i,ig) = cap(i,ig) +
     &      xs*pc_mbarn*grndlev(i,ii,iz)
         enddo !ie
        enddo !iz
!       write(6,*) 'wl done:',ig !DEBUG
!       write(6,*) cap(:,ig) !DEBUG
       enddo !ig
c$omp end parallel do
c!}}}
      endif !in_nobfopac
c
      call time(t2)
c
c-- free-free
      if(.not. in_noffopac) then
c!{{{
c-- simple variant: nearest data grid point
       hlparr = (dd_natom/dd_vol)**2*dd_nelec
       lwarn = .true.
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,iu,help,cap1,gg,igg,gff,yend,dydx,dy)
c$omp& firstprivate(lwarn)
c$omp& shared(hckt,hlparr,cap)
       do ig=1,gas_ng
        wl = gas_wl(ig)  !in cm
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        cap1 = 0d0
        do i=1,dd_ncell
         u = hckt(i)*wlinv
         iu = nint(10d0*(log10(u) + 4d0)) + 1
c
         help = c1*sqrt(hckt(i))*(1d0 - exp(-u))*wl**3*hlparr(i)
         if(iu<1 .or. iu>ff_nu) then
          if(lwarn) then
           lwarn = .false.
           call warn('opacity_calc','ff: iu out of data limit')
          endif
          iu = min(iu,ff_nu)
          iu = max(iu,1)
         endif
c-- element loop
         do iz=1,gas_nelem
          gg = iz**2*pc_rydberg*hckt(i)
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
          cap1(i) = cap1(i) +
     &      help*gff*iz**2*dd_natom1fr(iz,i)
         enddo !iz
        enddo !i
        cap(:,ig) = cap(:,ig) + cap1
       enddo !ig
c$omp end parallel do
c!}}}
      endif !in_noffopac
c
      call time(t3)
c
      dd_cap = transpose(cap)
c
c-- sanity check
      m = 0
      do i=1,dd_ncell
       do ig=1,gas_ng
        if(dd_cap(ig,i)<=0d0) m = ior(m,1)
        if(dd_cap(ig,i)/=dd_cap(ig,i)) m = ior(m,2)
        if(dd_cap(ig,i)>huge(help)) m = ior(m,4)
       enddo !ig
      enddo !i
      if(m/=iand(m,1)) call warn('opacity_calc','some cap<=0')
      if(m/=iand(m,2)) call warn('opacity_calc','some cap==NaN')
      if(m/=iand(m,4)) call warn('opacity_calc','some cap==inf')
c
      call time(t4)
c-- register timing
      call timereg(t_opac,t4-t0)
      call timereg(t_bb,t1-t0)
      call timereg(t_bf,t2-t1)
      call timereg(t_ff,t3-t2)
c
      end subroutine physical_opacity
