      subroutine physical_opacity
c     ---------------------------
c$    use omp_lib
      use physconstmod
      use inputparmod
      use ffxsmod
      use bfxsmod, only:bfxs
      use bbxsmod, only:bb_xs,bb_nline
      use ionsmod
      use gasmod
      use groupmod
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
      real*8 :: grndlev(gas_ncell,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(gas_ncell)
      real*8 :: hlparr(gas_ncell)
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: yend,dydx,dy !extrapolation
      integer :: iu,igg
c-- bfxs
      integer :: il,iig,ig,iz,ii,ie
      logical :: dirty
      real*8 :: en,xs,wl
c-- bbxs
      real*8 :: phi,ocggrnd,wl0,dwl
      real*8 :: expfac(gas_ncell)
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8,allocatable :: cap(:,:)
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- warn once
      logical :: lwarn
c
c-- initialize
      allocate(cap(gas_ncell,grp_ng))
      cap = 0d0
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*gas_temp)
c
c-- thomson scattering
      if(in_nothmson) then
       gas_sig = 0d0
      else
       gas_sig = cthomson*gas_nelec*gas_natom/gas_vol
      endif
c
      t0 = t_time()
c
c-- bound-bound
      if(.not. in_nobbopac) then
       do iz=1,gas_nelem!{{{
        do i=1,gas_ncell
         if(gas_void(i)) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &     grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)*
     &     ion_grndlev(iz,i)%ginv(ii)
        enddo !i
       enddo !iz
c
c$omp parallel
c$omp& private(wl0,iz,ii,wl,wlinv,dwl,phi,ocggrnd,expfac,caphelp,ig,iig,
c$omp&   dirty)
c$omp& shared(gas_void,grndlev,hckt,cap)
       iig = 0
       dirty = .true.
       phi = 0d0
c$omp do schedule(static)
       do il=1,bb_nline
        wl0 = bb_xs(il)%wl0*pc_ang  !in cm
        iz = bb_xs(il)%iz
        ii = bb_xs(il)%ii
c-- ig pointer
        do ig=iig,grp_ng
         if(grp_wl(ig+1)>wl0) exit
         dirty = .true.
        enddo !ig
        iig = ig
c-- line in group
        if(ig<1) cycle
        if(ig>grp_ng) cycle !can't exit in omp
c
c-- update
        if(dirty) then
         dirty = .false.
         wl = .5*(grp_wl(ig) + grp_wl(ig+1))
         wlinv = 1d0/wl
         dwl = grp_wl(ig+1) - grp_wl(ig)  !in cm
c-- approximate expfac
         expfac = 1d0 - exp(-hckt*wlinv)
c-- approximate profile function
         phi = wl**2/(dwl*pc_c)
        endif
c
c-- profile function
*       phi = wl0**2/(dwl*pc_c)  !exact phi
c
c-- evaluate cap
        do i=1,gas_ncell
         if(gas_void(i)) cycle
         ocggrnd = grndlev(i,ii,iz)
         if(ocggrnd<=0d0) cycle
*        expfac = 1d0 - exp(-hckt(i)/wl0)  !exact expfac
         caphelp = phi*bb_xs(il)%gxs*ocggrnd*
     &     exp(-bb_xs(il)%chilw*hckt(i))*expfac(i)
!        if(caphelp==0.) write(6,*) 'cap0',cap(i,ig),phi,
!    &     bb_xs(il)%gxs,ocggrnd,exp(-bb_xs(il)%chilw*hckt(i)),expfac
         cap(i,ig) = cap(i,ig) + caphelp
        enddo !i
       enddo !il
c$omp end do
c$omp end parallel
c!}}}
      endif !in_nobbopac
c
      t1 = t_time()
c
c-- bound-free
      if(.not. in_nobfopac) then
c!{{{
       do iz=1,gas_nelem
        do i=1,gas_ncell
         if(gas_void(i)) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)
        enddo !i
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& shared(gas_void,grndlev,cap)
       do ig=1,grp_ng
        wl = grp_wl(ig)  !in cm
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(i=1:gas_ncell,.not.gas_void(i))
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
      t2 = t_time()
c
c-- free-free
      if(.not. in_noffopac) then
c!{{{
c-- simple variant: nearest data grid point
       hlparr = (gas_natom/gas_vol)**2*gas_nelec
       lwarn = .true.
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,iu,help,gg,igg,gff,yend,dydx,dy)
c$omp& shared(lwarn,hckt,hlparr,gas_void,cap)
       do ig=1,grp_ng
        wl = grp_wl(ig)  !in cm
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        do i=1,gas_ncell
         if(gas_void(i)) cycle
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
          cap(i,ig) = cap(i,ig) +
     &      help*gff*iz**2*gas_natom1fr(iz,i)
         enddo !iz
        enddo !i
       enddo !ig
c$omp end parallel do
c!}}}
      endif !in_noffopac
c
      t3 = t_time()
c
      gas_cap = transpose(sngl(cap))
c
c-- sanity check
      m = 0
      do i=1,gas_ncell
       if(gas_void(i)) cycle
       do ig=1,grp_ng
        if(gas_cap(ig,i)<=0.) m = ior(m,1)
        if(gas_cap(ig,i)/=gas_cap(ig,i)) m = ior(m,2)
        if(gas_cap(ig,i)>huge(gas_cap)) m = ior(m,4)
       enddo !ig
      enddo !i
      if(iand(m,1)/=0) call warn('opacity_calc','some cap<=0')
      if(iand(m,2)/=0) call warn('opacity_calc','some cap==NaN')
      if(iand(m,4)/=0) call warn('opacity_calc','some cap==inf')
c
      t4 = t_time()
c-- register timing
      call timereg(t_opac,t4-t0)
      call timereg(t_bb,t1-t0)
      call timereg(t_bf,t2-t1)
      call timereg(t_ff,t3-t2)
c
      end subroutine physical_opacity
