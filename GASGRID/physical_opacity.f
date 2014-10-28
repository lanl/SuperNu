      subroutine physical_opacity
c     ---------------------------
c$    use omp_lib
      use mpimod
      use physconstmod
      use inputparmod
      use ffxsmod
      use bfxsmod, only:bfxs
      use bbxsmod, only:bb_xs,bb_nline
      use ionsmod
      use gasgridmod, nx=>gas_nx,ny=>gas_ny,nz=>gas_nz
      use miscmod
      use timingmod
      implicit none
************************************************************************
* compute bound-free and bound-bound opacity.
************************************************************************
      integer :: i,j,k,l,n
      integer :: l0,l1,l2,i1,i2,j1,j2,k1,k2
      real*8 :: wlinv
c-- timing
      real*8 :: t0,t1,t2,t3,t4
c-- helper arrays
      real*8,allocatable :: grndlev(:,:,:) !(n,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(nx,ny,nz)
      real*8 :: hlparr(nx,ny,nz)
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: yend,dydx,dy !extrapolation
      integer :: iu,igg
      real*8 :: cap8
c-- bfxs
      integer :: il,ig,iz,ii,ie
      real*8 :: en,xs,wl
c-- bbxs
      real*8 :: phi,ocggrnd,expfac,wl0,dwl
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8,allocatable :: cap(:,:) !(n,gas_ng)
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- warn once
      logical :: lwarn
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*gas_temp)
c
c
c-- MPI task distribution
      n = nx*ny*nz!{{{
      do i=nmpi-1,impi,-1
       l = int(n/dble(i+1)) !divide remaining
       n = n - l !remaining
      enddo
      l1 = n+1 !start
      l2 = n+l !end
c
c-- translate range to index ranges
      l = l1-1
      i1 = nx
      j1 = ny
      k1 = nz
      i2 = 1
      j2 = 1
      k2 = 1
      do k=1,nz
      do j=1,ny
      do i=1,nx
       l = l + 1
       if(l<l1 .or. l>l2) cycle
       i1 = max(i1,i)
       j1 = max(j1,j)
       k1 = max(k1,k)
       i2 = min(i2,i)
       j2 = min(j2,j)
       k2 = min(k2,k)
      enddo !i
      enddo !j
      enddo !k
c-- find counter start
      l = 0
      do k=1,nz
      do j=1,ny
      do i=1,nx
       l = l + 1
       if(i==i1 .and. j==j1 .and. k==k1) l0 = l-1
      enddo !i
      enddo !j
      enddo !k!}}}
c
c
c-- allocate
      allocate(cap(l1:l2,gas_ng))
      allocate(grndlev(l1:l2,ion_iionmax-1,gas_nelem)) !(n,ion_iionmax-1,gas_nelem)
      cap = 0d0
c
c-- thomson scattering
      if(.not.in_nothmson) then
       gas_sig = cthomson*gas_vals2%nelec*gas_vals2%natom/gas_vals2%vol
      endif
c
      call time(t0)
c
c-- bound-bound
      if(.not. in_nobbopac) then
!{{{
       l = l0
       do k=k1,k2
       do j=j1,j2
       do i=i1,i1
        l = l + 1
        if(l<l1 .or. l>l2) cycle
        do iz=1,gas_nelem
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &     grndlev(l,ii,iz) = ion_grndlev(iz,i,j,k)%oc(ii)/
     &     ion_grndlev(iz,i,j,k)%g(ii)
        enddo !iz
       enddo !i
       enddo !j
       enddo !k
c
       ig = 0
c$omp parallel do
c$omp& schedule(static)
c$omp& private(iz,ii,wl0,dwl,wlinv,phi,caphelp,expfac,ocggrnd,l)
c$omp& firstprivate(ig)
c$omp& shared(grndlev,hckt,cap)
       do il=1,bb_nline
        wl0 = bb_xs(il)%wl0*pc_ang  !in cm
c-- ig pointer
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
        l = l0
        do k=k1,k2
        do j=j1,j2
        do i=i1,i1
         l = l + 1
         if(l<l1 .or. l>l2) cycle
         ocggrnd = grndlev(l,ii,iz)
c-- oc high enough to be significant?
*        if(ocggrnd<=1d-30) cycle !todo: is this _always_ low enoug? It is in the few tests I did.
         if(ocggrnd<=0d0) cycle
         expfac = 1d0 - exp(-hckt(i,j,k)*wlinv)
         caphelp = phi*bb_xs(il)%gxs*ocggrnd*
     &     exp(-bb_xs(il)%chilw*hckt(i,j,k))*expfac
!        if(caphelp==0.) write(6,*) 'cap0',cap(i,j,k,ig),phi,
!    &     bb_xs(il)%gxs,ocggrnd,exp(-bb_xs(il)%chilw*hckt(i,j,k)),expfac
         if(caphelp==0.) cycle
         cap(l,ig) = cap(l,ig) + caphelp
        enddo !i
        enddo !j
        enddo !k
       enddo !il
c$omp end parallel do!}}}
      endif !in_nobbopac
c
      call time(t1)
c
c-- bound-free
      if(.not. in_nobfopac) then
c!{{{
       l = l0
       do k=k1,k2
       do j=j1,j2
       do i=i1,i2
        l = l + 1
        if(l<l1 .or. l>l2) cycle
        do iz=1,gas_nelem
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &     grndlev(l,ii,iz) = ion_grndlev(iz,i,j,k)%oc(ii)
        enddo !iz
       enddo !i
       enddo !j
       enddo !k
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs,l)
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
          do k=k1,k2
          do j=j1,j2
          do i=i1,i2
           l = l + 1
           if(l<l1 .or. l>l2) cycle
           cap(l,ig) = cap(l,ig) + xs*pc_mbarn*grndlev(l,ii,iz)
          enddo !i
          enddo !j
          enddo !k
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
       hlparr = (gas_vals2%natom/gas_vals2%vol)**2*gas_vals2%nelec
       lwarn = .true.
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,wlinv,u,iu,help,cap8,gg,igg,gff,yend,dydx,dy,l)
c$omp& firstprivate(hckt,hlparr,lwarn)
c$omp& shared(cap)
       do ig=1,gas_ng
        wl = gas_wl(ig)  !in cm
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        l = l0
        do k=k1,k2
        do j=j1,j2
        do i=i1,i2
         l = l + 1
         if(l<l1 .or. l>l2) cycle
         u = hckt(i,j,k)*wlinv
         iu = nint(10d0*(log10(u) + 4d0)) + 1
c
         help = c1*sqrt(hckt(i,j,k))*(1d0 - exp(-u))*wl**3*hlparr(i,j,k)
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
          gg = iz**2*pc_rydberg*hckt(i,j,k)
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
          cap(l,ig) = cap(l,ig) +
     &      help*gff*iz**2*gas_vals2(i,j,k)%natom1fr(iz)
         enddo !iz
        enddo !i
        enddo !j
        enddo !k
       enddo !ig
c$omp end parallel do
c!}}}
      endif !in_noffopac
c
      call time(t3)
c
c-- insert into gas_cap
      gas_cap = 0.
      l = l0
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
       l = l + 1
       if(l<l1 .or. l>l2) cycle
       gas_cap(:,i,j,k) = sngl(cap(l,:))
      enddo !i
      enddo !j
      enddo !k
c
c-- sanity check
      m = 0
      l = l0
      do k=k1,k2
      do j=j1,j2
      do i=i1,i2
       l = l + 1
       if(l<l1 .or. l>l2) cycle
       do ig=1,gas_ng
        if(gas_cap(ig,i,j,k)<=0d0) m = ior(m,1)
        if(gas_cap(ig,i,j,k)/=gas_cap(ig,i,j,k)) m = ior(m,2)
        if(gas_cap(ig,i,j,k)>huge(help)) m = ior(m,4)
       enddo !ig
      enddo !i
      enddo !j
      enddo !k
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
