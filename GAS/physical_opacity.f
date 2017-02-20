*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2017 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
      integer :: i,j
      real*8 :: wlinv
c-- timing
      real*8 :: t0,t1,t2,t3,t4
c-- subgridded wavelength
      integer :: ngsub
      real*8 :: wlsub(grp_ng*grp_ngs+1)
c-- helper arrays
      real*8 :: grndlev(gas_ncell,ion_iionmax-1,gas_nelem)
      real*8 :: hckt(gas_ncell)
      real*8 :: hlparr(gas_ncell)
c-- ffxs
      real*8,parameter :: c1 = 4d0*pc_e**6/(3d0*pc_h*pc_me*pc_c**4)*
     &  sqrt(pc_pi2/(3*pc_me*pc_h*pc_c))
      real*8 :: gg,u,gff,help
      real*8 :: rgg,ru,dgg,du
      integer :: iu,igg
c-- bfxs
      integer :: il,ig,igs,iigs,iz,ii,ie
      logical :: dirty
      real*8 :: en,xs,wl
c-- bbxs
      real*8 :: phi,ocggrnd,wl0,dwl
      real*8 :: expfac(gas_ncell)
      real*8 :: caphelp
c-- temporary cap array in the right order
      real*8,allocatable :: cap(:,:) !(gas_ncell,grp_ng*grp_ngs)
      real*8 :: capp(gas_ncell),capr(gas_ncell)
c-- thomson scattering
      real*8,parameter :: cthomson = 8d0*pc_pi*pc_e**4/(3d0*pc_me**2
     &  *pc_c**4)
c-- warn once
      logical :: lwarn
c
c--
      ngsub = grp_ng*grp_ngs
c
c-- initialize
      allocate(cap(gas_ncell,ngsub))
      cap = 0d0
c
c-- subgridded wavelength
      igs = 1
      wlsub(igs) = grp_wl(1)
      if(grp_ngs==1) then
       wlsub = grp_wl
      else
       help = 1d0/grp_ngs
       do ig=1,grp_ng
        do i=1,grp_ngs
         igs = igs + 1
         wlsub(igs) = (1d0 - i*help)*grp_wl(ig) + i*help*grp_wl(ig+1)
        enddo
       enddo
      endif
c
c-- ion_grndlev helper array
      hckt = pc_h*pc_c/(pc_kb*gas_temp)
c
c-- thomson scattering
      if(in_nothmson) then
       gas_sig = 0d0
      else
       gas_sig = cthomson*gas_nelec*gas_natom/gas_vol
       where(gas_mass<=0d0) gas_sig = 0d0
      endif
c
      t0 = t_time()
c
c-- bound-bound
      if(.not. in_nobbopac) then
       do iz=1,gas_nelem!{{{
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &     grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)*
     &     ion_grndlev(iz,i)%ginv(ii)
        enddo !i
       enddo !iz
c
c$omp parallel
c$omp& private(wl0,iz,ii,wl,wlinv,dwl,phi,ocggrnd,expfac,caphelp,igs,
c$omp&   iigs,dirty)
c$omp& shared(ngsub,gas_mass,grndlev,hckt,wlsub,cap)
       iigs = 0
       dirty = .true.
       phi = 0d0
c$omp do schedule(static)
       do il=1,bb_nline
        wl0 = bb_xs(il)%wl0*pc_ang  !in cm
        iz = bb_xs(il)%iz
        ii = bb_xs(il)%ii
c-- igs pointer
        do igs=iigs,ngsub
         if(wlsub(igs+1)>wl0) exit
         dirty = .true.
        enddo !igs
        iigs = igs
c-- line in group
        if(igs<1) cycle
        if(igs>ngsub) cycle !can't exit in omp
c
c-- update
        if(dirty) then
         dirty = .false.
         wl = .5*(wlsub(igs) + wlsub(igs+1))
         wlinv = 1d0/wl
         dwl = wlsub(igs+1) - wlsub(igs)  !in cm
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
         if(gas_mass(i)<=0d0) cycle
         ocggrnd = grndlev(i,ii,iz)
         if(ocggrnd<=0d0) cycle
*        expfac = 1d0 - exp(-hckt(i)/wl0)  !exact expfac
         caphelp = phi*bb_xs(il)%gxs*ocggrnd*
     &     exp(-bb_xs(il)%chilw*hckt(i))*expfac(i)
!        if(caphelp==0.) write(6,*) 'cap0',cap(i,igs),phi,
!    &     bb_xs(il)%gxs,ocggrnd,exp(-bb_xs(il)%chilw*hckt(i)),expfac
         cap(i,igs) = cap(i,igs) + caphelp
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
         if(gas_mass(i)<=0d0) cycle
         forall(ii=1:min(iz,ion_el(iz)%ni - 1))
     &    grndlev(i,ii,iz) = ion_grndlev(iz,i)%oc(ii)
        enddo !i
       enddo !iz
c
c$omp parallel do
c$omp& schedule(static)
c$omp& private(wl,en,ie,xs)
c$omp& shared(ngsub,gas_mass,grndlev,wlsub,cap)
       do igs=1,ngsub
        wl = wlsub(igs)  !in cm
        en = pc_h*pc_c/(pc_ev*wl) !photon energy in eV
        do iz=1,gas_nelem
         do ii=1,min(iz,ion_el(iz)%ni - 1) !last stage is bare nucleus
          ie = iz - ii + 1
          xs = bfxs(iz,ie,en)
          if(xs==0d0) cycle
          forall(i=1:gas_ncell,.not.gas_mass(i)<=0d0)
     &      cap(i,igs) = cap(i,igs) +
     &      xs*pc_mbarn*grndlev(i,ii,iz)
         enddo !ie
        enddo !iz
!       write(6,*) 'wl done:',igs !DEBUG
!       write(6,*) cap(:,igs) !DEBUG
       enddo !igs
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
c$omp& private(wl,wlinv,u,ru,du,iu,help,gg,rgg,dgg,igg,gff)
c$omp& shared(ngsub,lwarn,hckt,hlparr,gas_mass,wlsub,cap)
       do igs=1,ngsub
        wl = wlsub(igs)  !in cm
        wlinv = 1d0/wl  !in cm
c-- gcell loop
        do i=1,gas_ncell
         if(gas_mass(i)<=0d0) cycle
         u = hckt(i)*wlinv
         ru = 10d0*(log10(u) + 4d0) + 1d0
         iu = floor(ru)
         iu = max(iu,1)
         iu = min(iu,ff_nu-1)
         du = ru - iu
c
c-- element loop
         help = c1*sqrt(hckt(i))*(1d0 - exp(-u))*wl**3*hlparr(i)
         do iz=1,gas_nelem
          gg = iz**2*pc_rydberg*hckt(i)
          rgg = 5d0*(log10(gg) + 4d0) + 1d0
          igg = floor(rgg)
          igg = max(igg,1)
          igg = min(igg,ff_ngg-1)
          dgg = rgg - igg
c
c-- bilinear inter-extrapolation
          gff = (1d0-du)*(1d0-dgg)*ff_gff(iu,igg) +
     &          du*(1d0-dgg)*ff_gff(iu+1,igg) +
     &          (1d0-du)*dgg*ff_gff(iu,igg+1) +
     &          du*dgg*ff_gff(iu+1,igg+1)
          if(rgg>ff_ngg) gff = max(gff,1d0)
c-- cross section
          cap(i,igs) = cap(i,igs) +
     &      help*gff*iz**2*gas_natom1fr(iz,i)
         enddo !iz
        enddo !i
       enddo !igs
c$omp end parallel do
c!}}}
      endif !in_noffopac
c
      t3 = t_time()
c
c-- collapse subgridding
      if(grp_ngs/=1) then
       help = 1d0/grp_ngs  !assume evenly spaced subgroups
       do ig=1,grp_ng
        i = (ig-1)*grp_ngs + 1
        j = i + grp_ngs - 1
c-- first read
        capp = sum(cap(:,i:j), dim=2, mask=cap(:,i:j)>0d0)*help
        capr = sum(1d0/cap(:,i:j), dim=2, mask=cap(:,i:j)>0d0)*help
        where(capr>0d0) capr = 1d0/capr
c-- then overwrite
        cap(:,ig) = (1d0-in_opacmixrossel)*capp + in_opacmixrossel*capr
       enddo !ig
      endif
c
c-- save
      gas_cap = transpose(sngl(cap(:,:grp_ng)))
c
c-- sanity check
      j = 0
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       do ig=1,grp_ng
        if(gas_cap(ig,i)==0.) j = ior(j,1)
        if(gas_cap(ig,i)<0.) j = ior(j,2)
        if(gas_cap(ig,i)/=gas_cap(ig,i)) j = ior(j,4)
        if(gas_cap(ig,i)>huge(gas_cap)) j = ior(j,8)
       enddo !ig
      enddo !i
      if(iand(j,1)/=0) write(0,*) 'opacity_calc: some cap==0'
      if(iand(j,2)/=0) write(0,*) 'opacity_calc: some cap<0'
      if(iand(j,4)/=0) write(0,*) 'opacity_calc: some cap==NaN'
      if(iand(j,8)/=0) write(0,*) 'opacity_calc: some cap==inf'
c-- refine sanity check
      if(j>0) then
       j = 0
       do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        do ig=1,grp_ng
         if(gas_cap(ig,i)/=0.) j = ior(j,1)
         if(gas_cap(ig,i)>=0.) j = ior(j,2)
         if(gas_cap(ig,i)==gas_cap(ig,i)) j = ior(j,4)
         if(gas_cap(ig,i)<=huge(gas_cap)) j = ior(j,8)
        enddo !ig
       enddo !i
       if(iand(j,1)==0) write(0,*) 'opacity_calc: all cap==0'
       if(iand(j,2)==0) write(0,*) 'opacity_calc: all cap<0'
       if(iand(j,4)==0) write(0,*) 'opacity_calc: all cap==NaN'
       if(iand(j,8)==0) write(0,*) 'opacity_calc: all cap==inf'
      endif
c
      deallocate(cap)
c
      t4 = t_time()
c-- register timing
      call timereg(t_opac,t4-t0)
      call timereg(t_bb,t1-t0)
      call timereg(t_bf,t2-t1)
      call timereg(t_ff,t3-t2)
c
      end subroutine physical_opacity
c vim: fdm=marker
