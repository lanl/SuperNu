* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
      subroutine tabular_opacity(lemiss)
c     --------------------------
      use tbxsmod
      use miscmod, only:binsrch
      use elemdatamod, only:elem_data
      use gasmod
      use groupmod
      use physconstmod
      implicit none
      logical, intent(in) :: lemiss
************************************************************************
* Interpolate Fontes opacity table
************************************************************************
      integer :: itemp,irho,i,iz,ig,j,l
      real*8 :: massfr, rhopart, temph, rhoh
      real*8 :: phi1t,phi2t,phi1r,phi2r
      real*8 :: sig1,sig2,sig3,sig4
      real*8,allocatable :: sig(:)
      real*4,allocatable :: cap1(:),cap2(:),cap3(:),cap4(:)
      real*4,allocatable :: cap(:,:)
c
c-- initialize
      allocate(sig(gas_ncell))
      allocate(cap1(grp_ng))
      allocate(cap2(grp_ng))
      allocate(cap3(grp_ng))
      allocate(cap4(grp_ng))
      allocate(cap(grp_ng,gas_ncell))
      sig = 0d0
      cap = 0.
c
c-- calculate opacities
      do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
c-- table abundance for cell
        do l=1,tb_nelem
          iz=tb_ielem(l)
          if(gas_natom1fr(iz,i)<=0d0) cycle
c-- mass fraction helper
          massfr=gas_natom1fr(iz,i)*gas_natom(i) *
     &       elem_data(iz)%m*pc_amu/gas_mass(i)
c-- partial density for table interpolation
          rhopart=massfr*gas_rho(i)
c-- search 1d tb arrays for temp-rho point
          irho=binsrch(rhopart,tb_rho,tb_nrho,.false.)
          itemp=binsrch(gas_temp(i),tb_temp,tb_ntemp,.false.)
c-- do not allow for extrapolation
          rhoh=max(rhopart,tb_rho(irho))
          rhoh=min(rhoh,tb_rho(irho+1))
          temph=max(gas_temp(i),tb_temp(itemp))
          temph=min(temph,tb_temp(itemp+1))
c-- 2d bilinear interpolation in temp-rho
c-- temperature basis functions
          phi1t=(tb_temp(itemp+1)-temph) /
     &       (tb_temp(itemp+1)-tb_temp(itemp))
          phi2t=(temph-tb_temp(itemp)) /
     &       (tb_temp(itemp+1)-tb_temp(itemp))
c-- density basis functions
          phi1r=log(tb_rho(irho+1)/rhoh) /
     &       log(tb_rho(irho+1)/tb_rho(irho))
          phi2r=log(rhoh/tb_rho(irho)) /
     &       log(tb_rho(irho+1)/tb_rho(irho))
c-- scattering
          sig1=tb_sig(itemp,irho,l)*phi1t*phi1r
          sig2=tb_sig(itemp,irho+1,l)*phi1t*phi2r
          sig3=tb_sig(itemp+1,irho+1,l)*phi2t*phi2r
          sig4=tb_sig(itemp+1,irho,l)*phi2t*phi1r
          sig(i)=sig1+sig2+sig3+sig4
c-- absorption
          cap1=tb_cap(:,itemp,irho,l)*sngl(phi1t*phi1r)
          cap2=tb_cap(:,itemp,irho+1,l)*sngl(phi1t*phi2r)
          cap3=tb_cap(:,itemp+1,irho+1,l)*sngl(phi2t*phi2r)
          cap4=tb_cap(:,itemp+1,irho,l)*sngl(phi2t*phi1r)
          cap(:,i)=cap1+cap2+cap3+cap4
c-- add macroscopic opacity to total
          gas_sig(i)=gas_sig(i)+rhopart*sig(i)
          gas_cap(:,i)=gas_cap(:,i)+sngl(rhopart)*cap(:,i)
        enddo
      enddo
c
c-- sanity check (duplicated from physical_opacity)
c-- scattering
      j = 0
      do i=1,gas_ncell
       if(gas_mass(i)<=0d0) cycle
       if(gas_sig(i)==0d0) j = ior(j,1)
       if(gas_sig(i)<0d0) j = ior(j,2)
       if(gas_sig(i)/=gas_sig(i)) j = ior(j,4)
       if(gas_sig(i)>huge(gas_sig)) j = ior(j,8)
      enddo !i
      if(iand(j,1)/=0) write(0,*) 'opacity_calc: some sig==0'
      if(iand(j,2)/=0) write(0,*) 'opacity_calc: some sig<0'
      if(iand(j,4)/=0) write(0,*) 'opacity_calc: some sig==NaN'
      if(iand(j,8)/=0) write(0,*) 'opacity_calc: some sig==inf'
c-- refine sanity check
      if(j>0) then
       j = 0
       do i=1,gas_ncell
        if(gas_mass(i)<=0d0) cycle
        do ig=1,grp_ng
         if(gas_sig(i)/=0d0) j = ior(j,1)
         if(gas_sig(i)>=0d0) j = ior(j,2)
         if(gas_sig(i)==gas_sig(i)) j = ior(j,4)
         if(gas_sig(i)<=huge(gas_sig)) j = ior(j,8)
        enddo !ig
       enddo !i
       if(iand(j,1)==0) write(0,*) 'opacity_calc: all sig==0'
       if(iand(j,2)==0) write(0,*) 'opacity_calc: all sig<0'
       if(iand(j,4)==0) write(0,*) 'opacity_calc: all sig==NaN'
       if(iand(j,8)==0) write(0,*) 'opacity_calc: all sig==inf'
      endif
c-- absorption
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
c
c
c--interpolate emission opacity
      if (lemiss) then
c
c--   calculate emission opacities
        do i=1,gas_ncell
          if(gas_mass(i)<=0d0) cycle
c--   table abundance for cell
          do l=1,tb_nelem
            iz=tb_ielem(l)
            if(gas_natom1fr(iz,i)<=0d0) cycle
c--   mass fraction helper
            massfr=gas_natom1fr(iz,i)*gas_natom(i) *
     &           elem_data(iz)%m*pc_amu/gas_mass(i)
c--   partial density for table interpolation
            rhopart=massfr*gas_rho(i)
c--   search 1d tb arrays for temp-rho point
            irho=binsrch(rhopart,tb_rho,tb_nrho,.false.)
            itemp=binsrch(gas_temp(i),tb_temp,tb_ntemp,.false.)
c--   do not allow for extrapolation
            rhoh=max(rhopart,tb_rho(irho))
            rhoh=min(rhoh,tb_rho(irho+1))
            temph=max(gas_temp(i),tb_temp(itemp))
            temph=min(temph,tb_temp(itemp+1))
c--   2d bilinear interpolation in temp-rho
c--   temperature basis functions
            phi1t=(tb_temp(itemp+1)-temph) /
     &           (tb_temp(itemp+1)-tb_temp(itemp))
            phi2t=(temph-tb_temp(itemp)) /
     &           (tb_temp(itemp+1)-tb_temp(itemp))
c--   density basis functions
            phi1r=log(tb_rho(irho+1)/rhoh) /
     &           log(tb_rho(irho+1)/tb_rho(irho))
            phi2r=log(rhoh/tb_rho(irho)) /
     &           log(tb_rho(irho+1)/tb_rho(irho))
c--   emission
            cap1=tb_em_cap(:,itemp,irho,l)*sngl(phi1t*phi1r)
            cap2=tb_em_cap(:,itemp,irho+1,l)*sngl(phi1t*phi2r)
            cap3=tb_em_cap(:,itemp+1,irho+1,l)*sngl(phi2t*phi2r)
            cap4=tb_em_cap(:,itemp+1,irho,l)*sngl(phi2t*phi1r)
            cap(:,i)=cap1+cap2+cap3+cap4
c--   add macroscopic emission opacity to total
            gas_em_cap(:,i)=gas_em_cap(:,i)+sngl(rhopart)*cap(:,i)
          enddo
        enddo
c
c--   sanity check (duplicated from physical_opacity)
c--   emission
        j = 0
        do i=1,gas_ncell
          if(gas_mass(i)<=0d0) cycle
          do ig=1,grp_ng
            if(gas_em_cap(ig,i)==0.) j = ior(j,1)
            if(gas_em_cap(ig,i)<0.) j = ior(j,2)
            if(gas_em_cap(ig,i)/=gas_em_cap(ig,i)) j = ior(j,4)
            if(gas_em_cap(ig,i)>huge(gas_em_cap)) j = ior(j,8)
          enddo                 !ig
        enddo                   !i
        if(iand(j,1)/=0) write(0,*) 'opacity_calc: some em_cap==0'
        if(iand(j,2)/=0) write(0,*) 'opacity_calc: some em_cap<0'
        if(iand(j,4)/=0) write(0,*) 'opacity_calc: some em_cap==NaN'
        if(iand(j,8)/=0) write(0,*) 'opacity_calc: some em_cap==inf'
c--   refine sanity check
        if(j>0) then
          j = 0
          do i=1,gas_ncell
            if(gas_mass(i)<=0d0) cycle
            do ig=1,grp_ng
              if(gas_em_cap(ig,i)/=0.) j = ior(j,1)
              if(gas_em_cap(ig,i)>=0.) j = ior(j,2)
              if(gas_em_cap(ig,i)==gas_em_cap(ig,i)) j = ior(j,4)
              if(gas_em_cap(ig,i)<=huge(gas_em_cap)) j = ior(j,8)
            enddo               !ig
          enddo                 !i
          if(iand(j,1)==0) write(0,*) 'opacity_calc: all em_cap==0'
          if(iand(j,2)==0) write(0,*) 'opacity_calc: all em_cap<0'
          if(iand(j,4)==0) write(0,*) 'opacity_calc: all em_cap==NaN'
          if(iand(j,8)==0) write(0,*) 'opacity_calc: all em_cap==inf'
        endif
      endif
c
c-- remove opacity helpers from heap
      deallocate(sig)
      deallocate(cap,cap1,cap2,cap3,cap4)
c
      end subroutine tabular_opacity
