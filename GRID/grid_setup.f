*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine grid_setup
c     ---------------------
      use gridmod
      use inputparmod
      use inputstrmod
      use physconstmod
      implicit none
************************************************************************
* Setup the grid on the computational domain
************************************************************************
      logical :: lexist
      integer :: i,j,k,l,idcell
      real*8 :: help
c-- statment functions
      integer :: ll
      real*8 :: xm,ym,zm
      xm(ll) = max(abs(grd_xarr(ll)),abs(grd_xarr(ll+1)))
      ym(ll) = max(abs(grd_yarr(ll)),abs(grd_yarr(ll+1)))
      zm(ll) = max(abs(grd_zarr(ll)),abs(grd_zarr(ll+1)))
c
c-- agnostic grid setup
      grd_xarr = str_xleft
      grd_yarr = str_yleft
      grd_zarr = str_zleft
c-- polar angles
      if(grd_igeom==1) grd_yacos = acos(grd_yarr)
c
c-- cell pointers
c-- initialize void cells
      grd_icell = grd_ivoid
c-- pointers into compressed grid
      l = 1
      idcell = 0
      loop_k: do k=1,grd_nz
       do j=1,grd_ny
       do i=1,grd_nx
        idcell = idcell + 1
        if(idcell == str_idcell(l)) then
         grd_icell(i,j,k) = l
         l = l + 1
        endif
        if(l>grd_ncell) exit loop_k
       enddo
       enddo
      enddo loop_k
      if(grd_ivoid>0) l = l + 1 !one dummy cell
      if(l/=grd_ncell+1) stop 'grid_setup: l/=grd_ncell+1'
c
c-- void-corner radius: cell criterium
      if(.not.in_voidcorners) then
        grd_rvoid = huge(help)
      else
       select case(grd_igeom)
       case(1,11)
        grd_rvoid = grd_xarr(grd_nx+1)
       case(2)
        grd_rvoid = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1))
       case(3)
        grd_rvoid = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1),
     &    grd_zarr(grd_nz+1))
       endselect
      endif
c
c-- maximum grid velocity
      select case(grd_igeom)
      case(1,11)
       grd_rout = grd_xarr(grd_nx+1)
c-- cylindrical
      case(2)
c-- box
       if(.not.in_voidcorners) then
        grd_rout = sqrt(grd_xarr(grd_nx+1)**2 +
     &    max(-grd_yarr(1),grd_yarr(grd_ny+1))**2)
c-- sphere
       else
        grd_rout = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1))**2 !squared
        do k=1,grd_nz
        do j=1,grd_ny
        do i=1,grd_nx
         if(grd_icell(i,j,k)==grd_ivoid) cycle !void
         help = grd_xarr(i+1)**2 + ym(j)**2
         grd_rout = max(grd_rout,help)
        enddo
        enddo
        enddo
        grd_rout = sqrt(grd_rout)
       endif
c-- cartesian
      case(3)
c-- box
       if(.not.in_voidcorners) then
        grd_rout = sqrt(
     &    max(-grd_xarr(1),grd_xarr(grd_nx+1))**2 +
     &    max(-grd_yarr(1),grd_yarr(grd_ny+1))**2 +
     &    max(-grd_zarr(1),grd_zarr(grd_nz+1))**2)
c-- sphere
       else
        grd_rout = min(grd_xarr(grd_nx+1),grd_yarr(grd_ny+1),
     &    grd_zarr(grd_nz+1))**2 !squared
        do k=1,grd_nz
        do j=1,grd_ny
        do i=1,grd_nx
         if(grd_icell(i,j,k)==grd_ivoid) cycle !void
         help = xm(i)**2 + ym(j)**2 + zm(k)**2
         grd_rout = max(grd_rout,help)
        enddo
        enddo
        enddo
        grd_rout = sqrt(grd_rout)
       endif
      endselect
c
c-- sanity check
      if(grd_isvelocity) then
       if(maxval(abs(grd_xarr))>pc_c) stop 'grid_setup: grd_xarr > pc_c'
       if(maxval(abs(grd_yarr))>pc_c) stop 'grid_setup: grd_yarr > pc_c'
       if(maxval(abs(grd_zarr))>pc_c) stop 'grid_setup: grd_zarr > pc_c'
      endif
c-- sanity check
      select case(grd_igeom)
      case(1,11)
       if(minval(grd_xarr)<0d0) stop 'grid_setup: grd_xarr < 0'
      case(2)
       if(minval(grd_xarr)<0d0) stop 'grid_setup: grd_xarr < 0'
      endselect
c
c-- zero amplification-factor energy to begin with
      grd_eamp = 0d0
c
c-- read preset temperature profiles
      inquire(file='input.temp',exist=lexist)
      if(lexist) call read_temp_preset
c
      end subroutine grid_setup
c vim: fdm=marker
