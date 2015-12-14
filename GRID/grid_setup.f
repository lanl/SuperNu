      subroutine grid_setup
c     ---------------------
      use gridmod
      use inputstrmod
      use physconstmod
      implicit none
************************************************************************
* Setup the grid on the computational domain
************************************************************************
      logical :: lexist
      integer :: i,j,k,l,idcell
c
c-- agnostic grid setup
      grd_xarr = str_xleft
      grd_yarr = str_yleft
      grd_zarr = str_zleft
c-- polar angles
      if(grd_igeom==1) grd_yacos = acos(grd_yarr)
c
c-- maximum grid velocity
      select case(grd_igeom)
      case(1,11)
       grd_voc = grd_xarr(grd_nx+1)
      case(2)
       grd_voc = sqrt(grd_xarr(grd_nx+1)**2 +
     &   max(-grd_yarr(1),grd_yarr(grd_ny+1))**2)
      case(3)
       grd_voc = sqrt(
     &   max(-grd_xarr(1),grd_xarr(grd_nx+1))**2 +
     &   max(-grd_yarr(1),grd_yarr(grd_ny+1))**2 +
     &   max(-grd_zarr(1),grd_zarr(grd_nz+1))**2)
      endselect
      grd_voc = grd_voc/pc_c
c
c-- cell pointers
c-- if dummy cell exists initialize void cells
      if(grd_lvoid) grd_icell = grd_ncell
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
      if(grd_lvoid) l = l + 1 !one dummy cell
      if(l/=grd_ncell+1) stop 'grid_setup: l/=grd_ncell+1'
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
