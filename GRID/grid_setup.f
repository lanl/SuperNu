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
      if(grd_igeom==1) then
       grd_yacos = acos(grd_yarr)
      else
       grd_yacos = 0d0
      endif
c
c-- cell pointers
c-- if padding exists initialize void cells
      if(grd_ncp/=grd_nc) grd_icell = grd_ncp
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
        if(l>grd_nc) exit loop_k
       enddo
       enddo
      enddo loop_k
      if(l/=grd_nc+1) stop 'grid_setup: l/=grd_nc+1'
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
