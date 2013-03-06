      subroutine ggrid_setup(ncg_in,ntim_in)
c     --------------------------------------
      use physconstmod
      use inputparmod
      use gasgridmod
      use miscmod, only:warn
      implicit none
      integer,intent(in) :: ncg_in,ntim_in
************************************************************************
* Initialize the gas grid, the part that is constant with time and
* temperature. The part that changes is done in gas_grid_update.
************************************************************************
      integer :: i
      real*8 :: help
c
c--
      write(6,*)
      write(6,*) 'setup gas grid:'     
      write(6,*) '==========================='
c
c-- allocate ggrid
      call ggrid_alloc(ncg_in,ntim_in)
c
c-- init temp structure
      ggrid2_temp(:,1) = 1d4  !initial guess, may be overwritten by read_temp_str
      ggrid2_temp(:,2:) = 0d0 !either overwritten by read_temp_str or copy from previous time step
c
c-- convert mass fractions to # atoms
      call massfr2natomfr
!c
!c-- output
!      write(6,*) 'mass fractions'
!      write(6,'(1p,33i12)') (i,i=-2,30)
!      write(6,'(1p,33e12.4)') (ggrid2(i)%mass0fr,i=1,gg_ncg)
!      write(6,*) 'number fractions'
!      write(6,'(1p,33i12)') (i,i=-2,30)
!      write(6,'(1p,33e12.4)') (ggrid2(i)%natom1fr,i=1,gg_ncg)
c
c-- gas wavelength grid
      forall(i=1:in_nwlg) ggrid_wl(i) =
     &  in_wlmin*(in_wlmax/dble(in_wlmin))**((i-1d0)/(in_nwlg-1d0))
      ggrid2_dwl = pc_ang*ggrid_wl*log(in_wlmax/dble(in_wlmin)) /
     &  (in_nwlg-1d0) !wl grid bin width
c-- sanity test
      help = sum(ggrid2_dwl)/pc_ang
      if(abs(help/(in_wlmax-in_wlmin) - 1d0) .gt. 1d-3) then
       call warn('ggrid_setup','ggrid_dwl not accurate')
       write(6,*) help,in_wlmax-in_wlmin
      endif
c
      end subroutine ggrid_setup
c
c
c
      subroutine massfr2natomfr
c     -------------------------
      use physconstmod
      use elemdatamod, only:elem_data
      use ggridmod
      implicit none
************************************************************************
* convert mass fractions to natom fractions, and mass to natom.
************************************************************************
      integer :: i,j
      real*8 :: help
c
      do i=1,gg_ncg
c
c-- renormalize (the container fraction (unused elements) is taken out)
       ggrid2(i)%mass0fr(:) = ggrid2(i)%mass0fr(:)/
     &   sum(ggrid2(i)%mass0fr(1:))
c
c-- partial mass
       ggrid2(i)%natom1fr = ggrid2(i)%mass0fr*ggrid2(i)%mass
c-- only stable nickel and cobalt
       ggrid2(i)%natom1fr(28) = ggrid2(i)%natom1fr(28) -
     &   ggrid2(i)%natom1fr(gg_ini56)
       ggrid2(i)%natom1fr(27) = ggrid2(i)%natom1fr(27) -
     &   ggrid2(i)%natom1fr(gg_ico56)
c
c-- convert to natoms
       do j=1,gg_nelem
        ggrid2(i)%natom1fr(j) = ggrid2(i)%natom1fr(j)/
     &    (elem_data(j)%m*pc_amu) 
       enddo !j
c-- special care for ni56 and co56
       help = elem_data(26)%m*pc_amu
!      help = elem_data(28)%m*pc_amu !phoenix compatible
       ggrid2(i)%natom1fr(gg_ini56) = ggrid2(i)%natom1fr(gg_ini56)/help
!      help = elem_data(27)%m*pc_amu !phoenix compatible
       ggrid2(i)%natom1fr(gg_ico56) = ggrid2(i)%natom1fr(gg_ico56)/help
c-- store initial fe/co/ni
       ggrid2(i)%natom0fr(-2:-1) = ggrid2(i)%natom1fr(-2:-1) !unstable
       ggrid2(i)%natom0fr(0:2) = ggrid2(i)%natom1fr(26:28) !stable
c-- add unstable to stable again
       ggrid2(i)%natom1fr(28) = ggrid2(i)%natom1fr(28) +
     &   ggrid2(i)%natom1fr(gg_ini56)
       ggrid2(i)%natom1fr(27) = ggrid2(i)%natom1fr(27) +
     &   ggrid2(i)%natom1fr(gg_ico56)
c
c-- total natom
       ggrid2(i)%natom = sum(ggrid2(i)%natom1fr(1:))
c
c-- convert natoms to natom fractions
       ggrid2(i)%natom1fr = ggrid2(i)%natom1fr/ggrid2(i)%natom
       ggrid2(i)%natom0fr = ggrid2(i)%natom0fr/ggrid2(i)%natom
c
      enddo !i
c
      end subroutine massfr2natomfr
