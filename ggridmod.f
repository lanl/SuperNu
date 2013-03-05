      module ggridmod
c     ---------------
      implicit none
************************************************************************
* gas grid structure
************************************************************************
      integer,parameter :: gg_nelem=30
      integer,parameter :: gg_ini56=-1, gg_ico56=-2 !positions in mass0fr and natom1fr arrays
c
c-- conversion factors and constants
      real*8 :: gg_vout      !outer boundary velocity
      real*8 :: gg_xi2beta   !converts position in rcell length units to v/c
      real*8 :: gg_cellength !converts rcell length units to cm
      real*8 :: gg_dxwin     !travel 'time' window
c
      integer :: gg_ncg=0    !number of gas_grid cells
c
      integer :: gg_npacket  !total # packets to be generated
      integer :: gg_mpacket  !# packets to be generated on the current mpi rank
c
c
c-- primary gas grid, available on all ranks
      type gas_grid
c-- energy
       real*8 :: enabs_c    !counted absorbed energy
       real*8 :: enabs_e    !estimated absorbed energy
c-- energy reservoir
       real*8 :: engdep     !energy deposited by gamma rays
c-- scattering
       real*8 :: sig        !Thomson scattering coeff
c-- gamma opacity
       real*8 :: capgam     !Thomson scattering coeff
      end type gas_grid
      type(gas_grid),allocatable :: ggrid(:)  !(gg_ncg)
c
c-- line opacity
      real*8,allocatable :: ggrid_wl(:) !(in_nwlg) wavelength grid
      real*8,allocatable :: ggrid2_dwl(:) !(in_nwlg) wavelength grid bin width
      real*4,allocatable :: ggrid_cap(:,:) !(gg_ncg,in_nwlg) Line+Cont extinction coeff
c
c
c-- secondary gas grid, available on master rank only
      type gas_grid2
       real*8 :: temp       !gcell temperature
       real*8 :: volr       !gcell volume [rout=1 units]
       real*8 :: vol        !gcell volume [cm^3]
       real*8 :: volcrp     !effective volume (of linked rgrid cells) [cm^3]
       real*8 :: mass       !gcell mass
       real*8 :: mass0fr(-2:gg_nelem) = 0d0  !initial mass fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       real*8 :: natom      !gcell # atoms
       real*8 :: natom1fr(-2:gg_nelem) = 0d0 !current natom fractions (>0:stable+unstable, -1:ni56, -2:co56, 0:container for unused elements)
       real*8 :: natom0fr(-2:2) = 0d0     !initial natom fractions (0,1,2:stable fe/co/ni, -1:ni56, -2:co56)
       real*8 :: nelec=1d0  !gcell # electrons per atom
c-- opacity invalidity flag
       logical :: opdirty=.true. !opacity needs recalculation
      end type gas_grid2
      type(gas_grid2),allocatable :: ggrid2(:) !(gg_ncg)
c
c-- temperature structure history
      real*8,allocatable :: ggrid2_temp(:,:) !(gg_ncg,tim_ntim)
c
      save
c
      contains
c
c
      subroutine ggrid_alloc(ncg_in,ntim_in)
c     --------------------------------------
      use inputparmod, only:in_nwlg,in_niwlem,in_ndim
      implicit none
      integer,intent(in) :: ncg_in,ntim_in
************************************************************************
* allocate ggrid variables
************************************************************************
      integer :: icgbyte
      character(28) :: labl
c
c--
      gg_ncg = ncg_in
      write(6,*) '# cells in ggrid          :',gg_ncg
c
c-- gcell size
      allocate(ggrid(1))
      icgbyte = sizeof(ggrid) + 4*in_nwlg + 8*3 !ggrid + ggrid_cap + enostor+enabs_e+enabs_c
      deallocate(ggrid)
c
c-- print used memory size
      if(in_ndim==1) then
       labl = 'allocate 1D sphericl ggrid:'
      elseif(in_ndim==2) then
       labl = 'allocate 2D cylindr ggrid :'
      endif
      write(6,'(1x,a,i10,"kB",i7,"MB")') labl,
     &   nint((icgbyte*gg_ncg)/1024.),nint((icgbyte*gg_ncg)/1024.**2)
c
c-- allocate
      allocate(ggrid(gg_ncg))       !primary gas grid
      allocate(ggrid_cap(gg_ncg,in_nwlg))
      allocate(ggrid_wl(in_nwlg))
      if(in_nwlg>huge(ggrid_icapbb)) then
       stop 'ggrid_icapbb type range too small for in_nwlg'
      endif
      allocate(ggrid_icapbb(0:in_niwlem,gg_ncg))
      allocate(ggrid2(gg_ncg))      !secondary gas grid
      allocate(ggrid2_dwl(in_nwlg)) !wavelength grid bin width
      allocate(ggrid2_temp(gg_ncg,ntim_in))
      end subroutine ggrid_alloc
c
      end module ggridmod
