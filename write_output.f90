subroutine write_output(nmpi)

  use inputparmod
  use timingmod
  use timestepmod
  use gridmod
  use totalsmod
  use fluxmod
  implicit none
  integer,intent(in) :: nmpi

  integer :: i,j,k
  integer :: ncpr !number of cells per rank
  integer :: reclenf, recleng
  logical,save :: lfirst=.true.
  character(16), save :: pos='rewind', fstat='replace'
!
  integer,allocatable :: iarr(:,:)
  real*8,allocatable :: arr(:,:)
!
!-- helper arrays
  ncpr = grd_ncp/nmpi
  allocate(iarr(ncpr,nmpi),arr(ncpr,nmpi))
!
  reclenf = (flx_ng+1)*12
  recleng = (ncpr+1)*12

!
!-- write once
!=============
  if(lfirst) then
     open(unit=4,file='output.grd_grid',status='replace',recl=max(recleng,reclenf))
!-- header: dimension
     write(4,*) "#",grd_igeom
     write(4,*) "#",grd_nx,grd_ny,grd_nz
     write(4,*) "#",grd_ncp,nmpi
!-- body
     write(4,'(1p,10000e12.4)') grd_xarr(:)
     write(4,'(1p,10000e12.4)') grd_yarr(:)
     write(4,'(1p,10000e12.4)') grd_zarr(:)
!-- cell indices
     iarr = reshape(grd_icell,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(10000i12)') iarr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.flx_grid',status='replace',recl=reclenf)
!-- header: dimension
     write(4,*) "#",flx_ng,flx_nmu,flx_nom
!-- body
     write(4,'(1p,10000e12.4)') flx_wl(:)
     write(4,'(1p,100e12.4)') flx_mu(:)
     write(4,'(1p,100e12.4)') flx_om(:)
     close(4)
  endif

!
!-- scalars
!==========
  open(unit=4,file='output.tsp_time',status=fstat,position=pos)
  if(lfirst) write(4,*) "#",tsp_nt
  write(4,*) tsp_t
  close(4)

  open(unit=4,file='output.tot_energy',status=fstat,position=pos,recl=(6+4)*12)
  if(lfirst) write(4,'("#",10a12)') 'error','erad','emat','eext','eout','evelo'
  write(4,'(1p,10e12.4)') tot_eerror,tot_erad,tot_emat,tot_eext,tot_eout,tot_evelo
  close(4)

!
!-- flux arrays
!==============
  open(unit=4,file='output.flx_luminos',status=fstat,position='append',recl=reclenf)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(1p,10000e12.4)') merge(flx_luminos(:,j,k),0d0,flx_luminos(:,j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
  enddo
  enddo
  close(4)

  open(unit=4,file='output.flx_lumnum',status=fstat,position='append',recl=reclenf)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(10000i12)') flx_lumnum(:,j,k)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.flx_lumdev',status=fstat,position='append',recl=reclenf)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(1p,10000e12.4)') flx_lumdev(:,j,k)/flx_luminos(:,j,k)
  enddo
  enddo
  close(4)
!
!-- gamma flux
  open(unit=4,file='output.flx_gamluminos',status=fstat,position='append',recl=reclenf)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(1p,10000e12.4)') merge(flx_gamluminos(j,k),0d0,flx_gamluminos(j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
  enddo
  enddo
  close(4)

!
!-- grid arrays
!==============
  if(.not.in_nogriddump) then
     iarr = reshape(grd_methodswap,[ncpr,nmpi])
     open(unit=4,file='output.grd_methodswap',status=fstat,position='append',recl=recleng)
     do i=1,nmpi
        write(4,'(10000i12)') iarr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.grd_temp',status=fstat,position='append',recl=recleng)
     arr = reshape(grd_temp,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(1p,10000e12.4)') arr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.grd_fcoef',status=fstat,position='append',recl=recleng)
     arr = reshape(grd_fcoef,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(1p,10000e12.4)') arr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.grd_eraddens',status=fstat,position='append',recl=recleng)
     arr = reshape(grd_eraddens/grd_vol,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(1p,10000e12.4)') arr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.grd_capgrey',status=fstat,position='append',recl=recleng)
     arr = reshape(grd_capgrey,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(1p,10000e12.4)') arr(:,i)
     enddo
     close(4)

     open(unit=4,file='output.grd_sig',status=fstat,position='append',recl=recleng)
     arr = reshape(grd_sig,[ncpr,nmpi])
     do i=1,nmpi
        write(4,'(1p,10000e12.4)') arr(:,i)
     enddo
     close(4)
  endif

!
!-- after the first iteration open files in append mode
  lfirst = .false.
  pos = 'append'
  fstat = 'old'


end subroutine write_output
