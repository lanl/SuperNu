subroutine write_output

  use inputparmod
  use timingmod
  use timestepmod
  use gridmod
  use totalsmod
  use fluxmod
  implicit none

  integer :: j,k
  integer :: reclenf, recleng
  logical,save :: lfirst=.true.
  character(16), save :: pos='rewind', fstat='replace'
!
  reclenf = (flx_ng+1)*12
  recleng = (max(grd_nx,grd_ny,grd_nz) + 1)*12

!
!-- write once
!=============
  if(lfirst) then
     open(unit=4,file='output.grd_grid',status='replace',recl=recleng)
!-- header: dimension
     write(4,*) "#",grd_igeom
     write(4,*) "#",grd_nx,grd_ny,grd_nz
!-- body
     write(4,'(1p,10000e12.4)') grd_xarr(:)
     write(4,'(1p,10000e12.4)') grd_yarr(:)
     write(4,'(1p,10000e12.4)') grd_zarr(:)
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
     open(unit=4,file='output.grd_methodswap',status=fstat,position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(10000i12)') grd_methodswap(:,j,k)
     enddo
     enddo
     close(4)

     open(unit=4,file='output.grd_temp',status=fstat,position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(1p,10000e12.4)') grd_temp(:,j,k)
     enddo
     enddo
     close(4)

     open(unit=4,file='output.grd_fcoef',position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(1p,10000e12.4)') grd_fcoef(:,j,k)
     enddo
     enddo
     close(4)

     open(unit=4,file='output.grd_eraddens',position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(1p,10000e12.4)') grd_eraddens(:,j,k)/grd_vol(:,j,k)
     enddo
     enddo
     close(4)

     open(unit=4,file='output.grd_capgrey',position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(1p,10000e12.4)') grd_capgrey(:,j,k)
     enddo
     enddo
     close(4)

     open(unit=4,file='output.grd_sig',position='append',recl=recleng)
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(1p,10000e12.4)') grd_sig(:,j,k)
     enddo
     enddo
     close(4)
  endif

!
!-- after the first iteration open files in append mode
  lfirst = .false.
  pos = 'append'
  fstat = 'old'


end subroutine write_output
