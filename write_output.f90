subroutine write_output

  use inputparmod
  use timingmod
  use timestepmod
  use gridmod
  use totalsmod
  use fluxmod
  implicit none

  integer :: i,j,k
  integer :: np,npbot
  integer,save :: nrow=0 !number of rows for printing grid variables
  integer,save :: ncpr=0 !number of cells per rank
  integer :: reclen
  logical,save :: lfirst=.true.
  character(16),save :: pos='rewind', fstat='replace'
!
  integer,allocatable :: iarr(:)
  real*8,allocatable :: arr(:)
!
  if(lfirst) then
!
!-- shape of volume quantities
     if(.not.grd_lpad) then
!-- no pad cells
       ncpr = grd_nx
       nrow = grd_nc/ncpr
     else
!-- find optimal row size
       npbot = grd_nc !init with worst case
       do i=grd_nx,max(grd_nx,grd_ny,grd_nz,grd_nc/2)+1
          nrow = ceiling(float(grd_nc)/i)
          np = nrow*i - grd_nc !number of pad cells
          if(np<npbot) then
             npbot = np
             ncpr = i !optimal row size
          endif
       enddo
       nrow = ceiling(float(grd_nc)/ncpr)
     endif

!
!-- write once
!=============
     reclen = max(4,grd_nx+1,grd_ny+1,grd_nz+1,grd_ny*grd_nz)*12
     open(unit=4,file='output.grd_grid',status='replace',recl=reclen)
!-- header: dimension
     write(4,*) "#",grd_igeom
     write(4,*) "#",grd_nx,grd_ny,grd_nz
     write(4,*) "#",nrow*ncpr,nrow,ncpr
!-- body
     write(4,'(1p,10000e12.4)') grd_xarr(:)
     write(4,'(1p,10000e12.4)') grd_yarr(:)
     write(4,'(1p,10000e12.4)') grd_zarr(:)
!-- cell indices
     do k=1,grd_nz
     do j=1,grd_ny
        write(4,'(10000i12)') grd_icell(:,j,k)
     enddo
     enddo
     close(4)

     reclen = max(4,flx_ng+1,flx_nmu+1,flx_nom+1)*12
     open(unit=4,file='output.flx_grid',status='replace',recl=reclen)
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

  open(unit=4,file='output.tot_energy',status=fstat,position=pos,recl=7*12)
  if(lfirst) write(4,'("#",6a12)') 'error','erad','emat','eext','eout','evelo'
  write(4,'(1p,10e12.4)') tot_eerror,tot_erad,tot_emat,tot_eext,tot_eout,tot_evelo
  close(4)

!
!-- flux arrays
!==============
  reclen = flx_ng*12
  open(unit=4,file='output.flx_luminos',status=fstat,position='append',recl=reclen)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(1p,10000e12.4)') merge(flx_luminos(:,j,k),0d0,flx_luminos(:,j,k)>1d-99) !prevent fortran number truncation, e.g. 1.1234-123
  enddo
  enddo
  close(4)

  open(unit=4,file='output.flx_lumnum',status=fstat,position='append',recl=reclen)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(10000i12)') flx_lumnum(:,j,k)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.flx_lumdev',status=fstat,position='append',recl=reclen)
  do k=1,flx_nom
  do j=1,flx_nmu
     write(4,'(1p,10000e12.4)') flx_lumdev(:,j,k)/flx_luminos(:,j,k)
  enddo
  enddo
  close(4)
!
!-- gamma flux
  open(unit=4,file='output.flx_gamluminos',status=fstat,position='append',recl=reclen)
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
     reclen = ncpr*12
!-- alloc
     allocate(iarr(ncpr*nrow))
     iarr(grd_nc+1:) = 0
!
     iarr(:grd_nc) = grd_methodswap(:grd_nc)
     open(unit=4,file='output.grd_methodswap',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
!-- alloc
     deallocate(iarr)
     allocate(arr(ncpr*nrow))
     arr(grd_nc+1:) = 0d0

     arr(:grd_nc) = grd_temp(:grd_nc)
     open(unit=4,file='output.grd_temp',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)

     arr(:grd_nc) = grd_fcoef(:grd_nc)
     open(unit=4,file='output.grd_fcoef',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)

     arr(:grd_nc) = grd_eraddens(:grd_nc)/grd_vol(:grd_nc)
     open(unit=4,file='output.grd_eraddens',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)

     arr(:grd_nc) = grd_capgrey(:grd_nc)
     open(unit=4,file='output.grd_capgrey',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)

     arr(:grd_nc) = grd_sig(:grd_nc)
     open(unit=4,file='output.grd_sig',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     deallocate(arr)
  endif

!
!-- after the first iteration open files in append mode
  lfirst = .false.
  pos = 'append'
  fstat = 'old'


end subroutine write_output
