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
  real*8 :: t0,t1
!
  t0 = t_time()
!
  if(lfirst) then
!
!-- shape of volume quantities
     if(.not.grd_lvoid) then
!-- no void cells
       ncpr = grd_nx
       nrow = grd_ncell/ncpr
     else
!-- find optimal row size
       npbot = grd_ncell !init with worst case
       do i=grd_nx,max(grd_nx,grd_ny,grd_nz,grd_ncell/2)+1
          nrow = ceiling(float(grd_ncell)/i)
          np = nrow*i - grd_ncell !number of pad cells
          if(np<npbot) then
             npbot = np
             ncpr = i !optimal row size
          endif
       enddo
       nrow = ceiling(float(grd_ncell)/ncpr)
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
  if(lfirst) write(4,*) tsp_t  !beginning of first time step
  write(4,*) tsp_t + tsp_dt  !end of time step
  close(4)

  open(unit=4,file='output.tot_energy',status=fstat,position=pos,recl=17*12)
  if(lfirst) then
     write(4,'("#",i5)') 16 !number of columns
     write(4,'("#",16a12)') 'eerror','erad','emat','eext','eout','evelo', &
        'sfluxgamma','sflux', &
        'sthermal','smanufac','sanalvol','sanalsurf','samp', &
        'sdecaygamma','sdecaybeta','sdeposgamma'
  endif
  write(4,'(1x,1p,16e12.4)') tot_eerror,tot_erad,tot_emat,tot_eext,tot_eout,tot_evelo, &
      tot_sfluxgamma,tot_sflux, &
      tot_sthermal,tot_smanufac,tot_sanalvol,tot_sanalsurf,tot_samp, &
      tot_sdecaygamma,tot_sdecaybeta,tot_sdeposgamma
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
     write(4,'(1p,10000e12.4)') flx_lumdev(:,j,k)
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
     allocate(iarr(ncpr*nrow),arr(ncpr*nrow))
     iarr(grd_ncell+1:) = 0
     arr(grd_ncell+1:) = 0d0
!
!-- grd input values
     iarr(:grd_ncell) = grd_nvol
     open(unit=4,file='output.grd_nvol',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     arr(:grd_ncell) = 1d0/grd_tempinv
     open(unit=4,file='output.grd_temp',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     arr(:grd_ncell) = grd_fcoef
     open(unit=4,file='output.grd_fcoef',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     arr(:grd_ncell) = grd_capgrey
     open(unit=4,file='output.grd_capgrey',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     arr(:grd_ncell) = grd_sig
     open(unit=4,file='output.grd_sig',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
     arr(:grd_ncell) = grd_tally(2,:)/grd_vol
     open(unit=4,file='output.grd_eraddens',status=fstat,position='append',recl=reclen)
     do i=1,nrow
        write(4,'(1p,10000e12.4)') arr((i-1)*ncpr+1:i*ncpr)
     enddo
     close(4)
!
!-- grd tally values
     if(in_io_dogrdtally) then
        iarr(:grd_ncell) = grd_methodswap
        open(unit=4,file='output.grd_methodswap',status=fstat,position='append',recl=reclen)
        do i=1,nrow
           write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
!
        iarr(:grd_ncell) = grd_numcensimc
        open(unit=4,file='output.grd_numcensimc',status=fstat,position='append',recl=reclen)
        do i=1,nrow
           write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
!
        iarr(:grd_ncell) = grd_numcensddmc
        open(unit=4,file='output.grd_numcensddmc',status=fstat,position='append',recl=reclen)
        do i=1,nrow
           write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
!
        iarr(:grd_ncell) = grd_numfluxorig
        open(unit=4,file='output.grd_numfluxorig',status=fstat,position='append',recl=reclen)
        do i=1,nrow
           write(4,'(10000i12)') iarr((i-1)*ncpr+1:i*ncpr)
        enddo
        close(4)
     endif !in_io_dogrdtally
!
     deallocate(iarr,arr)
  endif

!
!-- after the first iteration open files in append mode
  lfirst = .false.
  pos = 'append'
  fstat = 'old'
!
!-- timing
  t1 = t_time()
  call timereg(t_output, t1-t0)

end subroutine write_output
