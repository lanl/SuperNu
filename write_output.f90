subroutine write_output

  use timingmod
  use timestepmod
  use gridmod
  use totalsmod
  use fluxmod
  implicit none

  integer :: j
  logical :: lexist
  integer :: reclenf, reclen2
  character(16), save :: pos='rewind', fstat='replace'
!
  reclenf = flx_ng*12
  reclen2 = grd_nx*12

  inquire(file='output.wlgrid',exist=lexist)
  if(.not.lexist) then
   open(unit=4,file='output.wlgrid',status='unknown',position=pos)
!-- header: dimension
   write(4,*) "#",flx_ng,flx_nmu,flx_nom
!-- body
   write(4,*) flx_wl
   write(4,*) flx_mu
   write(4,*) flx_om
   close(4)
  endif

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_t
  close(4)

  open(unit=4,file='output.pcktstat',status='unknown',position=pos)
  write(4,'(3es16.8)') t_pckt_stat
  close(4)

  open(unit=4,file='output.conserve',status='unknown',position=pos)
  write(4,*) tot_eerror
  close(4)

!
!-- flux arrays
  open(unit=4,file='output.LumNum',status=fstat,position='append',recl=reclenf)
  do i=1,flx_nmu
  do j=1,flx_nom
     write(4,'(10000i12)') flx_lumnum(:,i,j)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.Lum',status=fstat,position='append',recl=reclenf)
  where(flx_luminos<1d-99) flx_luminos = 1d-99  !prevent non-universal number representation, e.g. 1.1234-123
  do i=1,flx_nmu
  do j=1,flx_nom
     write(4,'(1p,10000e12.4)') flx_luminos
  enddo
  enddo
  close(4)

  open(unit=4,file='output.devLum',status=fstat,position='append',recl=reclenf)
  do i=1,flx_nmu
  do j=1,flx_nom
     write(4,'(1p,10000e12.4)') flx_lumdev/flx_luminos
  enddo
  enddo
  close(4)

!
!-- grid arrays
  open(unit=4,file='output.methodswap',status=fstat,position='append',recl=reclen2)
  do j=1,grd_ny
  do k=1,grd_nz
     write(4,'(10000i12)') grd_methodswap
  enddo
  enddo
  close(4)

  open(unit=4,file='output.temp',status=fstat,position='append',recl=reclen2)
  do j=1,grd_ny
  do k=1,grd_nz
     write(4,'(1p,10000e12.4)') grd_temp(:,j,k)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.grd_fcoef',position='append',recl=reclen2)
  do j=1,grd_ny
  do k=1,grd_nz
     write(4,'(1p,10000e12.4)') grd_fcoef(:,j,k)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.eraddens',position='append',recl=reclen2)
  do j=1,grd_ny
  do k=1,grd_nz
     write(4,'(1p,10000e12.4)') grd_eraddens(:,j,k)/grd_vol(:,j,k)
  enddo
  enddo
  close(4)

  open(unit=4,file='output.siggrey',position='append',recl=reclen2)
  do j=1,grd_ny
  do k=1,grd_nz
     write(4,'(1p,10000e12.4)') grd_capgrey(:,j,k)
  enddo
  enddo
  close(4)

!
!-- after the first iteration open files in append mode
  pos='append'
  fstat='old'


end subroutine write_output
