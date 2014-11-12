subroutine write_output

  use timingmod
  use timestepmod
  use gridmod
  use totalsmod
  use fluxmod
  implicit none

  integer :: ig, j
  logical :: lexist
  integer :: reclen, reclen2
  character(16), save :: pos='rewind', fstat='replace'
!
  reclen = flx_ng*flx_nmu*flx_nom*12
  reclen2 = grd_nx*grd_ny*grd_nz*12

  inquire(file='output.wlgrid',exist=lexist)

  if(.not.lexist) then
   open(unit=4,file='output.wlgrid',status='unknown',position=pos)
!-- header: dimension
   write(4,*) "#",size(grd_wl)
!-- body
   write(4,*) grd_wl
   close(4)
  endif

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_t
  close(4)

  open(unit=4,file='output.LumNum',status=fstat,position='append',recl=reclen)
  write(4,'(10000i12)') flx_lumnum
  close(4)

  open(unit=4,file='output.methodswap',status=fstat,position='append',recl=reclen2)
  write(4,'(10000i12)') grd_methodswap
  close(4)

  open(unit=4,file='output.Lum',status=fstat,position='append',recl=reclen)
  where(flx_luminos<1d-99) flx_luminos = 1d-99  !prevent non-universal number representation, e.g. 1.1234-123
  write(4,'(1p,10000e12.4)') flx_luminos
  close(4)

  open(unit=4,file='output.devLum',status=fstat,position='append',recl=reclen)
  write(4,'(1p,10000e12.4)') flx_lumdev/flx_luminos
  close(4)

  open(unit=4,file='output.temp',status=fstat,position='append',recl=reclen2)
  write(4,'(1p,10000e12.4)') grd_temp
  close(4)

  open(unit=4,file='output.grd_fcoef',position='append',recl=reclen2)
  write(4,'(1p,10000e12.4)') grd_fcoef
  close(4)

  open(unit=4,file='output.eraddens',position='append',recl=reclen2)
  write(4,'(1p,10000e12.4)') grd_eraddens/grd_vol
  close(4)

  open(unit=4,file='output.conserve',status='unknown',position=pos)
  write(4,*) tot_eerror
  close(4)

  open(unit=4,file='output.siggrey',position='append',recl=reclen2)
  write(4,'(1p,10000e12.4)') grd_siggrey
  close(4)

  open(unit=4,file='output.pcktstat',status='unknown',position=pos)
  do j = 1, 3
     write(4,'(es16.8)',advance='no') t_pckt_stat(j)
  enddo
  close(4)

  pos='append'
  fstat='old'


end subroutine write_output
