subroutine write_output

  use inputparmod
  use timingmod
  use timestepmod
  use gasgridmod
  use particlemod
  implicit none

  integer :: ir, ig, j
  integer :: reclen
  character(16), save :: pos='rewind', fstat='replace'
!
  reclen = gas_ng*12

  if(.not.in_isbdf2.or.(in_isbdf2.and.tsp_it>1)) then

  if(tsp_it==1) then
   open(unit=4,file='output.wlgrid',status='unknown',position=pos)
!-- header: dimension
   write(4,*) "#",size(gas_wl)
!-- body
   write(4,*) gas_wl
   close(4)
  endif

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_t
  close(4)

  open(unit=4,file='output.LumNum',status=fstat,position='append',recl=reclen)
  write(4,'(10000i12)') gas_lumnum
  close(4)

  open(unit=4,file='output.Lum',status=fstat,position='append',recl=reclen)
  where(gas_luminos<1d-99) gas_luminos = 1d-99  !prevent non-universal number representation, e.g. 1.1234-123
  write(4,'(1p,10000e12.4)') gas_luminos
  close(4)

  open(unit=4,file='output.devLum',status='unknown',position=pos)
  do ig = 1, gas_ng
     if(gas_luminos(ig)>0d0) then
        write(4,'(es16.8)',advance='no') gas_lumdev(ig)/gas_luminos(ig)
     else
        write(4,'(es16.8)',advance='no') 0d0
     endif
  enddo
  close(4)

  open(unit=4,file='output.temp',status=fstat,position='append',recl=reclen)
  write(4,'(1p,10000e12.4)') gas_tempold
  close(4)

  open(unit=4,file='output.gas_fcoef',status='unknown',position=pos)
  do ir = 1, gas_nr
    write(4,'(es16.8)',advance='no') gas_fcoef(ir)
  enddo
  close(4)

  open(unit=4,file='output.eraddens',status='unknown',position=pos)
  do ir = 1, gas_nr
    write(4,'(es16.8)',advance='no') gas_vals2(ir)%eraddens
  enddo
  close(4)

  open(unit=4,file='output.conserve',status='unknown',position=pos)
  write(4,*) gas_eerror
  close(4)

  open(unit=4,file='output.siggrey',status='unknown',position=pos)
  do ir = 1, gas_nr
    write(4,'(es16.8)',advance='no') gas_siggrey(ir)
  enddo
  close(4)

  open(unit=4,file='output.pcktstat',status='unknown',position=pos)
  do j = 1, 3
     write(4,'(es16.8)',advance='no') t_pckt_stat(j)
  enddo
  close(4)

  open(unit=4,file='output.opac1',status='unknown',position=pos)
  do ig = 1, gas_ng
     write(4,'(es16.8)',advance='no') gas_cap(ig,gas_nr/4+1)
  enddo
  close(4)
  open(unit=4,file='output.opac2',status='unknown',position=pos)
  do ig = 1, gas_ng
     write(4,'(es16.8)',advance='no') gas_cap(ig,gas_nr/2+1)
  enddo
  close(4)
  open(unit=4,file='output.opac3',status='unknown',position=pos)
  do ig = 1, gas_ng
     write(4,'(es16.8)',advance='no') gas_cap(ig,3*gas_nr/4+1)
  enddo
  close(4)

  pos='append'
  fstat='old'

  endif

end subroutine write_output
