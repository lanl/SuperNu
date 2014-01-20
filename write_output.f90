subroutine write_output

  use timingmod
  use timestepmod
  use gasgridmod
  use particlemod
  implicit none

  integer :: ir, ig, j
  character(16), save :: pos='rewind'

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_t
  close(4)

  open(unit=4,file='output.Lum',status='unknown',position=pos)
  do ig = 1, gas_ng
     write(4,'(es16.8)',advance='no') gas_luminos(ig)
  enddo
  close(4)

  open(unit=4,file='output.temp',status='unknown',position=pos)
  do ir = 1, gas_nr
    write(4,'(es16.8)',advance='no') gas_temp(ir)
  enddo
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

  pos='append'

end subroutine write_output
