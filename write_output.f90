subroutine write_output

  use timestepmod
  use gasgridmod
  use particlemod
  implicit none

  integer :: ir, ig
  character(16), save :: pos='rewind'

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_texp
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

  pos='append'

end subroutine write_output
