subroutine write_output

  use timestepmod
  use gasgridmod
  use particlemod
  implicit none

  integer :: ir
  character(16), save :: pos='rewind'

  open(unit=4,file='output.tsp_time',status='unknown',position=pos)
  write(4,*) tsp_time
  close(4)

  open(unit=4,file='output.Lum',status='unknown',position=pos)
  write(4,*) gas_eright/tsp_dt
  close(4)

  open(unit=4,file='output.temp',status='unknown',position=pos)
  do ir = 1, gas_nr
    write(4,'(es16.8)',advance='no') gas_vals2(ir)%temp
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
