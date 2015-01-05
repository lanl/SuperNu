      subroutine read_temp_preset
c     ---------------------------
      use gridmod
      use timestepmod
************************************************************************
* Read preset temperature history from file
* Used in 1D only, so far.
************************************************************************
      integer :: istat
c
      open(4,file='input.temp',status='old',iostat=istat)
      if(istat/=0) stop 'rd_temp_preset: no file: input.temp'
c-- alloc and read
      allocate(grd_temppreset(grd_nc,tsp_nt))
      read(4,*,iostat=istat) grd_temppreset
      if(istat/=0) stop 'rd_temp_preset: file too short: input.temp'
c-- check EOF
      read(4,*,iostat=istat)
      if(istat==0) stop 'rd_temp_preset: file too long: input.temp'
      close(4)
      write(6,*) 'rd_temp_preset: custom temp profiles read'
c
      end subroutine read_temp_preset
