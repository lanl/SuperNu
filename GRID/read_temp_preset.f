      subroutine read_temp_preset
!-    ---------------------------!{{{
      use timestepmod
************************************************************************
* Read preset temperature history from file
* Used in 1D only, so far.
************************************************************************
      integer :: istat
!
      open(4,file='input.temp',status='old',iostat=istat)
      if(istat/=0) stop 'rd_temp_preset: no file: input.temp'
!-- alloc and read
      allocate(grd_temppreset(grd_nx,grd_ny,grd_nz,tsp_nt))
      read(4,*,iostat=istat) grd_temppreset
      if(istat/=0) stop 'rd_temp_preset: file too short: input.temp'
!-- check EOF
      read(4,*,iostat=istat)
      if(istat==0) stop 'rd_temp_preset: file too long: input.temp'
      close(4)
      write(6,*) 'rd_temp_preset: custom temp profiles read'
!}}}
      end subroutine read_temp_preset
