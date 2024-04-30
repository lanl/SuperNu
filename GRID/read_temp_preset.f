* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
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
      allocate(grd_temppreset(grd_ncell,tsp_nt))
      read(4,*,iostat=istat) grd_temppreset
      if(istat/=0) stop 'rd_temp_preset: file too short: input.temp'
c-- check EOF
      read(4,*,iostat=istat)
      if(istat==0) stop 'rd_temp_preset: file too long: input.temp'
      close(4)
      write(6,*) 'rd_temp_preset: custom temp profiles read'
c
      end subroutine read_temp_preset
c vim: fdm=marker
