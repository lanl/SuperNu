*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2022 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      subroutine banner
c     -----------------
      use mpimod
      implicit none
************************************************************************
* Print banner, start date/time, and code revision.
************************************************************************
      character*8 :: t_startdate
      character*10 :: t_starttime
c
      character*(MPI_MAX_PROCESSOR_NAME) :: pname
      integer :: ilen,ierr
c
      character*8  :: coderev_nr=''   !value set in version.inc
      character*40 :: coderev_id=''   !value set in version.inc
      character*24 :: coderev_date='' !value set in version.inc
      character*28 :: build_date=''   !value set in version.inc
      include 'version.inc'
c
      call date_and_time(t_startdate,t_starttime)
      call mpi_get_processor_name(pname,ilen,ierr)
c
      write(6,'(13(1x,"===",a62,"===",/))')
     & "==============================================================",
     & "==============================================================",
     & "                                                              ",
     & "       #####                              #     #             ",
     & "      #     # #    # #####  ###### #####  ##    # #    #      ",
     & "      #       #    # #    # #      #    # # #   # #    #      ",
     & "       #####  #    # #    # #####  #    # #  #  # #    #      ",
     & "            # #    # #####  #      #####  #   # # #    #      ",
     & "      #     # #    # #      #      #   #  #    ## #    #      ",
     & "       #####   ####  #      ###### #    # #     #  ####       ",
     & "                                                              ",
     & "==============================================================",
     & "=============================================================="
c
      write(6,*) "Radiation Transport Code"
      write(6,*) "by Ryan T. Wollaeger and Daniel R. van Rossum"
      write(6,*) "(2013/mar/05)"
      write(6,*)
c
      write(6,*) "code revision nr: ", trim(coderev_nr)
      write(6,*) "code revision id: ", trim(coderev_id)
      write(6,*) "code revis. date: ", trim(coderev_date)
      write(6,*) "build date:       ", trim(build_date)
c
c     write(6,*) "build hostname  : ", trim(build_hostname)
c     write(6,*) "compiler        : ", trim(compiler_version)
c     write(6,*) "compiler flags  : ", trim(compiler_flags)
c
      write(6,*) 'simulation date : ', t_startdate//" / "//t_starttime
      write(6,*) 'processor name  : ', trim(pname)
      write(6,*)
c
      end subroutine banner
c vim: fdm=marker
