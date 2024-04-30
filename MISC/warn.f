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
      subroutine warn(source,mesg,sunt)
c     ---------------------------------
      implicit none
      character(*),intent(in) :: source
      character(*),intent(in) :: mesg
      character(*),intent(in),optional :: sunt
************************************************************************
* write warning message to multiple file units
************************************************************************
      integer :: i,iunt
      integer,save :: nwarn=0 !warning counter
      character(27),parameter :: fmt='(" WARNING #",i0,2(": ",a))'
c
      nwarn = nwarn + 1 !update counter
c
      if(.not. present(sunt)) then
       write(0,fmt) nwarn,source,mesg
       write(6,fmt) nwarn,source,mesg
      else
       do i=1,len(sunt)
        read(sunt(i:i),*) iunt
        write(iunt,fmt) nwarn,source,mesg
       enddo
      endif
      end subroutine warn
c vim: fdm=marker
