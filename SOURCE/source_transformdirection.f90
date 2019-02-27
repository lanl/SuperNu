!This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
!Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
subroutine source_transformdirection

  use particlemod
  use sourcemod
  use inputparmod
  use physconstmod
  use transportmod

  implicit none

  integer :: ipart,ivac
!  real*8 :: om0, mu0, x0, y0, z0
!  real*8 :: om, mu
!  real*8 :: cmffact,mu1,mu2,gm

  type(packet),pointer :: p

!-- transform particle direction into lab frame
  do ipart=1,src_nnew
     ivac = src_ivacant(ipart)
     p => prt_particles(ivac)
     call direction2lab(p%x,p%y,p%z,p%mu,p%om)
  enddo

!!
!!-- transform particle direction into lab frame
!  do ipart=1,src_nnew
!     ivac = src_ivacant(ipart)
!!
!!-- particle properties from the array
!     x0 = prt_particles(ivac)%x
!     y0 = prt_particles(ivac)%y
!     z0 = prt_particles(ivac)%z
!     mu0 = prt_particles(ivac)%mu
!     om0 = prt_particles(ivac)%om
!!
!     call direction2lab(x0,y0,z0,mu0,om0)
!!
!!-- calculate transformation
!     select case(grd_igeom)
!     case(1,11)
!        cmffact = 1d0+mu0*x0/pc_c
!        mu = (mu0+x0/pc_c)/cmffact
!
!     case(2)
!        cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)/pc_c
!        gm = 1d0/sqrt(1d0-(x0**2+y0**2)/pc_c**2)
!!-- om
!        om = atan2(sqrt(1d0-mu0**2)*sin(om0), &
!             sqrt(1d0-mu0**2)*cos(om0)+(gm*x0/pc_c) * &
!             (1d0+gm*(cmffact-1d0)/(gm+1d0)))
!        if(om<0d0) om=om+pc_pi2
!!-- mu
!        mu = (mu0+(gm*y0/pc_c)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) / &
!             (gm*cmffact)
!
!     case(3)
!        mu1 = sqrt(1d0-mu0**2)*cos(om0)
!        mu2 = sqrt(1d0-mu0**2)*sin(om0)
!        cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)/pc_c
!!-- mu
!        mu = (mu0+z0/pc_c)/cmffact
!!-- om
!        om = atan2(mu2+y0/pc_c,mu1+x0/pc_c)
!        if(om<0d0) om = om+pc_pi2
!     endselect
!
!!-- in bounds
!     mu = min(mu,1d0)
!     mu = max(mu,-1d0)
!
!!-- save transformed values
!     prt_particles(ivac)%mu = mu
!     prt_particles(ivac)%om = om
!
!  enddo

end subroutine source_transformdirection
! vim: fdm=marker
