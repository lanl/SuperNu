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
      module transportmod
c     -------------------
      use physconstmod, only:pc_c,pc_pi2
      implicit none
c
      private pc_c,pc_pi2
      real*8,private,parameter :: cinv=1d0/pc_c
c
      logical :: trn_nolumpshortcut !disable approximation for large emitlump that sampling outside the lump collapses to the single most likely group
      logical :: trn_errorfatal     !stop on transport error, disable for production runs
      logical :: trn_noampfact      !don't use the particle amplification factor
c
      logical :: trn_isimcanlog  !sets flux tally and energy deposition to analog in IMC
      logical :: trn_isddmcanlog !sets flux tally and energy deposition to analog in DDMC
c
      real*8 :: trn_tauddmc
      real*8 :: trn_taulump
      character(4) :: trn_tauvtime ! unif|incr
c
c-- explicit interfaces
      interface
c!{{{
c
c-- transport_gamgrey
      pure subroutine transport1_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep
      integer,intent(out) :: ierr
      end subroutine transport1_gamgrey
c
      pure subroutine transport2_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep
      integer,intent(out) :: ierr
      end subroutine transport2_gamgrey
c
      pure subroutine transport3_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep
      integer,intent(out) :: ierr
      end subroutine transport3_gamgrey
c
      pure subroutine transport11_gamgrey(ptcl,ptcl2,rndstate,edep,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep
      integer,intent(out) :: ierr
      end subroutine transport11_gamgrey
c
c-- transport
      pure subroutine transport1(ptcl,ptcl2,rndstate,
     &  edep,eraddens,eamp,totevelo,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens, eamp
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport1
c
      pure subroutine transport2(ptcl,ptcl2,rndstate,
     &  edep,eraddens,eamp,totevelo,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens, eamp
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport2
c
      pure subroutine transport3(ptcl,ptcl2,rndstate,
     &  edep,eraddens,eamp,totevelo,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens, eamp
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport3
c
      pure subroutine transport11(ptcl,ptcl2,rndstate,
     &  edep,eraddens,eamp,totevelo,ierr)
      use randommod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens, eamp
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport11
c
c-- diffusion
      pure subroutine diffusion1(ptcl,ptcl2,cache,rndstate,
     &  edep,eraddens,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(grp_t_cache),target,intent(inout) :: cache
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion1
c
      pure subroutine diffusion2(ptcl,ptcl2,cache,rndstate,
     &  edep,eraddens,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(grp_t_cache),target,intent(inout) :: cache
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion2
c
      pure subroutine diffusion3(ptcl,ptcl2,cache,rndstate,
     &  edep,eraddens,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(grp_t_cache),target,intent(inout) :: cache
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion3
c
      pure subroutine diffusion11(ptcl,ptcl2,cache,rndstate,
     &  edep,eraddens,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(grp_t_cache),target,intent(inout) :: cache
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion11
!}}}
      end interface
c
c-- abstract interfaces
      abstract interface
      subroutine direction2lab_(x0,y0,z0,mu0,om0)!{{{
c     -------------------------------------------
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      end subroutine direction2lab_
c
      pure subroutine transport_gamgrey_(ptcl,ptcl2,rndstate,edep,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep
      integer,intent(out) :: ierr
      end subroutine transport_gamgrey_
c
      pure subroutine transport_(ptcl,ptcl2,rndstate,
     &  edep,eraddens,eamp,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens, eamp
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine transport_
c
      pure subroutine diffusion_(ptcl,ptcl2,cache,rndstate,
     &  edep,eraddens,totevelo,ierr)
      use randommod
      use groupmod
      use particlemod
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      type(grp_t_cache),target,intent(inout) :: cache
      type(rnd_t),intent(inout) :: rndstate
      real*8,intent(out) :: edep, eraddens
      real*8,intent(inout) :: totevelo
      integer,intent(out) :: ierr
      end subroutine diffusion_
!}}}
      end interface
c
c-- procedure pointers
      procedure(direction2lab_),pointer :: direction2lab => null()
      procedure(transport_gamgrey_),pointer ::
     &  transport_gamgrey => null()
      procedure(transport_),pointer :: transport => null()
      procedure(diffusion_),pointer :: diffusion => null()
c
c-- private interfaces
      private diffusion11,diffusion1,diffusion2,diffusion3
      private transport11,transport1,transport2,transport3
      private transport11_gamgrey,transport1_gamgrey,transport2_gamgrey,
     &  transport3_gamgrey
c
      save
c
      contains
c
c
c
      subroutine transportmod_init(igeom)
c     -----------------------------------------------
      integer,intent(in) :: igeom
c
c-- set procedure pointers
      select case(igeom)
      case(1)
!      labfact => labfact1
!      cmffact => cmffact1
       direction2lab => direction2lab1
       transport_gamgrey => transport1_gamgrey
       transport => transport1
       diffusion => diffusion1
      case(2)
!      labfact => labfact2
!      cmffact => cmffact2
       direction2lab => direction2lab2
       transport_gamgrey => transport2_gamgrey
       transport => transport2
       diffusion => diffusion2
      case(3)
!      labfact => labfact3
!      cmffact => cmffact3
       direction2lab => direction2lab3
       transport_gamgrey => transport3_gamgrey
       transport => transport3
       diffusion => diffusion3
      case(11)
!      labfact => labfact1
!      cmffact => cmffact1
       direction2lab => direction2lab1
       transport_gamgrey => transport11_gamgrey
       transport => transport11
       diffusion => diffusion11
      case default
       stop 'transportmod_init: invalid igeom'
      end select
      end subroutine transportmod_init
c
c
c
      subroutine direction2lab1(x0,y0,z0,mu0,om0)
c     -------------------------------------------
      implicit none
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      real*8 :: cmffact,mu,dummy
c
      dummy = y0
      dummy = z0
      dummy = om0
c
      cmffact = 1d0+mu0*x0*cinv
      mu = (mu0+x0*cinv)/cmffact
      mu = min(mu,1d0)
      mu = max(mu,-1d0)
      mu0 = mu
      end subroutine direction2lab1
c
c
      subroutine direction2lab2(x0,y0,z0,mu0,om0)
c     -------------------------------------------
      implicit none
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      real*8 :: cmffact,gm,mu,om,dummy
c
      dummy = z0
c
      cmffact = 1d0+(mu0*y0+sqrt(1d0-mu0**2)*cos(om0)*x0)*cinv
      gm = 1d0/sqrt(1d0-(x0**2+y0**2)*cinv**2)
!-- om
      om = atan2(sqrt(1d0-mu0**2)*sin(om0),
     &     sqrt(1d0-mu0**2)*cos(om0)+(gm*x0*cinv) *
     &     (1d0+gm*(cmffact-1d0)/(gm+1d0)))
      if(om<0d0) om=om+pc_pi2
!-- mu
      mu = (mu0+(gm*y0*cinv)*(1d0+gm*(cmffact-1d0)/(1d0+gm))) /
     &     (gm*cmffact)
      mu = min(mu,1d0)
      mu = max(mu,-1d0)
      mu0 = mu
      om0 = om
      end subroutine direction2lab2
c
c
      subroutine direction2lab3(x0,y0,z0,mu0,om0)
c     -------------------------------------------
      implicit none
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      real*8 :: cmffact,mu1,mu2,mu,om
c
      mu2 = sqrt(1d0-mu0**2)
      mu1 = mu2*cos(om0)
      mu2 = mu2*sin(om0)
      cmffact = 1d0+(mu0*z0+mu1*x0+mu2*y0)*cinv
!-- mu
      mu = (mu0+z0*cinv)/cmffact
!-- om
      om = atan2(mu2+y0*cinv,mu1+x0*cinv)
      if(om<0d0) om = om+pc_pi2
!-- in bounds
      mu = min(mu,1d0)
      mu = max(mu,-1d0)
      mu0 = mu
      om0 = om
      end subroutine direction2lab3
c
c
      elemental function ddmc_emiss_bc(dx, f, cap, sig, dext)
c     -------------------------------------------------------
      implicit none
      real*8 :: ddmc_emiss_bc
      real*4,intent(in) :: cap
      real*8,intent(in) :: dx, f, sig, dext
************************************************************************
* See Densmore, Davidson and Carrington (2006)
************************************************************************
      real*8 :: help, alb, eps, beta
      help = (cap+sig)*dx
      alb = f*cap / (cap+sig)
      eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+dext*sqrt(3d0*alb))
      beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 +
     &     2.25d0*alb**2*help**4)
      ddmc_emiss_bc = 0.5d0*eps*beta/(beta-0.75*eps*help)
c
      end function ddmc_emiss_bc
c
c
      end module transportmod
c vim: fdm=marker
