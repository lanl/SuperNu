      module transportmod
c     -------------------
      use physconstmod, only:pc_c,pc_pi2
      implicit none
c
      real*8,private,parameter :: cinv=1d0/pc_c
c
      private pc_c,pc_pi2
c
      interface
      pure subroutine advection1(pretrans,ptcl,ptcl2)!{{{
c     -----------------------------------------------
      use particlemod
      logical,intent(in) :: pretrans
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      end subroutine advection1
c
      pure subroutine advection2(pretrans,ptcl,ptcl2)
c     -----------------------------------------------
      use particlemod
      logical,intent(in) :: pretrans
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      end subroutine advection2
c
      pure subroutine advection3(pretrans,ptcl,ptcl2)
c     -----------------------------------------------
      use particlemod
      logical,intent(in) :: pretrans
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      end subroutine advection3!}}}
      end interface
c
c
      abstract interface
      subroutine direction2lab_(x0,y0,z0,mu0,om0)!{{{
c     -------------------------------------------
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      end subroutine direction2lab_
c
      pure subroutine advection_(pretrans,ptcl,ptcl2)
c     -----------------------------------------------
      use particlemod
      logical,intent(in) :: pretrans
      type(packet),target,intent(inout) :: ptcl
      type(packet2),target,intent(inout) :: ptcl2
      end subroutine advection_!}}}
      end interface
c
c-- procedure pointers
      procedure(direction2lab_),pointer :: direction2lab => null()
      procedure(advection_),pointer :: advection => null()
c
      contains
c
c
c
      subroutine transportmod_init(igeom)
c     --------------------------------
      integer,intent(in) :: igeom
      select case(igeom)
      case(1,11)
!      labfact => labfact1
!      cmffact => cmffact1
       direction2lab => direction2lab1
       advection => advection1
      case(2)
!      labfact => labfact2
!      cmffact => cmffact2
       direction2lab => direction2lab2
       advection => advection2
      case(3)
!      labfact => labfact3
!      cmffact => cmffact3
       direction2lab => direction2lab3
       advection => advection3
      case default
       stop 'transportmod_init: invalid igeom'
      end select
      end subroutine transportmod_init
c
c
c
      subroutine direction2lab1(x0,y0,z0,mu0,om0)
c     ----------------------------------------
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
c     ----------------------------------------
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
c     ----------------------------------------
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
      end module transportmod
