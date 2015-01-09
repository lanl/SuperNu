      module transportmod
c     -------------------
      use physconstmod, only:pc_c,pc_pi2
      implicit none
c
      real*8,private,parameter :: cinv=1d0/pc_c
c
      private pc_c,pc_pi2
c
      abstract interface
!       pure function labfact_func(x,y,mu,om)
!       real*8 :: labfact_func
!       real*8,intent(in) :: x,y,mu,om
!       end function labfact_func
!c
!       pure function cmffact_func(x0,y0,z0,mu0,om0)
!       real*8 :: cmffact_func
!       real*8,intent(in) :: x0,y0,z0,mu0,om0
!       end function cmffact_func
c
      subroutine direction2lab_(x0,y0,z0,mu0,om0)
c     -------------------------------------------
      real*8,intent(in) :: x0,y0,z0
      real*8,intent(inout) :: mu0,om0
      end subroutine direction2lab_
c
      end interface
c
c-- procedure pointers
!      procedure(labfact_func),pointer :: labfact => null()
!      procedure(cmffact_func),pointer :: cmffact => null()
c
      procedure(direction2lab_),pointer :: direction2lab => null()
c
      contains
c
      subroutine transport_init(igeom)
c     --------------------------------
      integer,intent(in) :: igeom
      select case(igeom)
      case(1,11)
!      labfact => labfact1
!      cmffact => cmffact1
       direction2lab => direction2lab1
      case(2)
!      labfact => labfact2
!      cmffact => cmffact2
       direction2lab => direction2lab2
      case(3)
!      labfact => labfact3
!      cmffact => cmffact3
       direction2lab => direction2lab3
      case default
       stop 'transport_init: invalid igeom'
      end select
      end subroutine transport_init
c
c
!      pure function labfact1(x,y,mu,om) result(labfact)
!c     -------------------------------------------------
!      real*8 labfact
!      real*8,intent(in) :: x,y,mu,om
!      labfact = 1d0-x*mu/pc_c
!      end function labfact1
!c
!      pure function labfact2(x,y,mu,om) result(labfact)
!c     -------------------------------------------------
!      real*8 labfact
!      real*8,intent(in) :: x,y,mu,om
!      labfact = 1d0 - (mu*y + sqrt(1d0-mu**2) * cos(om)*x)/pc_c
!      end function labfact2
!c
!      pure function labfact3(x,y,mu,om) result(labfact)
!c     -------------------------------------------------
!      real*8 labfact
!      real*8,intent(in) :: x,y,mu,om
!      real*8 :: mu1,mu2
!      mu2 = sqrt(1d0-mu**2)
!      mu1 = mu2*cos(om)
!      mu2 = mu2*sin(om)
!      labfact = 1d0 - (mu*z + mu1*x + mu2*y)/pc_c
!      end function labfact3
!c
!c
!      pure function cmffact1(x0,y0,z0,mu0,om0) result(cmffact)
!c     -------------------------------------------------
!      real*8 cmffact
!      real*8,intent(in) :: x0,y0,z0,mu0,om0
!      cmffact = 1d0-x0*mu0/pc_c
!      end function cmffact1
!c
!      pure function cmffact2(x0,y0,z0,mu0,om0) result(cmffact)
!c     -------------------------------------------------
!      real*8 cmffact
!      real*8,intent(in) :: x0,y0,z0,mu0,om0
!      cmffact = 1d0 - (mu0*y0 + sqrt(1d0-mu0**2) * cos(om0)*x0)/pc_c
!      end function cmffact2
!c
!      pure function cmffact3(x0,y0,z0,mu0,om0) result(cmffact)
!c     -------------------------------------------------
!      real*8 cmffact
!      real*8,intent(in) :: x0,y0,z0,mu0,om0
!      real*8 :: mu1,mu2
!      mu2 = sqrt(1d0-mu0**2)
!      mu1 = mu2*cos(om0)
!      mu2 = mu2*sin(om0)
!      cmffact = 1d0 - (mu0*z + mu1*x0 + mu2*y0)/pc_c
!      end function cmffact3
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
