      subroutine leakage_opacity
c     --------------------------
      use gridmod
      use timingmod
      implicit none
************************************************************************
* wrapper
************************************************************************
      real*8 :: t0,t1
c
      call time(t0)
c
      select case(grd_igeom)
      case(1)
         call leakage_opacity1
      case(2)
         call leakage_opacity2
      case(3)
         call leakage_opacity3
      case(4)
         call leakage_opacity11
      case default
         stop 'leakage_opacity: invalid igeom'
      endselect
c
      call time(t1)
      call timereg(t_opacleak,t1-t0) 
c
      end subroutine leakage_opacity
