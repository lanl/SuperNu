      subroutine leakage_opacity(icell1,ncell)
c     ----------------------------------------
      use gridmod
      use timingmod
      implicit none
      integer,intent(in) :: icell1,ncell
************************************************************************
* wrapper
************************************************************************
      real*8 :: t0,t1
c
      t0 = t_time()
c
      select case(grd_igeom)
      case(1)
         call leakage_opacity1(icell1,ncell)
      case(2)
         call leakage_opacity2(icell1,ncell)
      case(3)
         call leakage_opacity3(icell1,ncell)
      case(11)
         call leakage_opacity11(icell1,ncell)
      case default
         stop 'leakage_opacity: invalid igeom'
      endselect
c
      t1 = t_time()
      call timereg(t_opacleak,t1-t0) 
c
      end subroutine leakage_opacity
