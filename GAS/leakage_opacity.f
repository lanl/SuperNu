      subroutine leakage_opacity
c     --------------------------
      use gridmod
      implicit none
************************************************************************
* wrapper
************************************************************************
      select case(grd_igeom)
      case(1)
         call leakage_opacity1
      case(2)
         call leakage_opacity2
      case(3)
         call leakage_opacity3
      case default
         stop 'leakage_opacity: invalid igeom'
      endselect
c
      end subroutine leakage_opacity
