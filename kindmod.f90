MODULE kindmod

  INTEGER, PARAMETER :: rknd = SELECTED_REAL_KIND(14,100)
  INTEGER, PARAMETER :: iknd = SELECTED_INT_KIND(8)

  ! Derived data types
  TYPE packet
     INTEGER(iknd) :: zsrc, gsrc, rtsrc
     REAL(rknd) :: rsrc, musrc, tsrc
     REAL(rknd) :: Esrc, Ebirth
     LOGICAL :: isvacant
  END TYPE packet

END MODULE kindmod
