module particlemod

  implicit none

  !Ryan W.: Changing group attribute to continuous wavelength (rev. 120)
  type packet
     integer :: zsrc, rtsrc !,gsrc
     real*8 :: rsrc, musrc, tsrc
     real*8 :: Esrc, Ebirth, wlsrc
     logical :: isvacant
  end type packet
  TYPE(packet), dimension(:), pointer :: prt_particles  !(prt_npartmax)

  integer :: prt_npartmax, prt_ns, prt_ninit
  integer :: prt_nsurf, prt_nexsrc, prt_nnew, prt_ninitnew
!-- rtw: random number counter added (rev. 262). associated with particle routines
  integer :: prt_tlyrnd
  integer, dimension(:), allocatable :: prt_vacantarr

  logical :: prt_done
  logical :: prt_isimcanlog !sets flux tally and energy deposition ...
  !to analog in IMC
  logical :: prt_isddmcanlog !sets flux tally and energy deposition ...
  !to analog in DDMC

  real*8 :: prt_tauddmc

  save

  contains

  subroutine particle_init(npartmax,ns,ninit,isimcanlog,isddmcanlog,tauddmc)
!--------------------------------------
    integer,intent(in) :: npartmax, ns, ninit
    logical,intent(in) :: isimcanlog, isddmcanlog
    real*8,intent(in) :: tauddmc
!***********************************************************************
! init particle module
!***********************************************************************
!
!-- adopt input values in module internal storage
    prt_npartmax = npartmax
    prt_ns = ns
    prt_ninit = ninit
    prt_isimcanlog = isimcanlog
    prt_isddmcanlog = isddmcanlog
    prt_tauddmc = tauddmc
!
!-- allocate permanent storage (dealloc in dealloc_all.f)
    allocate(prt_particles(prt_npartmax))
    prt_particles%isvacant = .true.
!
  end subroutine particle_init

end module particlemod
