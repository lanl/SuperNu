      module physconstmod
c     -------------------
      implicit none
************************************************************************
* physical constants and conversion factors
************************************************************************
c
c-- LOOK OUT KEV UNITS!!!
      real*8,parameter :: pc_acoef = 1.371e14  !radiation constant (erg/Kev^4/cm^3)
c
c
c-- physical constants
      real*8,parameter :: pc_pi = 3.1415926535897932385d0
      real*8,parameter :: pc_pi2 = 2d0*pc_pi
      real*8,parameter :: pc_pi4 = 4d0*pc_pi
c
      real*8,parameter :: pc_c = 2.997924562d10 !cm/s
      real*8,parameter :: pc_gc = 6.6742d-8
      real*8,parameter :: pc_rg = 8.314472d7
      real*8,parameter :: pc_kb = 1.3806505d-16 !erg/K
      real*8,parameter :: pc_me = 9.1093826d-28 !g
      real*8,parameter :: pc_mh = 1.67333d-24 !g
      real*8,parameter :: pc_amu = 1.66053886d-24 !g
      real*8,parameter :: pc_navo = 1d0/pc_amu  ! Avogadro's number
      real*8,parameter :: pc_h = 6.6260693d-27 !erg*s
      real*8,parameter :: pc_e = 4.80325d-10
      real*8,parameter :: pc_abohr = 0.5291772108d-8 !cm
      real*8,parameter :: pc_rydberg = 1.0973731534d5 !cm^-1
c
      real*8,parameter :: pc_sb = 2*pc_pi**5*pc_kb**4/
     &  (15d0*pc_h**3*pc_c**2) !5.670400d-5 !erg/cm^2/s/K^4
      real*8,parameter :: pc_hsun = 3.826d33/(pc_pi4*pc_pi4) !erg/s/cm^2/sr
      real*8,parameter :: pc_msun = 1.989d33 !g
c
c-- conversion factors
      real*8,parameter :: pc_day = 24d0*60d0**2
      real*8,parameter :: pc_year = 31556925.9747d0 !sec
      real*8,parameter :: pc_ev = 1.60217653d-12 !erg
      real*8,parameter :: pc_km = 1d5 !cm
      real*8,parameter :: pc_ang = 1d-8 !cm
      real*8,parameter :: pc_mbarn = 1d-18 !cm^2 (Mega barn)
c
c-- nuclear data: http://ie.lbl.gov/
c-- nuclear decay half-life times
cbad  real*8,parameter :: pc_thl_ni56 = 5.25053d5 !sec, 6.077 days
cbad  real*8,parameter :: pc_thl_co56 = 6.67613d6 !sec, 77.27 days
      real*8,parameter :: pc_thl_ni56 = 7.575d5 !sec, 6.077 days, converted to 1/exp life
      real*8,parameter :: pc_thl_co56 = 9.632d6 !sec, 77.27 days, converted to 1/exp life
c-- nuclear decay gamma emission, 1MeV=1.602176d-6 erg
c     real*8,parameter :: pc_qhl_ni56 = 3.421d-6 !ergs, 2135 keV, total decay energy
c     real*8,parameter :: pc_qhl_co56 = 7.316d-6 !ergs, 4566 keV
      real*8,parameter :: pc_qhl_ni56 = 2.756d-6 !ergs, 1720 keV, total gamma ray production
      real*8,parameter :: pc_qhl_co56 = 5.815d-6 !ergs, 3440 + .19*2*511 kev !direct gamma's plus positron annihilation
c-- average kinetic energy of co56->fe56 emergent positron
      real*8,parameter :: pc_q_poskin = 0.035d-6 !ergs, .19*116 kev
c-- gamma ray gray absorption cross section
      real*8,parameter :: pc_capgam = .06d0 !cm^2/g
c
c
      contains
c
c
c
      pure function planck(wl,temp) result(b)
c     ---------------------------------------
      implicit none
      real*8,intent(in) :: wl(:)
      real*8,intent(in) :: temp
      real*8 :: b(size(wl))
************************************************************************
* planck vector function
* note that wl input is in Anstrom
************************************************************************
      real*8,parameter :: c1 = 2d0*pc_h*pc_c**2
      real*8,parameter :: c2 = pc_h*pc_c/pc_kb
c-- standard form
c     b = wl*pc_ang
c     b = 2*pc_h*pc_c**2/(b**5*(exp(pc_h*pc_c/(b*pc_kb*temp)) - 1d0)) !in erg/cm^2/s/cm/ster
c-- optimized form
      b = 1d0/(wl*pc_ang*temp)
      b = c1*(temp*b)**5/(exp(c2*b) - 1d0)
      end function planck
c
c
c
      pure function dplanckdtemp(wl,temp) result(b)
c     ---------------------------------------------
      implicit none
      real*8,intent(in) :: wl(:)
      real*8,intent(in) :: temp
      real*8 :: b(size(wl))
************************************************************************
* dB/dT vector function
* note that wl input is in Anstrom
************************************************************************
      real*8,parameter :: c1 = 2*pc_h*pc_c**2
      real*8,parameter :: c2 = pc_h*pc_c/pc_kb
c
c-- normal form
c     b = pc_h*pc_c/(wl*pc_ang*pc_kb*temp)
c     b = 2*pc_kb**5*temp**4/(pc_h**4*pc_c**3) * b**6*exp(b)/(exp(b) - 1d0)**2   !in erg/cm^2/s/cm/ster/K
c-- optimized form
      b = 1d0/(wl*pc_ang*temp)
      b = c1*c2*(temp*b)**6*exp(c2*b)/(temp*(exp(c2*b) - 1d0))**2
      end function dplanckdtemp
c
      end module physconstmod
