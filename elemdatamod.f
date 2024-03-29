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
      module elemdatamod
c     ------------------
      implicit none
************************************************************************
*  important chemical element data
************************************************************************
      integer,parameter :: elem_neldata = 111
      integer,private :: i
c
      type element_data
       real*8       :: m
       character*2  :: sym
       character*14 :: nam
      end type element_data
      type(element_data) :: elem_data(elem_neldata)
c
      data (elem_data(i)%m,elem_data(i)%sym,elem_data(i)%nam,
     &  i=1,elem_neldata)
     &/   1.0079, 'H ', 'Hydrogen'      !  1,
     &,   4.0026, 'He', 'Helium'        !  2,
     &,   6.941,  'Li', 'Lithium'       !  3,
     &,   9.0122, 'Be', 'Beryllium'     !  4,
     &,  10.811,  'B ', 'Boron'         !  5,
     &,  12.0107, 'C ', 'Carbon'        !  6,
     &,  14.0067, 'N ', 'Nitrogen'      !  7,
     &,  15.9994, 'O ', 'Oxygen'        !  8,
     &,  18.9984, 'F ', 'Fluorine'      !  9,
     &,  20.1797, 'Ne', 'Neon'          ! 10,
     &,  22.9897, 'Na', 'Sodium'        ! 11,
     &,  24.305,  'Mg', 'Magnesium'     ! 12,
     &,  26.9815, 'Al', 'Aluminum'      ! 13,
     &,  28.0855, 'Si', 'Silicon'       ! 14,
     &,  30.9738, 'P ', 'Phosphorus'    ! 15,
     &,  32.065,  'S ', 'Sulfur'        ! 16,
     &,  35.453,  'Cl', 'Chlorine'      ! 17,
     &,  39.948,  'Ar', 'Argon'         ! 18,
     &,  39.0983, 'K ', 'Potassium'     ! 19,
     &,  40.078,  'Ca', 'Calcium'       ! 20,
     &,  44.9559, 'Sc', 'Scandium'      ! 21,
     &,  47.867,  'Ti', 'Titanium'      ! 22,
     &,  50.9415, 'V ', 'Vanadium'      ! 23,
     &,  51.9961, 'Cr', 'Chromium'      ! 24,
     &,  54.938,  'Mn', 'Manganese'     ! 25,
     &,  55.845,  'Fe', 'Iron'          ! 26,
     &,  58.9332, 'Co', 'Cobalt'        ! 27,
     &,  58.6934, 'Ni', 'Nickel'        ! 28,
     &,  63.546,  'Cu', 'Copper'        ! 29,
     &,  65.39,   'Zn', 'Zinc'          ! 30,
     &,  69.723,  'Ga', 'Gallium'       ! 31,
     &,  72.64,   'Ge', 'Germanium'     ! 32,
     &,  74.9216, 'As', 'Arsenic'       ! 33,
     &,  78.96,   'Se', 'Selenium'      ! 34,
     &,  79.904,  'Br', 'Bromine'       ! 35,
     &,  83.8,    'Kr', 'Krypton'       ! 36,
     &,  85.4678, 'Rb', 'Rubidium'      ! 37,
     &,  87.62,   'Sr', 'Strontium'     ! 38,
     &,  88.9059, 'Y ', 'Yttrium'       ! 39,
     &,  91.224,  'Zr', 'Zirconium'     ! 40,
     &,  92.9064, 'Nb', 'Niobium'       ! 41,
     &,  95.94,   'Mo', 'Molybdenum'    ! 42,
     &,  98,      'Tc', 'Technetium'    ! 43,
     &, 101.07,   'Ru', 'Ruthenium'     ! 44,
     &, 102.9055, 'Rh', 'Rhodium'       ! 45,
     &, 106.42,   'Pd', 'Palladium'     ! 46,
     &, 107.8682, 'Ag', 'Silver'        ! 47,
     &, 112.411,  'Cd', 'Cadmium'       ! 48,
     &, 114.818,  'In', 'Indium'        ! 49,
     &, 118.71,   'Sn', 'Tin'           ! 50,
     &, 121.76,   'Sb', 'Antimony'      ! 51,
     &, 127.6,    'Te', 'Tellurium'     ! 52,
     &, 126.9045, 'I ', 'Iodine'        ! 53,
     &, 131.293,  'Xe', 'Xenon'         ! 54,
     &, 132.9055, 'Cs', 'Cesium'        ! 55,
     &, 137.327,  'Ba', 'Barium'        ! 56,
     &, 138.9055, 'La', 'Lanthanum'     ! 57,
     &, 140.116,  'Ce', 'Cerium'        ! 58,
     &, 140.9077, 'Pr', 'Praseodymium'  ! 59,
     &, 144.24,   'Nd', 'Neodymium'     ! 60,
     &, 145,      'Pm', 'Promethium'    ! 61,
     &, 150.36,   'Sm', 'Samarium'      ! 62,
     &, 151.964,  'Eu', 'Europium'      ! 63,
     &, 157.25,   'Gd', 'Gadolinium'    ! 64,
     &, 158.9253, 'Tb', 'Terbium'       ! 65,
     &, 162.5,    'Dy', 'Dysprosium'    ! 66,
     &, 164.9303, 'Ho', 'Holmium'       ! 67,
     &, 167.259,  'Er', 'Erbium'        ! 68,
     &, 168.9342, 'Tm', 'Thulium'       ! 69,
     &, 173.04,   'Yb', 'Ytterbium'     ! 70,
     &, 174.967 , 'Lu', 'Lutetium'      ! 71,
     &, 178.49,   'Hf', 'Hafnium'       ! 72,
     &, 180.9479, 'Ta', 'Tantalum'      ! 73,
     &, 183.84,   'W ', 'Tungsten'      ! 74,
     &, 186.207,  'Re', 'Rhenium'       ! 75,
     &, 190.23,   'Os', 'Osmium'        ! 76,
     &, 192.217,  'Ir', 'Iridium'       ! 77,
     &, 195.078,  'Pt', 'Platinum'      ! 78,
     &, 196.9665, 'Au', 'Gold'          ! 79,
     &, 200.59,   'Hg', 'Mercury'       ! 80,
     &, 204.3833, 'Tl', 'Thallium'      ! 81,
     &, 207.2,    'Pb', 'Lead'          ! 82,
     &, 208.9804, 'Bi', 'Bismuth'       ! 83,
     &, 209,      'Po', 'Polonium'      ! 84,
     &, 210,      'At', 'Astatine'      ! 85,
     &, 222,      'Rn', 'Radon'         ! 86,
     &, 223,      'Fr', 'Francium'      ! 87,
     &, 226,      'Ra', 'Radium'        ! 88,
     &, 227,      'Ac', 'Actinium'      ! 89,
     &, 232.0381, 'Th', 'Thorium'       ! 90,
     &, 231.0359, 'Pa', 'Protactinium'  ! 91,
     &, 238.0289, 'U ', 'Uranium'       ! 92,
     &, 237,      'Np', 'Neptunium'     ! 93,
     &, 244,      'Pu', 'Plutonium'     ! 94,
     &, 243,      'Am', 'Americium'     ! 95,
     &, 247,      'Cm', 'Curium'        ! 96,
     &, 247,      'Bk', 'Berkelium'     ! 97,
     &, 251,      'Cf', 'Californium'   ! 98,
     &, 252,      'Es', 'Einsteinium'   ! 99,
     &, 257,      'Fm', 'Fermium'       !100,
     &, 258,      'Md', 'Mendelevium'   !101,
     &, 259,      'No', 'Nobelium'      !102,
     &, 262,      'Lr', 'Lawrencium'    !103,
     &, 261,      'Rf', 'Rutherfordium' !104,
     &, 262,      'Db', 'Dubnium'       !105,
     &, 266,      'Sg', 'Seaborgium'    !106,
     &, 264,      'Bh', 'Bohrium'       !107,
     &, 277,      'Hs', 'Hassium'       !108,
     &, 268,      'Mt', 'Meitnerium'    !109,
     &, 275,      'Ds', 'Darmstadtium'  !110,
     &, 272,      'Rg', 'Roentgenium'   !111,
     &/
c
      save
      end module elemdatamod
c vim: fdm=marker
