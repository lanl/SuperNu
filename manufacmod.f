*This file is part of SuperNu.  SuperNu is released under the terms of the GNU GPLv3, see COPYING.
*Copyright (c) 2013-2019 Ryan T. Wollaeger and Daniel R. van Rossum.  All rights reserved.
      module manufacmod
c     -----------------
      implicit none
************************************************************************
* any constants or factors from particular manufactured solutions
* to be used in source
************************************************************************
c
c
c-- max and min values for constant linear profile function
      real*8,parameter,private :: mnf_aa11 = 1.371d14*2.997924562d10 !erg/cm^2/s
      real*8,parameter,private :: mnf_aa22 = 1.371d12*2.997924562d10 !erg/cm^2/s
c
c-- a uniform temperature value (or possibly nominal temperature value)
!      real*8,parameter :: mnf_temp0 = 1.1602621d7 !K
      real*8,parameter,private :: mnf_temp0 = 1.160237998048407d7 !K
c
      save
c
      contains
c
c
c
      subroutine check_manufacpars
c     ----------------------------!{{{
      use miscmod, only:warn
      use inputparmod
      implicit none
************************************************************************
* check relevant input parameters for consistency with manufactured
* solution constraints
************************************************************************
c
      if(in_opacanaltype=='none')
     &     stop 'check_manufacpars: invalid in_opacanaltype'
      if(.not.in_nobbopac)
     &     stop 'check_manufacpars: invalid in_nobbopac'
      if(.not.in_nobfopac)
     &     stop 'check_manufacpars: invalid in_nobfopac'
      if(.not.in_noffopac)
     &     stop 'check_manufacpars: invalid in_noffopac'
      if(.not.in_nothmson)
     &     stop 'check_manufacpars: invalid in_nothmson'
      if(in_str_dentype/='unif')
     &     call warn('check_manufacpars','in_str_dentype/=unif')
c!}}}
      end subroutine check_manufacpars
c
c
      subroutine generate_manuradsrc(totmass,sigcoef,texp,dt)
c     ------------------------------!{{{
      use miscmod, only:warn
      use physconstmod
      use gridmod
      use gasmod
      use inputparmod
      implicit none
      real*8,intent(in) :: totmass,sigcoef,texp,dt
************************************************************************
* calculate finite volume manufactured radiation source per cell
* per group in ergs with manufactured parameters
************************************************************************
      integer :: i,l
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(in_isvelocity) then
c
c-- implement/modify velocity dependent manufactured radiation source
         select case (in_opacanaltype)
         case ('grey')
c-- grey solution
            do i=1,grd_nx
               l = grd_icell(i,1,1)
               grd_emitex(l)= (1d0/dt)*(
     &            log((texp+dt)/texp)
     &            *(4d0*mnf_aa11/pc_c)+
     &            (3d0*totmass*sigcoef/
     &            (8d0*pc_pi*in_str_velout))*
     &            ((in_str_velout*texp)**(-2d0)-
     &            (in_str_velout*(texp+dt))**(-2d0))*
     &            (mnf_aa11-pc_acoef*pc_c*mnf_temp0**4)
     &            )
!
            enddo
            grd_emitex = grd_emitex*grd_vol*dt
c--
         case ('mono')
            stop 'generate_manuradsrc: in_opacanaltype=mono'
         case ('pick')
            stop 'generate_manuradsrc: in_opacanaltype=pick'
         case ('line')
c-- line solution
            stop 'generate_manuradsrc: in_opacanaltype=line'
         case default
            stop 'in_opacanaltype unknown'
         end select
c
c
      else
c
c-- implement/modify static manufactured radiation source
         stop 'generate_manuradsrc: no static sources'
c
c
      endif
c!}}}
      end subroutine generate_manuradsrc
c
c
      subroutine generate_manutempsrc(totmass,sigcoef,texp,dt)
c     -------------------------------!{{{
      use physconstmod
      use gasmod
      use inputparmod
      implicit none
      real*8,intent(in) :: totmass,sigcoef,texp,dt
************************************************************************
* calculate finite volume manufactured temperature source
* (gas_matsrc) in ergs/cm^3/s with manufactured parameters
************************************************************************
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(in_isvelocity) then
c
c-- implement/modify velocity dependent manufactured temperature source
         select case (in_opacanaltype)
         case ('grey')
c--   grey solution
            gas_matsrc = (1d0/dt)*
     &           (3d0*totmass*sigcoef/(8d0*pc_pi*in_str_velout))*
     &           ((in_str_velout*texp)**(-2d0)-
     &           (in_str_velout*(texp+dt))**(-2d0))*
     &           (pc_acoef*pc_c*mnf_temp0**4d0-mnf_aa11)
c
         case ('mono')
            stop 'generate_manutempsrc: in_opacanaltype=mono'
         case ('pick')
            stop 'generate_manutempsrc: in_opacanaltype=pick'
         case ('line')
c--   line solution
            gas_matsrc = 0d0 !already set zero in gasmod
c
         case default
            stop 'in_opacanaltype unknown'
         end select
c
c
      else
c
c-- implement/modify static manufactured temperature source
         stop 'generate_manutempsrc: no static sources'
c
c
      endif
c!}}}
      end subroutine generate_manutempsrc
c
c
      subroutine init_manuprofile
c     ---------------------------
      use physconstmod
      use inputparmod
      use gridmod
      implicit none
************************************************************************
* calculate finite volume manufactured initial energy per cell per group
* in ergs with manufactured parameters
************************************************************************
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(.not.in_isvelocity) then
c-- implement/modify static manufactured temperature source
         stop 'init_manuprofile: no static sources'
      else
c
c-- implement/modify velocity dependent manufactured initial profile
         select case (in_opacanaltype)
         case ('grey')
c-- grey solution
           grd_evolinit = (mnf_aa11/pc_c) * grd_vol
c
         case ('mono')
            stop 'init_manuprofile: in_opacanaltype=mono'
         case ('pick')
            stop 'init_manuprofile: in_opacanaltype=pick'
         case ('line')
c-- line solution
           stop 'init_manuprofile: in_opacanaltype=line ! implemented'
         case default
            stop 'in_opacanaltype unknown'
         end select
      endif
c
      end subroutine init_manuprofile
c
c
      subroutine init_manutemp
c     ---------------------------
      use physconstmod
      use gridmod
      use inputparmod
      use gasmod
      use groupmod
      implicit none
************************************************************************
* calculate initial temperature in Kelvin
* with manufactured parameters
************************************************************************
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(in_isvelocity) then
c
c-- implement/modify velocity dependent manufactured initial profile
         select case (in_opacanaltype)
         case ('grey')
c-- grey solution (uniform nominal temperature)
            gas_temp = mnf_temp0
c
         case ('mono')
            stop 'init_manutemp: in_opacanaltype=mono'
         case ('pick')
            stop 'init_manutemp: in_opacanaltype=pick'
         case ('line')
c-- line solution
            if(grp_ng/=2)
     &           stop 'init_manutemp: in_opacanaltype=line'
            if(in_ldisp1/in_ldisp2>=1d-3)
     &           stop 'init_manutemp: in_ldisp1/in_ldisp2>=1d-3'
c
            gas_temp = mnf_temp0
c
         case default
            stop 'in_opacanaltype unknown'
         end select
c
c
      else
c
c-- implement/modify static manufactured initial temperature
         stop 'init_manutemp: no static solution'
c
      endif
c
      end subroutine init_manutemp
c
      end module manufacmod
c vim: fdm=marker
