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
      real*8,parameter :: man_aa11 = 1.371d14*2.997924562d10 !erg/cm^2/s
      real*8,parameter :: man_aa22 = 1.371d12*2.997924562d10 !erg/cm^2/s
c
c-- a uniform temperature value (or possibly nominal temperature value)
!      real*8,parameter :: man_temp0 = 1.1602621d7 !K
      real*8,parameter :: man_temp0 = 1.160237998048407d7 !K
c
c
      contains
c
c
c
      subroutine check_manufacpars
c     ----------------------------
      use miscmod, only:warn
      use inputparmod, only: in_nobbopac,in_nobfopac,in_noffopac,
     &     in_nothmson,in_dentype,in_opacanaltype
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
      if(in_dentype/='unif')
     &     call warn('check_manufacpars','in_dentype/=unif')
c
      end subroutine check_manufacpars
c
c
      subroutine generate_manuradsrc(totmass,sigcoef,texp,dt)
c     ------------------------------
      use miscmod, only:warn
      use physconstmod
      use gasgridmod
      implicit none
      real*8,intent(in) :: totmass,sigcoef,texp,dt
************************************************************************
* calculate finite volume manufactured radiation source per cell
* per group in ergs with manufactured parameters
************************************************************************
      integer :: ir, ig
      real*8 :: x1,x2,x3,x4,xx3,xx4
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(gas_isvelocity) then
c         
c-- implement/modify velocity dependent manufactured radiation source
         select case (gas_opacanaltype)
         case ('grey')
c-- grey solution
            do ir = 1, gas_nr
               do ig = 1, gas_ng
                  gas_emitex(ig,ir)= (1d0/dt)*(
     &                 log((texp+dt)/texp)
     &                 *(3d0*man_aa11/pc_c)+
     &                 (3d0*totmass*sigcoef/
     &                 (8d0*pc_pi*gas_velout))*
     &                 ((gas_velout*texp)**(-2d0)-
     &                 (gas_velout*(texp+dt))**(-2d0))*
     &                 (man_aa11-pc_acoef*pc_c*man_temp0**4)
     &                 )*(x4-x3)/(x2-x1)
               enddo
!     
               gas_emitex(:,ir) = gas_emitex(:,ir)*
     &              gas_vals2(ir)%vol*dt
!     
            enddo
c--
         case ('mono')
            stop 'generate_manuradsrc: gas_opacanaltype=mono'
         case ('pick')
            stop 'generate_manuradsrc: gas_opacanaltype=pick'
         case ('line')
c-- line solution
            if(gas_ng/=2)
     &           stop 'generate_manuradsrc: gas_opacanaltype=line'
            if(gas_ldisp1/gas_ldisp2>=1d-3)
     &           stop 'generate_manuradsrc: gas_ldisp1/gas_ldisp2>=1d-3'
c
            do ir = 1, gas_nr           
!-- thin lines
               do ig = 1, gas_ng, 2
                  x3 = 1d0/gas_wl(ig+1)
                  x4 = 1d0/gas_wl(ig)
!-- calculating manufactured source
!     xx3 = x3*pc_h*pc_c/(pc_kb*man_temp0)
!     xx4 = x4*pc_h*pc_c/(pc_kb*man_temp0)              
!     bspeced = 15d0*specint(xx3,xx4,3)/pc_pi**4
!     
                  gas_emitex(ig,ir)=(1d0/dt)*
     &                 log((texp+dt)/texp)*
     &                 (man_aa11/pc_c)*
     &                 (1.5d0)
!     write(*,*) x3/(x4-x3)
!
!     gas_emitex(ig,ir)=0d0
               enddo
!
!-- thick lines
               do ig = 2, gas_ng, 2
                  x3 = 1d0/gas_wl(ig)
                  x4 = 1d0/gas_wl(ig-1)
                  xx3 = x3*pc_h*pc_c/(pc_kb*man_temp0)
                  gas_emitex(ig,ir)=(1d0/dt)*
     &                 log((texp+dt)/texp)*(man_aa11/pc_c)*
     &                 (2d0-0.5d0*x3/(x4-x3))
!     write(*,*) x3/(x4-x3)
               enddo
!     
!
               gas_emitex(:,ir) = gas_emitex(:,ir)*
     &              gas_vals2(ir)%vol*dt
!     
!     
            enddo
c     
         case default
            stop 'gas_opacanaltype unknown'
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
c
      end subroutine generate_manuradsrc
c
c
      subroutine generate_manutempsrc(totmass,sigcoef,texp,dt)
c     -------------------------------
      use physconstmod
      use gasgridmod
      implicit none
      real*8,intent(in) :: totmass,sigcoef,texp,dt
************************************************************************
* calculate finite volume manufactured temperature source
* (gas_vals2%matsrc) in ergs/cm^3/s with manufactured parameters
************************************************************************
      integer :: ir, ig
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(gas_isvelocity) then
c
c-- implement/modify velocity dependent manufactured temperature source
         select case (gas_opacanaltype)
         case ('grey')
c--   grey solution
            gas_vals2%matsrc = (1d0/dt)*
     &           (3d0*totmass*sigcoef/(8d0*pc_pi*gas_velout))*
     &           ((gas_velout*texp)**(-2d0)-
     &           (gas_velout*(texp+dt))**(-2d0))*
     &           (pc_acoef*pc_c*man_temp0**4d0-man_aa11)
c
         case ('mono')
            stop 'generate_manutempsrc: gas_opacanaltype=mono'
         case ('pick')
            stop 'generate_manutempsrc: gas_opacanaltype=pick'
         case ('line')
c--   line solution
            gas_vals2%matsrc = 0d0 !already set zero in gasgridmod
c
         case default
            stop 'gas_opacanaltype unknown'
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
c
      end subroutine generate_manutempsrc
c
c
      subroutine init_manuprofile
c     ---------------------------
      use physconstmod
      use gasgridmod      
      implicit none
************************************************************************
* calculate finite volume manufactured initial energy per cell per group
* in ergs with manufactured parameters
************************************************************************
c
c-- verify applicable input pars
      call check_manufacpars
c
c-- determine manufacture type
      if(gas_isvelocity) then
c
c-- implement/modify velocity dependent manufactured temperature source
         select case (gas_opacanaltype)
         case ('grey')
c-- grey solution
         case ('mono')
            stop 'init_manuprofile: gas_opacanaltype=mono'
         case ('pick')
            stop 'init_manuprofile: gas_opacanaltype=pick'
         case ('line')
c-- line solution
         case default
            stop 'gas_opacanaltype unknown'
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
c
      end subroutine init_manuprofile
c
      end module manufacmod
