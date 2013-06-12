subroutine initialnumbers

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  implicit none

!##################################################
  !This subroutine computes the distribution of initial particles before
  !the first time step.  A fraction of the total initial particle number
  !is given to each cell based on the amount of inital radiative energy
  !profile.
!##################################################
!-------------
!Ryan W.(4/22/2013): currently only for manufactured solution
!-------------
  
  integer :: ir, ig, iir, iig, ipart
  integer,dimension(gas_nr) :: nvolinit, iirused
  integer :: nvolinittot, nvolinitapp
  real*8 :: wl0, mu0, Ep0, r0
  real*8 :: help
  real*8 :: r1,r2,r3
  real*8 :: exsumg, rrcenter, etotinit,denom2
  real*8 :: x1,x2,x3,x4
  real*8 :: aa11 = 1.371e14
  real*8 :: aa22 = 1.371e2
  real*8 :: uudd = 2.5d8
  logical :: isnotvacnt
  !
  nvolinittot = 5000*gas_nr
  nvolinitapp = 0
  etotinit = 0d0
  !
  gas_vals2(:)%eraddens=0d0
  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)

  if(gas_isvelocity) then
     help = gas_velout*tsp_texp
  else
     help = gas_l0+gas_lr
  endif

  if(gas_srctype=='manu') then
!{{{
     do ir = 1, gas_nr
        !
        rrcenter=(gas_rarr(ir+1)+gas_rarr(ir))/2d0
        do ig = 1, gas_ng, 2
           x3 = 1d0/gas_wl(ig+1)
           x4 = 1d0/gas_wl(ig)
           gas_vals2(ir)%eraddens=gas_vals2(ir)%eraddens+&
                ((x4-x3)/(x2-x1))*((aa11*(in_velout-rrcenter)+ &
                aa22*rrcenter)/in_velout)* &
                (in_tfirst*pc_day/tsp_texp)**4
        enddo

        do ig = 2, gas_ng, 2
           x3 = 1d0/gas_wl(ig+1)
           x4 = 1d0/gas_wl(ig)
           gas_vals2(ir)%eraddens=gas_vals2(ir)%eraddens+&
                ((x4-x3)/(x2-x1))*pc_acoef*gas_temp(ir)**4
        enddo
        
        etotinit = etotinit + gas_vals2(ir)%eraddens* &
             gas_vals2(ir)%volr*help**3
     enddo
     do ir = 1, gas_nr
        nvolinit(ir)=nint(gas_vals2(ir)%eraddens*gas_vals2(ir)%volr*help**3 &
             *nvolinittot/etotinit)
        nvolinitapp = nvolinitapp+nvolinit(ir)
     enddo!}}}
  endif

  !insantiating initial particles
  iir = 1
  iirused(1:gas_nr) = 0
  do ipart = 1, nvolinitapp
     isnotvacnt=.false.
     do while(.not.isnotvacnt)
        if(iirused(iir)<nvolinit(iir)) then
           iirused(iir)=iirused(iir)+1
           denom2 = 0d0
           r1=rand()
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1)
              x4 = 1d0/gas_wl(ig)
              iig = ig
              if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
              denom2 = denom2+(x4-x3)/(x2-x1)
           enddo
           !calculating wavelegth unformly
           r1 = rand()
           wl0 = (1d0-r1)*gas_wl(iig)+r1*gas_wl(iig+1)
           !calculating radial position
           r3 = rand()
           prt_particles(ipart)%rsrc = (r3*gas_rarr(iir+1)**3 + &
                (1.0-r3)*gas_rarr(iir)**3)**(1.0/3.0)
           r0 = prt_particles(ipart)%rsrc
           !calculating direction cosine (comoving)
           if(mod(iig,2)==0) then
              r1 = rand()
              mu0 = 1d0-2d0*r1
           else
              mu0=1d0
              !r1 = rand()
              !mu0 = 1d0-2d0*r1
           endif
           !calculating particle time
           prt_particles(ipart)%tsrc = tsp_time
           !calculating particle energy
           Ep0 = gas_vals2(iir)%eraddens*gas_vals2(iir)%volr*help**3 &
                /real(nvolinit(iir))
           !write(*,*) Ep0, gas_vals2(iir)%eraddens, nvolinit(iir)
           
           if(gas_isvelocity) then
              prt_particles(ipart)%Esrc = Ep0*(1.0+r0*mu0/pc_c)
              prt_particles(ipart)%Ebirth = Ep0*(1.0+r0*mu0/pc_c)
           
              prt_particles(ipart)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
           !
              prt_particles(ipart)%musrc = (mu0+r0/pc_c)/&
                   (1.0+r0*mu0/pc_c)
           else
              prt_particles(ipart)%Esrc = Ep0
              prt_particles(ipart)%Ebirth = Ep0
           
              prt_particles(ipart)%wlsrc = wl0
           !
              prt_particles(ipart)%musrc = mu0
           endif
           prt_particles(ipart)%rtsrc = 1
           
           !Setting ir = zone of particle
           prt_particles(ipart)%zsrc = iir
           !Setting particle index to not vacant
           prt_particles(ipart)%isvacant = .false.

           isnotvacnt = .true.
        else
           iir=iir+1
        endif
     enddo
  enddo
  !
  gas_vals2(:)%eraddens=0d0

end subroutine initialnumbers
