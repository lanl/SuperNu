subroutine interior_source

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod
  use manufacmod

  implicit none

!##################################################
  !This subroutine instantiates new volume (cell) particle properties.
  !Composed of external source particle loop (1st) and thermal source
  !particle loop (2nd).
!##################################################

  integer :: ir,irl,irr, ipart, ivac, ig, iig
  integer, dimension(gas_nr) :: irused
  real*8 :: r1, r2, r3, r4, uul, uur, uumax, mu0, r0, Ep0, wl0
  real*8 :: denom2,x1,x2,x3,x4, xx0, bmax, help
  real*8, dimension(gas_nr) :: exsumg
  logical :: isnotvacnt !checks for available particle space to populate in cell

  if(gas_isvelocity) then
     help = tsp_t
  else
     help = 1d0
  endif

  ir = 1
  irused(1:gas_nr) = 0
  !Volume particle instantiation: loop
  !Loop run over the number of new particles that aren't surface source
  !particles.
  
  x1=1d0/gas_wl(gas_ng+1)
  x2=1d0/gas_wl(1)
  do ipart = prt_nsurf+1, prt_nsurf+prt_nexsrc
     ivac = prt_vacantarr(ipart)
     isnotvacnt = .false.
     !If adding particle ivac in current cell ir does not exceed nvolex, add ivac to ir: loop
     do while (.not.isnotvacnt)
        if (irused(ir)<gas_nvolex(ir)) then
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0d0
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           do ig = 1, gas_ng
              x3=1d0/gas_wl(ig+1)
              x4=1d0/gas_wl(ig)
              iig = ig
              if(r1>=denom2.and.r1<denom2+(x4-x3)/(x2-x1)) exit
              denom2 = denom2+(x4-x3)/(x2-x1)
           enddo
           
           !Ryan W.: particle group removed (rev. 120)
           !prt_particles(ivac)%gsrc = iig
           !Calculating comoving wavelength uniformly from group

           r1 = rand()
           prt_tlyrand = prt_tlyrand+1

           wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))

           !Calculating radial position

           r3 = rand()
           prt_tlyrand = prt_tlyrand+1

           prt_particles(ivac)%rsrc = (r3*gas_rarr(ir+1)**3 + &
                (1.0-r3)*gas_rarr(ir)**3)**(1.0/3.0)
           r0 = prt_particles(ivac)%rsrc

           !Calculating direction cosine (comoving)
           !mu0 = 1d0
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu0 = 1d0-2d0*r1

           !Calculating particle time
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1

           prt_particles(ivac)%tsrc = tsp_t+r1*tsp_dt

           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = exsumg(ir)/real(gas_nvolex(ir))
           gas_eext=gas_eext+Ep0
           if (((gas_sig(ir)+gas_cap(iig,ir))*gas_drarr(ir)* &
                help < prt_tauddmc*gas_curvcent(ir)) &
                .or.(in_puretran)) then
              if(gas_isvelocity) then
                 prt_particles(ivac)%esrc = Ep0*(1.0+r0*mu0/pc_c)
                 prt_particles(ivac)%ebirth = Ep0*(1.0+r0*mu0/pc_c)
!-- velocity effects accounting
                 gas_evelo=gas_evelo-Ep0*r0*mu0/pc_c
!
              !(rev 120)
                 prt_particles(ivac)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
              !
                 prt_particles(ivac)%musrc = (mu0+r0/pc_c)/(1.0+r0*mu0/pc_c)
              else
                 prt_particles(ivac)%esrc = Ep0
                 prt_particles(ivac)%ebirth = Ep0
              !
                 prt_particles(ivac)%wlsrc = wl0
              !
                 prt_particles(ivac)%musrc = mu0
              endif
              prt_particles(ivac)%rtsrc = 1
           else
              prt_particles(ivac)%esrc = Ep0
              prt_particles(ivac)%ebirth = Ep0
              !(rev 120)
              prt_particles(ivac)%wlsrc = wl0
              !
              prt_particles(ivac)%musrc = mu0
              prt_particles(ivac)%rtsrc = 2
           endif
           !Setting ir = zone of particle
           prt_particles(ivac)%zsrc = ir
           !Setting particle index to not vacant
           prt_particles(ivac)%isvacant = .false.

           isnotvacnt = .true.
           
           !source tally
           !if(prt_particles(ivac)%rtsrc==2) then
           !   gas_eraddens(iig,ir)=gas_eraddens(iig,ir)+Ep0
           !endif
        else
           ir = ir+1
        endif
     enddo
  enddo
  

!-- Thermal volume particle instantiation: loop
  ir = 1
  irused(1:gas_nr) = 0

  do ipart = prt_nsurf+prt_nexsrc+1, prt_nnew
     ivac = prt_vacantarr(ipart)
     isnotvacnt = .false.
     !If adding particle ivac in current cell ir does not exceed nvol, add ivac to ir: loop
     do while (.not.isnotvacnt)
        if (irused(ir)>=gas_nvol(ir)) then
!-- do what?
           ir = ir + 1
        else
!-- or what?
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0d0
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           
           do ig = 1, gas_ng
              iig = ig
              if (r1>=denom2.and.r1<denom2+gas_emitprob(ig,ir)) exit
              denom2 = denom2+gas_emitprob(ig,ir)
           enddo
           !write(*,*) 'here',ivac
           !Ryan W.: particle group removed (rev. 120)
           !prt_particles(ivac)%gsrc = iig
           !Calculating wavelength uniformly from group
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           wl0 = 1d0/((1d0-r1)/gas_wl(iig)+r1/gas_wl(iig+1))
           !wl0 = 0.5d0*(gas_wl(iig)+gas_wl(iig+1))
!            x1 = pc_h*pc_c/(gas_wl(iig+1)*pc_kb*gas_temp(ir))
!            x2 = pc_h*pc_c/(gas_wl(iig)*pc_kb*gas_temp(ir))
!            if (x2<pc_plkpk) then
!               bmax = x2**3/(exp(x2)-1d0)
!            elseif (x1>pc_plkpk) then
!               bmax = x1**3/(exp(x1)-1d0)
!            else
!               bmax = pc_plkpk
!            endif
!            r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!            r2 = rand()
!           prt_tlyrand = prt_tlyrand+1
!            xx0 = (1d0-r1)*x1+r1*x2
!            do while (r2>xx0**3/(exp(xx0)-1d0)/bmax)
!               r1 = rand()
!           prt_tlyrand = prt_tlyrand+1
!               r2 = rand()
!           prt_tlyrand = prt_tlyrand+1
!               xx0 = (1d0-r1)*x1+r1*x2
!            enddo
!            wl0 = pc_h*pc_c/(xx0*pc_kb*gas_temp(ir))

           !Calculating radial position
           r1 = 0d0
           r2 = 1d0
           irl = max(ir-1,1)  !-- left neighbor
           irr = min(ir+1,gas_nr)  !-- right neighbor
           uul = .5d0*(gas_temp(irl)**4 + gas_temp(ir)**4)
           uur = .5d0*(gas_temp(irr)**4 + gas_temp(ir)**4)
           uumax = max(uul,uur)
           do while (r2 > r1)
              r3 = rand()
              prt_tlyrand = prt_tlyrand+1
              r0 = (r3*gas_rarr(ir+1)**3+(1.0-r3)*gas_rarr(ir)**3)**(1.0/3.0)
              r3 = (r0-gas_rarr(ir))/gas_drarr(ir)
              r1 = (r3*uur+(1.0-r3)*uul)/uumax
              r2 = rand()
              prt_tlyrand = prt_tlyrand+1
           enddo
           prt_particles(ivac)%rsrc = r0

           !Calculating direction cosine (comoving)
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           mu0 = 1d0-2d0*r1
           if(abs(mu0)<0.0000001d0) then
              mu0=0.0000001d0
           endif
           !Calculating particle time
           r1 = rand()
           prt_tlyrand = prt_tlyrand+1
           prt_particles(ivac)%tsrc = tsp_t+r1*tsp_dt
           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = gas_emit(ir)/real(gas_nvol(ir))

           if (((gas_cap(iig,ir)+gas_sig(ir))*gas_drarr(ir)* &
                help < prt_tauddmc*gas_curvcent(ir)) &
                .or.(in_puretran)) then
              if(gas_isvelocity) then
                 prt_particles(ivac)%esrc = Ep0*(1.0+r0*mu0/pc_c)
                 prt_particles(ivac)%ebirth = Ep0*(1.0+r0*mu0/pc_c)
!-- velocity effects accounting
                 gas_evelo=gas_evelo-Ep0*r0*mu0/pc_c
!
              !(rev 120)
                 prt_particles(ivac)%wlsrc = wl0/(1.0+r0*mu0/pc_c)
              !
                 prt_particles(ivac)%musrc = (mu0+r0/pc_c)/(1.0+r0*mu0/pc_c)
              else
                 prt_particles(ivac)%esrc = Ep0
                 prt_particles(ivac)%ebirth = Ep0
              !
                 prt_particles(ivac)%wlsrc = wl0
              !
                 prt_particles(ivac)%musrc = mu0
              endif
              prt_particles(ivac)%rtsrc = 1
           else
              prt_particles(ivac)%esrc = Ep0
              prt_particles(ivac)%ebirth = Ep0
              !(rev 120)
              prt_particles(ivac)%wlsrc = wl0
              !
              prt_particles(ivac)%musrc = mu0
              prt_particles(ivac)%rtsrc = 2
           endif
           !Setting ir = zone of particle
           prt_particles(ivac)%zsrc = ir
           !Setting particle index to not vacant
           prt_particles(ivac)%isvacant = .false.

           isnotvacnt = .true.

        endif
     enddo
  enddo


end subroutine interior_source
