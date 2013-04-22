subroutine interior_source

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  use inputparmod

  implicit none

!##################################################
  !This subroutine instantiates new volume (cell) particle properties.
!##################################################

  integer :: ir,iir, ipart, ivac, ig, iig
  integer, dimension(gas_nr) :: irused
  real*8 :: r1, r2, r3, uul, uur, uumax, mu0, r0, Ep0, wl0
  real*8 :: denom2
  real*8, dimension(gas_nr) :: exsumg
  logical :: isnotvacnt !checks for available particle space to populate in cell

  ir = 1
  irused(1:gas_nr) = 0
  exsumg(1:gas_nr) = 0d0
  !Volume particle instantiation: loop
  !Loop run over the number of new particles that aren't surface source
  !particles.
  do iir = 1, gas_nr
     do ig = 1, gas_ng
        exsumg(iir)=exsumg(iir)+gas_exsource(ig,iir)
     enddo
     !write(*,*) exsumg(ir)
  enddo
  do ipart = prt_nsurf+1, prt_nsurf+prt_nexsrc
     ivac = prt_vacantarr(ipart)
     isnotvacnt = .false.
     !If adding particle ivac in current cell ir does not exceed gas_vals2, add ivac to ir: loop
     do while (isnotvacnt.eqv..false.)
        if (irused(ir)<gas_vals2(ir)%nvolex) then
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0d0
           r1 = rand()
           do ig = 1, gas_ng
              iig = ig
              if(r1>=denom2.and.r1<denom2+gas_exsource(ig,ir)/ &
                   exsumg(ir)) exit
              denom2 = denom2+gas_exsource(ig,ir)/exsumg(ir)
           enddo
           !Ryan W.: particle group removed (rev. 120)
           !prt_particles(ivac)%gsrc = iig
           !Calculating comoving wavelength uniformly from group
           r1 = rand()
           wl0 = (1d0-r1)*gas_wl(iig)+r1*gas_wl(iig+1)
           !write(*,*) wl0, gas_wl(iig)+r1*gas_wl(iig+1)
           !write(*,*) gas_wl
           !Calculating radial position
           r3 = rand()
           prt_particles(ivac)%rsrc = (r3*gas_rarr(ir+1)**3 + &
                (1.0-r3)*gas_rarr(ir)**3)**(1.0/3.0)
           r0 = prt_particles(ivac)%rsrc
           !Calculating direction cosine (comoving)
           if(gas_isanalsrc.and.gas_srctype=='manu'.and.mod(iig,2)==0) then
              mu0=1d0
           else
              r1 = rand()
              mu0 = 1d0-2d0*r1
           endif
           !Calculating particle tsp_time
           r1 = rand()
           prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt
           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = exsumg(ir)*tsp_dt* &
             (4.0*pc_pi*gas_vals2(ir)%dr3_34pi/3.0)/real(gas_vals2(ir)%nvolex)
           !write(*,*) Ep0, gas_vals2(ir)%nvolex
           if (((gas_sig(ir)+gas_sigmapg(iig,ir))*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0) &
                .or.(in_puretran.eqv..true.)) then
              prt_particles(ivac)%Esrc = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%Ebirth = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              !(rev 120)
              prt_particles(ivac)%wlsrc = wl0/(1.0+gas_velyes*r0*mu0/pc_c)
              !
              prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%rtsrc = 1
           else
              prt_particles(ivac)%Esrc = Ep0
              prt_particles(ivac)%Ebirth = Ep0
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
           !   gas_eraddensg(iig,ir)=gas_eraddensg(iig,ir)+Ep0
           !endif
        else
           ir = ir+1
        endif
     enddo
  enddo

  !Thermal volume particle instantiation: loop
  ir = 1
  irused(1:gas_nr) = 0
  do ipart = prt_nsurf+prt_nexsrc+1, prt_nnew
     ivac = prt_vacantarr(ipart)
     isnotvacnt = .false.
     !If adding particle ivac in current cell ir does not exceed gas_vals2, add ivac to ir: loop
     do while (isnotvacnt.eqv..false.)
        if (irused(ir)<gas_vals2(ir)%nvol) then
           irused(ir) = irused(ir)+1
           !Calculating Group
           denom2 = 0d0
           r1 = rand()
           do ig = 1, gas_ng
              iig = ig
              if (r1>=denom2.and.r1<denom2+gas_emitprobg(ig,ir)) exit
              denom2 = denom2+gas_emitprobg(ig,ir)
           enddo
           !Ryan W.: particle group removed (rev. 120)
           !prt_particles(ivac)%gsrc = iig
           !Calculating wavelength uniformly from group
           r1 = rand()
           wl0 = (1d0-r1)*gas_wl(iig)+r1*gas_wl(iig+1)
           !write(*,*) wl0, iig
           !write(*,*) gas_wl
           !Calculating radial position
           r1 = 0d0
           r2 = 1d0
           uul = gas_tempb(ir)**4
           uur = gas_tempb(ir+1)**4
           uumax = max(uul,uur)
           do while (r2 > r1)
              r3 = rand()
              r0 = (r3*gas_rarr(ir+1)**3+(1.0-r3)*gas_rarr(ir)**3)**(1.0/3.0)
              r3 = (r0-gas_rarr(ir))/gas_drarr(ir)
              r1 = (r3*uur+(1.0-r3)*uul)/uumax
              r2 = rand()
           enddo
           prt_particles(ivac)%rsrc = r0

           !Calculating direction cosine (comoving)
           r1 = rand()
           mu0 = 1d0-2d0*r1
           !Calculating particle tsp_time
           r1 = rand()
           prt_particles(ivac)%tsrc = tsp_time+r1*tsp_dt
           !Calculating particle energy, lab frame direction and propagation type
           Ep0 = gas_vals2(ir)%emit/real(gas_vals2(ir)%nvol)
           if ((gas_sigmapg(iig,ir)*gas_drarr(ir)*(gas_velno*1.0+gas_velyes*tsp_texp)<5.0d0).OR.(in_puretran.eqv..true.)) then
              prt_particles(ivac)%Esrc = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%Ebirth = Ep0*(1.0+gas_velyes*r0*mu0/pc_c)
              !(rev 120)
              prt_particles(ivac)%wlsrc = wl0/(1.0+gas_velyes*r0*mu0/pc_c)
              !
              prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%rtsrc = 1
           else
              prt_particles(ivac)%Esrc = Ep0
              prt_particles(ivac)%Ebirth = Ep0
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
           !   gas_eraddensg(iig,ir)=gas_eraddensg(iig,ir)+Ep0
           !endif

        else
           ir = ir + 1
        endif
     enddo
  enddo

  !deallocate(prt_vacantarr)

end subroutine interior_source
