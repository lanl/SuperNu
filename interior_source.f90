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

  integer :: ir, ipart, ivac, ig, iig
  integer, dimension(gas_nr) :: irused
  real*8 :: r1, r2, r3, uul, uur, uumax, mu0, r0, Ep0
  real*8 :: denom2
  logical :: isnotvacnt !checks for available particle space to populate in cell

  ir = 1
  irused(1:gas_nr) = 0
  !Volume particle instantiation: loop
  !Loop run over the number of new particles that aren't surface source
  !particles.
  do ipart = prt_nsurf+1, prt_nnew
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
           prt_particles(ivac)%gsrc = iig
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
              prt_particles(ivac)%musrc = (mu0+gas_velyes*r0/pc_c)/(1.0+gas_velyes*r0*mu0/pc_c)
              prt_particles(ivac)%rtsrc = 1
           else
              prt_particles(ivac)%Esrc = Ep0
              prt_particles(ivac)%Ebirth = Ep0
              prt_particles(ivac)%musrc = mu0
              prt_particles(ivac)%rtsrc = 2
           endif
           !Setting ir = zone of particle
           prt_particles(ivac)%zsrc = ir
           !Setting particle index to not vacant
           prt_particles(ivac)%isvacant = .false.

           isnotvacnt = .true.

        else
           ir = ir + 1
        endif
     enddo
  enddo

  !deallocate(prt_vacantarr)

end subroutine interior_source
