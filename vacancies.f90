subroutine vacancies

  use particlemod
  implicit none

!##################################################
  !This subroutine creates an array of vacant index locations
  !in the particle array each time step.
!##################################################

  integer :: ipart, ivac

  !Initializing index counters and full array checking boolean
  ipart = 0
  ivac = 0

  !Filling prt_vacantarr with particle index of vacant particles: loop
  do ipart=1,prt_npartmax
     if (prt_particles(ipart)%isvacant) then
        ivac = ivac+1
        prt_vacantarr(ivac) = ipart
     endif
     if (ivac == prt_nnew) exit
  enddo
  if(ipart > prt_npartmax) write(6,*) 'Maximum number of prt_particles reached'

end subroutine vacancies
