subroutine vacancies

  use particlemod
  implicit none

!##################################################
  !This subroutine creates an array of vacant index locations
  !in the particle array each time step.
!##################################################

  integer :: ipart, ivac
  logical :: isfull

  !Initializing index counters and full array checking boolean
  isfull = .false.
  ipart = 0
  ivac = 0

  !Filling prt_vacantarr with particle index of vacant particles: loop
  do while (isfull.eqv..false.)
     ipart = ipart+1
     if (prt_particles(ipart)%isvacant.eqv..true.) then
        ivac = ivac+1
        prt_vacantarr(ivac) = ipart
     endif
     if (ivac == prt_nnew) then
        isfull = .true.
     elseif (ipart == prt_npartmax) then
        write(6,*) 'Maximum number of prt_particles reached'
        isfull = .true.
     endif
  enddo
  
end subroutine vacancies
