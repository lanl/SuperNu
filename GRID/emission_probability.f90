subroutine emission_probability

  use miscmod
  use inputparmod
  use timingmod
  use gridmod
  use groupmod
  use physconstmod
  implicit none

!-----------------------
!multigroup volume emission probabilities
!-----------------------

  integer :: i,j,k,l,ig,igp1,iep,nepg
  real*8 :: t0,t1
  real*8 :: help
  real*8 :: specarr(grp_ng)

  t0 = t_time()

!-- grouped volume emission probabilities:
  if(in_opacanaltype=='pick') then
     stop 'emission_probability: not implemented'
     do k=1,grd_nz
     do j=1,grd_ny
     do i=1,grd_nx
        l = grd_icell(i,j,k)
        grd_emitprob(1,l) = in_suolpick1*grd_cap(1,l)
        grd_emitprob(2,l) = (1d0 - in_suolpick1)*grd_cap(2,l)
        grd_emitprob(3:grp_ng,l) = 0d0  !-- not necessary
     enddo !i
     enddo !j
     enddo !k
     return
  endif

!-- one group
  if(grp_ng==1) then
     grd_emitprob = 1d0
     return
  endif

!-- multi-group
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
     l = grd_icell(i,j,k)
!-- piecewise integration of planck function
     specarr = specintv(1d0/grd_temp(l),0)
!-- cumulative sum of unnormalized emission probability
     ig = 1
     help = 0d0
     do iep=1,grd_nep
        nepg = min(iep*grd_nepg,grp_ng) - ig + 1
        igp1 = ig + nepg - 1
        help = help + sum(specarr(ig:igp1)*grd_cap(ig:igp1,l))
        grd_emitprob(iep,l) = help
        ig = igp1 + 1
     enddo !iep
  enddo !i
  enddo !j
  enddo !k

  t1 = t_time()
  call timereg(t_emitp,t1-t0)

end subroutine emission_probability
