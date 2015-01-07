subroutine emission_probability(icell1,ncell)

  use miscmod
  use inputparmod
  use timingmod
  use gridmod
  use groupmod
  use physconstmod
  implicit none
  integer,intent(in) :: icell1, ncell

!-----------------------
!multigroup volume emission probabilities
!-----------------------

  integer :: l,ig,igp1,iep,nepg
  real*8 :: t0,t1
  real*8 :: help
  real*8 :: specarr(grp_ng)

  t0 = t_time()

!-- grouped volume emission probabilities:
  if(in_opacanaltype=='pick') then
     stop 'emission_probability: not implemented'
     do l=icell1,icell1+ncell-1
        grd_emitprob(1,l) = in_suolpick1*grd_cap(1,l)
        grd_emitprob(2,l) = (1d0 - in_suolpick1)*grd_cap(2,l)
        grd_emitprob(3:grp_ng,l) = 0d0  !-- not necessary
     enddo !l
     return
  endif

!-- one group
  if(grp_ng==1) then
     grd_emitprob(:,icell1:icell1+ncell-1) = 1d0
     return
  endif

!-- multi-group
  do l=icell1,icell1+ncell-1
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
  enddo !l

  t1 = t_time()
  call timereg(t_emitp,t1-t0)

end subroutine emission_probability
