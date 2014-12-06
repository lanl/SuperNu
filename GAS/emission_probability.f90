subroutine emission_probability

  use miscmod
  use inputparmod
  use timingmod
  use gridmod
  use physconstmod
  use miscmod, only:specint
  implicit none

!-----------------------
!multigroup volume emission probabilities
!-----------------------

  integer :: i,ig,igp1,j,k,iep,nepg
  real*8 :: t0,t1
  real*8 :: help
  real*8 :: specarr(grd_ng)

  call time(t0)

!-- grouped volume emission probabilities:
  if(in_opacanaltype=='pick') then
     stop 'emission_probability: not implemented'
     do k=1,grd_nz
     do j=1,grd_ny
     do i=1,grd_nx
        grd_emitprob(1,i,j,k) = in_suolpick1*grd_cap(1,i,j,k)
        grd_emitprob(2,i,j,k) = (1d0 - in_suolpick1)*grd_cap(2,i,j,k)
        grd_emitprob(3:grd_ng,i,j,k) = 0d0  !-- not necessary
     enddo !i
     enddo !j
     enddo !k
     return
  endif

!-- one group
  if(grd_ng==1) then
     grd_emitprob = 1d0
     return
  endif

!-- multi-group
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
!-- piecewise integration of planck function
     specarr = specint3(pc_h*pc_c*grd_wlinv/(pc_kb*grd_temp(i,j,k)),grd_ng+1,1)
!-- cumulative sum of unnormalized emission probability
     ig = 1
     help = 0d0
     do iep=1,grd_nep
        nepg = min(iep*grd_nepg,grd_ng) - ig + 1
        igp1 = ig + nepg - 1
        help = help + sum(specarr(ig:igp1)*grd_cap(ig:igp1,i,j,k))
        grd_emitprob(iep,i,j,k) = help
        ig = igp1 + 1
     enddo !iep
     if(ig/=grd_ng+1) write(0,*) ig,grd_ng
     if(ig/=grd_ng+1) stop 'emission-probability: ig/=ng'
  enddo !i
  enddo !j
  enddo !k

  call time(t1)
  call timereg(t_emitp,t1-t0)

end subroutine emission_probability
