subroutine tau_update(t,tfirst,tlast)

  use physconstmod
  use transportmod
  implicit none
  real*8,intent(in) :: t,tfirst,tlast
!-------------------------------------------------
! Update mfp thresholds for DDMC:
! - IMC-DDMC transition,
! - DDMC group lumping.
!-------------------------------------------------
  real*8 :: slp1,slp2

  if(trn_tauvtime=='unif') then
!-- constant thresholds
     return
  elseif(trn_tauvtime=='incr') then
!-- linear increase in thresholds (max mfp thresh = 5)
     slp1=0.667d0*trn_tauddmc/(tlast-tfirst)
     slp2=0.667d0*trn_taulump/(tlast-tfirst)

     trn_tauddmc = trn_tauddmc+(t-tfirst)*slp1
     trn_taulump = trn_taulump+(t-tfirst)*slp2
  else
     stop 'tau_update: trn_tauvtime invalid'
  endif

end subroutine tau_update
