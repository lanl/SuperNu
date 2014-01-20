subroutine tau_update

  use inputparmod
  use particlemod
  use timestepmod
  use physconstmod
  implicit none
!-------------------------------------------------
! Update mfp thresholds for DDMC:
! - IMC-DDMC transition,
! - DDMC group lumping.
!-------------------------------------------------
  real*8 :: slp1,slp2

  if(prt_tauvtime=='unif') then
!-- constant thresholds
     return
  elseif(prt_tauvtime=='incr') then
!-- linear increase in thresholds (max mfp thresh = 5)
     slp1=0.667d0*in_tauddmc/(pc_day*(in_tlast-in_tfirst))
     slp2=0.667d0*in_taulump/(pc_day*(in_tlast-in_tfirst))

     prt_tauddmc = in_tauddmc+(tsp_t-pc_day*in_tfirst)*slp1
     prt_taulump = in_taulump+(tsp_t-pc_day*in_tfirst)*slp2


!      if(tsp_it<=20) then
!         return
!      else
!         prt_tauddmc=15d0
!         prt_taulump=15d0
!      endif

  else
     stop 'tau_update: prt_tauvtime invalid'
  endif


end subroutine tau_update
