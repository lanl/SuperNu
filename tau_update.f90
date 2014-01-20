subroutine tau_update

  use inputparmod
  use particlemod
  use timestepmod
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
!-- linear increase in thresholds
!      slp1=5d0*in_tauddmc/(tsp_nt*tsp_dt)
!      slp2=5d0*in_taulump/(tsp_nt*tsp_dt)

!      prt_tauddmc = in_tauddmc+(tsp_it-1)*tsp_dt*slp1
!      prt_taulump = in_taulump+(tsp_it-1)*tsp_dt*slp2



     if(tsp_it<=20) then
        return
     else
        prt_tauddmc=15d0
        prt_taulump=15d0
     endif

  else
     stop 'tau_update: prt_tauvtime invalid'
  endif


end subroutine tau_update
