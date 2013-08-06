subroutine leakage_opacity

  use gasgridmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC grouped boundary probabilities (albedo) and leakage opacities.
!##################################################

  integer :: ir, ig
  real*8 :: tt, gg, ggg, eps, bb, sigtot
  real*8 :: curvleft, curvright, help

  logical :: missive = .false.
  ! Here: left=>toward r=0 and right=>outward

  
  !Calculating IMC-to-DDMC leakage albedo coefficients (Densmore, 2007): loop
  !These quantities may not need to be stored directly (pending further analysis)
  if(missive) then
     do ir = 1, gas_nr
        !calculating left cell curvature
        curvleft = gas_rarr(ir)**2/ &
             (gas_rarr(ir+1)**2+gas_rarr(ir)*gas_rarr(ir+1)+gas_rarr(ir)**2)
        !calculating right cell curvature
        curvright = gas_rarr(ir+1)**2/ &
             (gas_rarr(ir+1)**2+gas_rarr(ir)*gas_rarr(ir+1)+gas_rarr(ir)**2)
        !write(*,*) curvleft, curvright
        do ig = 1, gas_ng
           !calculating left albedo
           sigtot=gas_caprosl(ig,ir)+gas_sigbl(ir)
           gg = sqrt(3d0*gas_fcoef(ir)*gas_caprosl(ig,ir)/sigtot)
           !calculating left optical depth
           tt = sigtot*gas_drarr(ir)
           if(gas_isvelocity) then
              tt = tt*tsp_texp
           endif
           !calculating left discretization eigenvalue
           ggg=0.5d0*(curvleft+curvright)/curvright+(gg*tt)**2/(6d0*curvright) &
                -sqrt((0.5d0*(curvright-curvleft)/curvright)**2 &
                +(curvleft+curvright)*(gg*tt/curvright)**2/6d0 &
                +(gg*tt)**4/(36d0*curvright**2))
           !write(*,*) 'left:', ggg, 0.5d0*(curvleft+curvright)/curvright !+(gg*tt)**2/(6d0*curvright)
           !calculating left conveniency coefficient
           bb = curvright*(1d0-ggg)/curvleft+(gg*tt)**2/(3.0*curvleft)
           !calculating left emissivity
           help = sigtot*gas_rarr(ir)
           if(gas_isvelocity) then
              help = help*tsp_texp
           endif
           eps = (4.0/3.0)*(gg+1d0/help) &
                /(1d0+0.7104d0*(gg+1d0/help))
           !write(*,*) 'here',eps
           if(eps>0d0.and.eps<=1.0001d0) then
              gas_ppl(ig,ir) = 0.5*eps*bb/(bb-3d0*eps*tt/4d0)
              if(gas_ppl(ig,ir)<0d0) then
                 gas_ppl(ig,ir) = 4d0/(3d0*tt+6d0*0.7104d0)
              endif
           else
              gas_ppl(ig,ir) = 4d0/(3d0*tt+6d0*0.7104d0)
           endif
           !
           !calculating right albedo
           sigtot=gas_caprosr(ig,ir)+gas_sigbr(ir)
           gg = sqrt(3d0*gas_fcoef(ir)*gas_caprosr(ig,ir)/sigtot)
           !calculating right optical depth
           tt = sigtot*gas_drarr(ir)
           if(gas_isvelocity) then
              tt = tt*tsp_texp
           endif
           !calculating right discretization eigenvalue
           ggg=0.5d0*(curvleft+curvright)/curvright+(gg*tt)**2/(6.0*curvright) &
                +sqrt((0.5*(curvright-curvleft)/curvright)**2 &
                +(curvleft+curvright)*(gg*tt/curvright)**2/6.0 &
                +(gg*tt)**4/(36.0*curvright**2))
           !write(*,*) 'right: ', ggg
           !calculating right conveniency coefficient
           bb = curvleft*(1d0-1d0/ggg)/curvright+(gg*tt)**2/(3.0*curvright)
           !calculating right emissivity
           help = sigtot*gas_rarr(ir+1)
           if(gas_isvelocity) then
              help = help*tsp_texp
           endif
           eps = (4.0/3.0)*(gg-1d0/help) &
                /(1d0+0.7104d0*(gg-1d0/help))
           !write(*,*) 'here', eps
           if(eps>0d0.and.eps<=1.0001d0) then
              !write(*,*) 'here'
              gas_ppr(ig,ir) = 0.5*eps*bb/(bb-3d0*eps*tt/4d0)
              if(gas_ppr(ig,ir)<0d0) then
                 gas_ppr(ig,ir) = 4d0/(3d0*tt+6d0*0.7104d0)
              endif
           else
              gas_ppr(ig,ir) = 4d0/(3d0*tt+6d0*0.7104d0)
           endif
        enddo
     enddo
  else
     do ir = 1, gas_nr
        do ig = 1, gas_ng
           tt = (gas_caprosl(ig,ir)+gas_sigbl(ir)) &
                *gas_drarr(ir)
           if(gas_isvelocity) then
              tt = tt*tsp_texp
           endif
           gas_ppl(ig,ir) = 4.0d0/(3d0*tt+6d0*0.7104d0)
           !
           tt = (gas_caprosr(ig,ir)+gas_sigbr(ir)) &
                *gas_drarr(ir)
           if(gas_isvelocity) then
              tt = tt*tsp_texp
           endif
           gas_ppr(ig,ir) = 4.0d0/(3d0*tt+6d0*0.7104d0)
           !
        enddo
     enddo
  endif


  !Calculating DDMC(-to-IMC) leakage opacities (Densmore, 2007, 2012): loop
  do ir = 1, gas_nr
     do ig = 1, gas_ng
        !Computing left-leakage opacities
        
        if (ir==1) then
        !
           gas_opacleakl(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
           gas_opacleakl(ig,ir)=gas_opacleakl(ig,ir)/ &
                (3d0*gas_vals2(ir)%vol/pc_pi4)
           if(gas_isvelocity) then
              gas_opacleakl(ig,ir) = gas_opacleakl(ig,ir)*tsp_texp**2
           endif
        !   
        else
           help = (gas_sig(ir-1)+gas_cap(ig,ir-1))*gas_drarr(ir-1)
           if(gas_isvelocity) then
              help = help*tsp_texp
           endif
           if(help < prt_tauddmc*gas_curvcent(ir-1)) then
        !   
              gas_opacleakl(ig,ir)=1.5*gas_ppl(ig,ir)*gas_rarr(ir)**2
              gas_opacleakl(ig,ir)=gas_opacleakl(ig,ir)/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
              if(gas_isvelocity) then
                 gas_opacleakl(ig,ir) = gas_opacleakl(ig,ir)*tsp_texp**2
              endif
        !   
           else
        !
              tt = (gas_sigbl(ir)+gas_caprosl(ig,ir))*gas_drarr(ir)+ &
                   (gas_sigbr(ir-1)+gas_caprosr(ig,ir-1))*gas_drarr(ir-1)
        !   
              gas_opacleakl(ig,ir)=(2.0*gas_rarr(ir)**2)/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
              
              gas_opacleakl(ig,ir) = gas_opacleakl(ig,ir)/tt
              
              if(gas_isvelocity) then
                 gas_opacleakl(ig,ir) = gas_opacleakl(ig,ir)*tsp_texp
              endif
           endif

        endif
        !Computing right-leakage opacities
        if (ir==gas_nr) then
        !
           gas_opacleakr(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
           gas_opacleakr(ig,ir)=gas_opacleakr(ig,ir)/ &
                (3d0*gas_vals2(ir)%vol/pc_pi4)
           if(gas_isvelocity) then
              gas_opacleakr(ig,ir) = gas_opacleakr(ig,ir)*tsp_texp**2
           endif
        !   
        else
           help = (gas_sig(ir+1)+gas_cap(ig,ir+1))*gas_drarr(ir+1)
           if(gas_isvelocity) then
              help = help*tsp_texp
           endif
           if(help < prt_tauddmc*gas_curvcent(ir+1)) then
        !   
              gas_opacleakr(ig,ir)=1.5*gas_ppr(ig,ir)*gas_rarr(ir+1)**2
              gas_opacleakr(ig,ir)=gas_opacleakr(ig,ir)/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
              if(gas_isvelocity) then
                 gas_opacleakr(ig,ir) = gas_opacleakr(ig,ir)*tsp_texp**2
              endif
        !   
           else
        !
              tt = (gas_sigbr(ir)+gas_caprosr(ig,ir))*gas_drarr(ir)+ &
                   (gas_sigbl(ir+1)+gas_caprosl(ig,ir+1))*gas_drarr(ir+1)
              gas_opacleakr(ig,ir) = (2.0*gas_rarr(ir+1)**2)/ &
                   (3d0*gas_vals2(ir)%vol/pc_pi4)
              gas_opacleakr(ig,ir) = gas_opacleakr(ig,ir)/tt
              if(gas_isvelocity) then
                 gas_opacleakr(ig,ir) = gas_opacleakr(ig,ir)*tsp_texp
              endif
           endif
        endif
     enddo
  enddo
  

end subroutine leakage_opacity
