subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use particlemod
  use inputparmod
  use manufacmod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2, x3, x4, srcren
  real*8 :: specint
  !
  real*8 :: uudd = 2.5d8  !velocity gauss width
  real*8 :: tmpgauss(gas_nr)  !manufactured gaussian temperature
  real*8 :: eradthin(gas_nr) !manufactured rad. en. density
  
  real*8 :: rrcenter, bspeced, xx3, xx4
  real*8 :: ddrr2, ddrr3, ddrr4, aleff1 = 1d0

  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)

  gas_emitex = 0d0

  if(gas_srctype=='none') then
    return
  elseif(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)!{{{
     if (tsp_time<=gas_theav*pc_day) then
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1) 
              x4 = 1d0/gas_wl(ig)
              gas_emitex(ig,ir)=gas_srcmax*(x4-x3)/(x2-x1)
              if(gas_isvelocity) then
                 gas_emitex(ig,ir)=gas_emitex(ig,ir)/tsp_texp**3
              endif
           enddo
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!
        enddo
        do ir = gas_nheav+1, gas_nr
           do ig = 1, gas_ng
              gas_emitex(ig,ir)=0d0
           enddo
        enddo
     else
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              gas_emitex(ig,ir)=0d0
           enddo
        enddo
     endif
     !!}}}
  elseif(gas_srctype=='strt') then
     !Linear source profile!{{{
     do ir = 1, gas_nr
        srcren = gas_srcmax*(gas_rarr(gas_nr+1)- &
             0.5d0*(gas_rarr(ir)+gas_rarr(ir+1)))/ & 
             (gas_rarr(gas_nr+1)-gas_rarr(1))
        do ig = 1, gas_ng
           x3 = 1d0/gas_wl(ig+1) 
           x4 = 1d0/gas_wl(ig)
           gas_emitex(ig,ir)=srcren*(x4-x3)/(x2-x1)
        enddo
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!
     enddo!}}}
  elseif(gas_srctype=='manu') then
     !!{{{
     !Manufactured Source (for gas_opacanaltype='line')
     !if(gas_opacanaltype.ne.'line') then
     !   stop 'analytic_source: gas_opacanaltype=line for gas_srctype=manu'
     !endif
     !

     if(gas_isvelocity) then
        !
        if(gas_ng>2) then
           stop 'analytic_source: gas_ng<=2 manufactured'
        endif
        do ir = 1, gas_nr           
           !Thin lines
           do ig = 1, gas_ng, 2
              x3 = 1d0/gas_wl(ig+1)
              x4 = 1d0/gas_wl(ig)
              !calculating manufactured source
              xx3 = x3*pc_h*pc_c/(pc_kb*man_temp0)
!              xx4 = x4*pc_h*pc_c/(pc_kb*man_temp0)
              
!              bspeced = 15d0*specint(xx3,xx4,3)/pc_pi**4
              !
!               gas_emitex(ig,ir)= (1d0/tsp_dt)*(&
!                    log((tsp_texp+tsp_dt)/tsp_texp)*(3d0*man_aa11/pc_c)+&
!                    (3d0*in_totmass*in_sigcoef/(8d0*pc_pi*gas_velout))* &
!                    ((gas_velout*tsp_texp)**(-2d0)-&
!                    (gas_velout*(tsp_texp+tsp_dt))**(-2d0))*&
!                    (man_aa11-pc_acoef*pc_c*man_temp0**4) &
!                    )*(x4-x3)/(x2-x1)
              gas_emitex(ig,ir)=(1d0/tsp_dt)*&
                   log((tsp_texp+tsp_dt)/tsp_texp)*(man_aa11/pc_c)*&
                   (1.5d0-(1d0-aleff1)*0.5d0*x3/(x4-x3))!-0.5*(15d0/pc_pi**4)*xx3**4/(exp(xx3)-1))
                   !2.0d0
              !write(*,*) x3/(x4-x3)
              !
              !gas_emitex(ig,ir)=0d0
           enddo
           !
           !Thick lines
            do ig = 2, gas_ng, 2
               x3 = 1d0/gas_wl(ig)
               x4 = 1d0/gas_wl(ig-1)
               xx3 = x3*pc_h*pc_c/(pc_kb*man_temp0)
               gas_emitex(ig,ir)=(1d0/tsp_dt)*&
                    log((tsp_texp+tsp_dt)/tsp_texp)*(man_aa11/pc_c)*&
                    (2d0-aleff1*0.5d0*x3/(x4-x3))
            enddo
           !
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!
        enddo
        !write(*,*) eradthin(1), eradthin(5), eradthin(gas_nr), gas_nr
     
     else
        !constant linear profile test (grey)
        if(gas_opacanaltype.ne.'grey') then
           stop 'analytic_source: grey for static manu'
        endif
        !
        do ir = 1, gas_nr
           ddrr3 = gas_rarr(ir+1)**3-gas_rarr(ir)**3
           ddrr4 = gas_rarr(ir+1)**4-gas_rarr(ir)**4           
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1)
              x4 = 1d0/gas_wl(ig)
              gas_emitex(ig,ir)=gas_cap(ig,ir)* &
                   (man_aa11-0.75d0*(man_aa11-man_aa22)* &
                   ddrr4/(gas_rarr(gas_nr+1)*ddrr3) - &
                   pc_c*pc_acoef*man_temp0**4)
              gas_emitex(ig,ir)=gas_emitex(ig,ir)*(x4-x3)/(x2-x1)
              !write(*,*) gas_emitex(ig,ir)
              !gas_emitex(ig,ir) = 0d0
           enddo
!
           gas_emitex(:,ir) = gas_emitex(:,ir)* &
                gas_vals2(ir)%vol*tsp_dt
!           
        enddo
        !
        !stop 'analytic_source: no static manufactured source'
     endif!}}}
  else
     stop 'analytic_source: gas_srctype invalid'
  endif

  !write(*,*) gas_emitex(1,1), gas_emitex(2,1), gas_emitex(3,1)  
  !write(*,*) gas_emitex(1,:)
  !write(*,*)
  !write(*,*) gas_emitex(2,:)
  !write(*,*) gas_siggrey(gas_nr), gas_cap(1,gas_nr)

end subroutine analytic_source
