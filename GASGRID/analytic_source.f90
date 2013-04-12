subroutine analytic_source

  use gasgridmod
  use physconstmod
  use timestepmod
  use particlemod
  use inputparmod
  implicit none

  integer :: ir, ig
  real*8 :: x1, x2, x3, x4, srcren
  real*8 :: specint
  !
  real*8 :: uudd = 2.5d8  !velocity gauss width
  real*8 :: tmpgauss(gas_nr)  !manufactured gaussian temperature
  real*8 :: eradmanu(gas_ng,gas_nr) !manufactured rad. en. density
  
  real*8 :: aa11=1.371d14
  real*8 :: aa22=1.371d11
  real*8 :: rrcenter, aa00, bspeced, xx3, xx4

  x1 = 1d0/gas_wl(gas_ng+1)
  x2 = 1d0/gas_wl(1)
  if(gas_srctype=='heav') then
     !Heaviside source (uniform source sphere)
     if (tsp_time<=gas_theav*pc_day) then
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              x3 = 1d0/gas_wl(ig+1) 
              x4 = 1d0/gas_wl(ig)
              gas_exsource(ig,ir)=gas_srcmax*(x4-x3)/(x2-x1)
           enddo
        enddo
        do ir = gas_nheav+1, gas_nr
           do ig = 1, gas_ng
              gas_exsource(ig,ir)=0d0
           enddo
        enddo
     else
        do ir = 1, min(gas_nheav,gas_nr)
           do ig = 1, gas_ng
              gas_exsource(ig,ir)=0d0
           enddo
        enddo
        !Ryan W.: Temporary fix (?) for over generation,
        !resetting prt_ns
        !if(gas_sigcoef==0d0) then
        !   prt_ns = 1
        !endif
     endif
     !
  elseif(gas_srctype=='strt') then
     !Linear source profile
     do ir = 1, gas_nr
        srcren = gas_srcmax*(gas_rarr(gas_nr+1)- &
             0.5d0*(gas_rarr(ir)+gas_rarr(ir+1)))/ & 
             (gas_rarr(gas_nr+1)-gas_rarr(1))
        do ig = 1, gas_ng
           x3 = 1d0/gas_wl(ig+1) 
           x4 = 1d0/gas_wl(ig)
           gas_exsource(ig,ir)=srcren*(x4-x3)/(x2-x1)
        enddo
     enddo
  elseif(gas_srctype=='manu') then
     !
     !Manufactured Source (for gas_grptype='line')
     if(gas_grptype.ne.'line') then
        stop 'analytic_source: gas_grptype=line for gas_srctype=manu'
     endif
     !

     if(gas_isvelocity) then
        !calculating false gaussian temperature
        do ir = 1, gas_nr
           !cell center
           rrcenter=(gas_rarr(ir+1)+gas_rarr(ir))/2d0
           !
           tmpgauss(ir)=in_templ0*exp(-0.5d0*(rrcenter/uudd)**2)* &
                in_tfirst*pc_day/tsp_texp
           !
           !Thin lines
           do ig = 1, gas_ng, 2
              x3 = 1d0/gas_wl(ig+1)
              x4 = 1d0/gas_wl(ig)
              !calculating false thin group rad. energy density
              eradmanu(ig,ir)=((aa11*(in_velout-rrcenter)+ &
                   aa22*rrcenter)/in_velout)* &
                   (in_tfirst*pc_day/tsp_texp)**4 &
                   *(x4-x3)/(x2-x1)
              !
              !calculating manufactured source
              xx3 = x3*pc_h*pc_c/(1d3*pc_ev*tmpgauss(ir))
              xx4 = x4*pc_h*pc_c/(1d3*pc_ev*tmpgauss(ir))
              aa00=in_tfirst*pc_day*in_velout
              bspeced = 15d0*specint(xx3,xx4,3)/pc_pi**4
              !
              gas_exsource(ig,ir) = &
                   (aa00*aa11/(pc_day*in_tfirst*in_velout)**5 - &
                   5.0*eradmanu(ig,ir)/(tsp_texp*in_velout))*in_velout &
                   +2.0*pc_c*eradmanu(ig,ir)/rrcenter + &
                   pc_c*(aa22-aa11)*aa00/(tsp_texp*in_velout)**5 - &
                   pc_c*bspeced*gas_sigmapg(ig,ir)*pc_acoef*tmpgauss(ir)**4 &
                   -gas_sigmapg(ig,ir)*rrcenter*(pc_c-4d0*rrcenter/3d0)* &
                   eradmanu(ig,ir)/pc_c
              !
              !write(*,*) 'thin: ',gas_exsource(ig,ir)
              !if(gas_exsource(ig,ir)<0d0) then
              !   gas_exsource(ig,ir)=0d0
              !endif
              !
           enddo
           !
           !Thick lines
           do ig = 2, gas_ng, 2
              x3 = 1d0/gas_wl(ig+1)
              x4 = 1d0/gas_wl(ig)
              !calculating false thick group rad. energy density
              eradmanu(ig,ir)=pc_acoef*tmpgauss(ir)**4*(x4-x3)/(x2-x1)
              !
              !calculating manufactured source
              xx3 = x3*pc_h*pc_c/(1d3*pc_ev*tmpgauss(ir))
              xx4 = x4*pc_h*pc_c/(1d3*pc_ev*tmpgauss(ir))
              bspeced = 15d0*specint(xx3,xx4,3)/pc_pi**4
              !
              gas_exsource(ig,ir) = eradmanu(ig,ir)* &
                   (4d0*((rrcenter*tsp_texp)**2*uudd/ &
                   (uudd*tsp_texp)**3-1d0/tsp_texp &
                   +pc_c*(3d0*(uudd*tsp_texp)**2- &
                   4.0*(rrcenter*tsp_texp)**2) &
                   /(3d0*gas_sigmapg(ig,ir)*(uudd*tsp_texp)**4)) &
                   -pc_c*bspeced*gas_sigmapg(ig,ir)*(x2-x1)/(x4-x3) &
                   +pc_c*gas_sigmapg(ig,ir) &
                   +4.0*gas_sigmapg(ig,ir)*rrcenter*(rrcenter &
                   -pc_c*rrcenter/(gas_sigmapg(ig,ir)*tsp_texp*uudd**2)))
              !
              !write(*,*) 'thick: ',gas_exsource(ig,ir)
              !if(gas_exsource(ig,ir)<0d0) then
              !   gas_exsource(ig,ir)=0d0
              !endif
              !
           enddo
           !
        enddo
     else
        stop 'analytic_source: no static manufactured source'
     endif
  else
     stop 'analytic_source: gas_srctype invalid'
  endif
  !write(*,*) gas_exsource(1,1), gas_exsource(2,1), gas_exsource(3,1)
  !write(*,*) gas_vals2(1)%tempkev, gas_vals2(2)%tempkev
  !write(*,*) gas_exsource(1,:)
  !write(*,*)
  !write(*,*) gas_exsource(2,:)

end subroutine analytic_source
