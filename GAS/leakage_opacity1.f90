subroutine leakage_opacity1

  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC 1D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: i,j,k, ig, khelp
  real*8 :: thelp, help
  real*8 :: speclump, specval
  real*8 :: specarr(grp_ng)
  real*8 :: pp
!-- statement functions
  integer :: l
  real*8 :: dx,dx2,dx3,xm,dy,dyac,ym,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dx2(l)= grd_xarr(l+1)**2 - grd_xarr(l)**2
  dx3(l)= grd_xarr(l+1)**3 - grd_xarr(l)**3
  xm(l) = 0.5*(grd_xarr(l+1) + grd_xarr(l))
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
  dyac(l) = grd_yacos(l) - grd_yacos(l+1)
  ym(l) = sqrt(1d0-0.25*(grd_yarr(l+1)+grd_yarr(l))**2)
  dz(l) = grd_zarr(l+1) - grd_zarr(l)
!
!-- setting vel-space helper
  if(grd_isvelocity) then
     thelp = tsp_t
  else
     thelp = 1d0
  endif

!-- init (necessary for domain decomposition)
  grd_opacleak = 0d0

!
!-- calculating leakage opacities
  do k=1,grd_nz
  do j=1,grd_ny
  do i=1,grd_nx
!
!-- initializing Planck integral vectorized
     specarr = specintv(1d0/grd_temp(i,j,k),0)
     help = min(dx(i),xm(i)*dyac(j),xm(i)*ym(j)*dz(k))*thelp
     speclump = 1d0/sum(specarr,grd_cap(:,i,j,k)*help>=prt_taulump)
!-- lumping opacity
     do ig=1,grp_ng
        if(grd_cap(ig,i,j,k)*min(dx(i),xm(i)*dyac(j),xm(i)*ym(j)*dz(k)) * &
           thelp < prt_taulump) cycle
!
!-- obtaining spectral weight
        specval = specarr(ig)
!
!-- calculating i->i-1 leakage opacity
        if(i==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i-1,j,k)+ &
              grd_sig(i-1,j,k))*min(dx(i-1),xm(i-1)*dyac(j), &
              xm(i-1)*ym(j)*dz(k))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                1.5d0*pp*grd_xarr(i)**2/(dx3(i)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i-1,j,k)+grd_cap(ig,i-1,j,k))*dx(i-1))*thelp
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                2.0d0*grd_xarr(i)**2/(dx3(i)*thelp*help)
        endif

!
!-- calculating i->i+1 leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i+1,j,k)+ &
              grd_sig(i+1,j,k))*min(dx(i+1),xm(i+1)*dyac(j), &
              xm(i+1)*ym(j)*dz(k))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                1.5d0*pp*grd_xarr(i+1)**2/(dx3(i)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i+1,j,k)+grd_cap(ig,i+1,j,k))*dx(i+1))*thelp
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                2.0d0*grd_xarr(i+1)**2/(dx3(i)*thelp*help)
        endif

!
!-- calculating j->j-1 leakage opacity
        if(j==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j-1,k)+ &
              grd_sig(i,j-1,k))*min(dx(i),xm(i)*dyac(j-1), &
              xm(i)*ym(j-1)*dz(k))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*xm(i)*dyac(j)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                0.75d0*pp*dx2(i)*sqrt(1d0-grd_yarr(j)**2)/(dy(j)*dx3(i)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dyac(j) + &
                (grd_sig(i,j-1,k)+grd_cap(ig,i,j-1,k))*dyac(j-1))
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                2.0d0*sqrt(1d0-grd_yarr(j)**2)*dx(i) / &
                (dy(j)*dx3(i)*help*thelp**2)
        endif

!
!-- calculating j->j+1 leakage opacity
        if(j==grd_ny) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j+1,k)+ &
              grd_sig(i,j+1,k))*min(dx(i),xm(i)*dyac(j+1), &
              xm(i)*ym(j+1)*dz(k))*thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*xm(i)*dyac(j)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                0.75d0*pp*dx2(i)*sqrt(1d0-grd_yarr(j+1)**2)/(dy(j)*dx3(i)*thelp)
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dyac(j) + &
                (grd_sig(i,j+1,k)+grd_cap(ig,i,j+1,k))*dyac(j+1))
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                2.0d0*sqrt(1d0-grd_yarr(j+1)**2)*dx(i) / &
                (dy(j)*dx3(i)*help*thelp**2)
        endif

!-- 0 azimuthal leakage opacities if nz=1
        if(grd_nz==1) cycle
!
!-- calculating k->k-1 leakage opacity
        if(k==1) then
           lhelp = (grd_cap(ig,i,j,grd_nz)+ &
              grd_sig(i,j,grd_nz))*min(dx(i),xm(i)*dyac(j), &
              xm(i)*ym(j)*dz(grd_nz))*thelp<prt_tauddmc
           khelp = grd_nz
        else
           lhelp = (grd_cap(ig,i,j,k-1)+ &
              grd_sig(i,j,k-1))*min(dx(i),xm(i)*dyac(j), &
              xm(i)*ym(j)*dz(k-1))*thelp<prt_tauddmc
           khelp = k-1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*xm(i)*ym(j) * &
              dz(k)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(5,i,j,k)=grd_opacleak(5,i,j,k)+(specval*speclump)*&
                0.75d0*pp*dx2(i)*dyac(j)/(dy(j)*dx3(i)*dz(k))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dz(k) + &
                (grd_sig(i,j,khelp)+grd_cap(ig,i,j,khelp))*dz(khelp))
           grd_opacleak(5,i,j,k)=grd_opacleak(5,i,j,k)+(specval*speclump)*&
                2.0d0*dyac(j)*dx(i) / &
                (ym(j)*dy(j)*dz(k)*dx3(i)*help*thelp**2)
        endif

!
!-- calculating k->k+1 leakage opacity
        if(k==grd_nz) then
           lhelp = (grd_cap(ig,i,j,1)+ &
              grd_sig(i,j,1))*min(dx(i),xm(i)*dyac(j), &
              xm(i)*ym(j)*dz(1))*thelp<prt_tauddmc
           khelp = 1
        else
           lhelp = (grd_cap(ig,i,j,k+1)+ &
              grd_sig(i,j,k+1))*min(dx(i),xm(i)*dyac(j), &
              xm(i)*ym(j)*dz(k+1))*thelp<prt_tauddmc
           khelp = k+1
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*xm(i)*ym(j) * &
              dz(k)*thelp
           pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(6,i,j,k)=grd_opacleak(6,i,j,k)+(specval*speclump)*&
                0.75d0*pp*dx2(i)*dyac(j)/(dy(j)*dx3(i)*dz(k))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dz(k) + &
                (grd_sig(i,j,khelp)+grd_cap(ig,i,j,khelp))*dz(khelp))
           grd_opacleak(6,i,j,k)=grd_opacleak(6,i,j,k)+(specval*speclump)*&
                2.0d0*dyac(j)*dx(i) / &
                (ym(j)*dy(j)*dz(k)*dx3(i)*help*thelp**2)
        endif

     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity1
