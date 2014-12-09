subroutine leakage_opacity3

  use miscmod
  use gridmod
  use groupmod
  use timestepmod
  use particlemod
  use physconstmod
  implicit none

!##################################################
  !This subroutine computes
  !DDMC 3D lumped leakage opacities.
!##################################################
  logical :: lhelp
  integer :: i,j,k, ig
  real*8 :: thelp, help
  real*8 :: speclump, specval
  real*8 :: specarr(grp_ng)
  real*8 :: pp, alb, eps, beta
!-- statement functions
  integer :: l
  real*8 :: dx,dy,dz
  dx(l) = grd_xarr(l+1) - grd_xarr(l)
  dy(l) = grd_yarr(l+1) - grd_yarr(l)
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
     specarr = specintv(1d0/grd_temp(i,j,k))
     help = min(dx(i),dy(j),dz(k))*thelp
     speclump = 1d0/sum(specarr,grd_cap(:,i,j,k)*help>=prt_taulump)
!-- lumping opacity
     do ig=1,grp_ng
        if(grd_cap(ig,i,j,k)*min(dx(i),dy(j),dz(k))*thelp < prt_taulump) cycle
!
!-- obtaining spectral weight
        specval = specarr(ig)
!
!-- calculating i->i-1 leakage opacity
        if(i==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i-1,j,k)+ &
              grd_sig(i-1,j,k))*min(dx(i-1),dy(j),dz(k)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dx(i))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i-1,j,k)+grd_cap(ig,i-1,j,k))*dx(i-1))*thelp
           grd_opacleak(1,i,j,k)=grd_opacleak(1,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dx(i)*thelp)
        endif

!
!-- calculating i->i+1 leakage opacity
        if(i==grd_nx) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i+1,j,k)+ &
              grd_sig(i+1,j,k))*min(dx(i+1),dy(j),dz(k)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dx(i)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dx(i))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dx(i)+&
                (grd_sig(i+1,j,k)+grd_cap(ig,i+1,j,k))*dx(i+1))*thelp
           grd_opacleak(2,i,j,k)=grd_opacleak(2,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dx(i)*thelp)
        endif

!
!-- calculating j->j-1 leakage opacity
        if(j==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j-1,k)+ &
              grd_sig(i,j-1,k))*min(dx(i),dy(j-1),dz(k)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dy(j)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dy(j)+&
                (grd_sig(i,j-1,k)+grd_cap(ig,i,j-1,k))*dy(j-1))*thelp
           grd_opacleak(3,i,j,k)=grd_opacleak(3,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!
!-- calculating j->j+1 leakage opacity
        if(j==grd_ny) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j+1,k)+ &
              grd_sig(i,j+1,k))*min(dx(i),dy(j+1),dz(k)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dy(j)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dy(j))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dy(j)+&
                (grd_sig(i,j+1,k)+grd_cap(ig,i,j+1,k))*dy(j+1))*thelp
           grd_opacleak(4,i,j,k)=grd_opacleak(4,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dy(j)*thelp)
        endif

!
!-- calculating k->k-1 leakage opacity
        if(k==1) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j,k-1)+ &
              grd_sig(i,j,k-1))*min(dx(i),dy(j),dz(k-1)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dz(k)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(5,i,j,k)=grd_opacleak(5,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dz(k))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dz(k)+&
                (grd_sig(i,j,k-1)+grd_cap(ig,i,j,k-1))*dz(k-1))*thelp
           grd_opacleak(5,i,j,k)=grd_opacleak(5,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dz(k)*thelp)
        endif

!
!-- calculating k->k+1 leakage opacity
        if(k==grd_nz) then
           lhelp = .true.
        else
           lhelp = (grd_cap(ig,i,j,k+1)+ &
              grd_sig(i,j,k+1))*min(dx(i),dy(j),dz(k+1)) * &
              thelp<prt_tauddmc
        endif
!
        if(lhelp) then
!-- DDMC interface
           help = (grd_cap(ig,i,j,k)+grd_sig(i,j,k))*dz(k)*thelp
           alb = grd_fcoef(i,j,k)*grd_cap(ig,i,j,k)/ &
                (grd_cap(ig,i,j,k)+grd_sig(i,j,k))
           eps = (4d0/3d0)*sqrt(3d0*alb)/(1d0+pc_dext*sqrt(3d0*alb))
           beta = 1.5d0*alb*help**2+sqrt(3d0*alb*help**2 + &
                2.25d0*alb**2*help**4)
           pp = 0.5d0*eps*beta/(beta-0.75*eps*help)
!              pp = 4d0/(3d0*help+6d0*pc_dext)
           grd_opacleak(6,i,j,k)=grd_opacleak(6,i,j,k)+(specval*speclump)*&
                0.5d0*pp/(thelp*dz(k))
        else
!-- DDMC interior
           help = ((grd_sig(i,j,k)+grd_cap(ig,i,j,k))*dz(k)+&
                (grd_sig(i,j,k+1)+grd_cap(ig,i,j,k+1))*dz(k+1))*thelp
           grd_opacleak(6,i,j,k)=grd_opacleak(6,i,j,k)+(specval*speclump)*&
                (2d0/3d0)/(help*dz(k)*thelp)
        endif

     enddo !ig
  enddo !i
  enddo !j
  enddo !k


end subroutine leakage_opacity3
