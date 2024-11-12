* © 2023. Triad National Security, LLC. All rights reserved.
* This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National
* Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of
* Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad
* National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration.
* The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
* irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute
* copies to the public, perform publicly and display publicly, and to permit others to do so.
      module tbxsmod
c     ---------------
      implicit none
************************************************************************
* Opacity tables by C. Fontes
************************************************************************
c-- number of tabulated elements, densities, temperatures
      integer,parameter :: tb_nrho=17
      integer,parameter :: tb_ntemp=27
      integer :: tb_nelem=0     ! number of elements
      integer :: tb_nelem_em=0  ! number of elements with emission opacity
      integer :: tb_nrho_em=0   ! number of densities with emission opacity
c-- density, temperature table points
      real*8 :: tb_rho(tb_nrho) !(nrho)
      real*8 :: tb_temp(tb_ntemp) !(ntemp)
c-- raw table input
      real*8,allocatable,private :: tb_raw(:,:,:,:,:) !(ncol,ngr,ntemp,nrho,nelem)
c-- raw table for emission (size may differ from tb_raw)
      real*8,allocatable,private :: tb_em_raw(:,:,:,:,:) !(ncol,ngr,ntemp,nrho,nelem)
c-- grey scattering opacity
      real*8,allocatable :: tb_sig(:,:,:)
c-- group-coarsened table
      real*4,allocatable :: tb_cap(:,:,:,:) !(ng,ntemp,nrho,nelem)
c-- group-coarsened emission opacity table
      real*4,allocatable :: tb_em_cap(:,:,:,:) !(ng,ntemp,nrho_em,nelem_em)
c-- number of energy points
      integer,parameter,private :: ngr=14900
c-- indicies of structure elements (from input.str)
      integer, allocatable :: tb_ielem(:) !(nelem)
c-- maps of emission element/density index to index of that of absoprtion
      integer, allocatable :: tb_ielem_em_map(:)
      integer, allocatable :: tb_irho_em_map(:)
c-- mask of tb_ielem that indicates emission opacity tables
      logical, allocatable :: tb_ielem_em_mask(:,:) !(nelem,nrho)
c
      save
      public
c
      contains
c
c
c
      subroutine tbxs_dealloc
c     -----------------------
      implicit none
      if(allocated(tb_ielem)) deallocate(tb_ielem)
      if(allocated(tb_sig)) deallocate(tb_sig)
      if(allocated(tb_cap)) deallocate(tb_cap)
c
      if(allocated(tb_ielem_em_map)) deallocate(tb_ielem_em_map)
      if(allocated(tb_irho_em_map)) deallocate(tb_irho_em_map)
      if(allocated(tb_ielem_em_mask)) deallocate(tb_ielem_em_mask)
      if(allocated(tb_em_cap)) deallocate(tb_em_cap)
c
      end subroutine tbxs_dealloc
c
c
c
      subroutine read_tbxs(lemiss)
c     --------------------
      use miscmod, only:lcase
      use elemdatamod, only: elem_neldata, elem_data
      use inputstrmod, only: str_nabund,str_iabund
      implicit none
      logical, intent(in) :: lemiss
************************************************************************
* Read Chris Fontes opacity data tables to total opacity array.
************************************************************************
      integer,parameter :: ncol=7
      character(3),parameter :: word1='op_'
      character(4),parameter :: word2='_1Em'
      character(9),parameter :: word3='gcc.table'
      character(3),parameter :: word4='em_'
      integer :: istat,ierr
      integer :: ielem,itemp,irho,iirho,l,ll,irho_em
      real*8 :: dmy
      character(6) :: sdmy
      character(2) :: fid,fnum
      character(20) :: fname
      character(23) :: fname_em
      logical :: lexist_em
      logical :: lrho_mask(tb_nrho)
      real*8 :: temp_em(tb_ntemp)
c
c-- sanity check
      if(str_nabund==0) stop 'read_tbxs: str_nabund=0'
c
c-- find number of tabulated elements from iabund
      do ielem=1,elem_neldata
        if(any(str_iabund==ielem)) tb_nelem=tb_nelem+1
      enddo
c-- need at least one tabulated element in input.str
      if(tb_nelem<=0) stop 'read_tbxs: tb_nelem<=0'
c-- store correct indices from elem_data
      allocate(tb_ielem(tb_nelem))
      l = 0
      do ielem=1,elem_neldata
        if(any(str_iabund==ielem)) then
          l=l+1
          tb_ielem(l)=ielem
        endif
      enddo
c
c-- raw table
      allocate(tb_raw(ncol,ngr,tb_ntemp,tb_nrho,tb_nelem))
c
c-- elements
      do l=1,tb_nelem
c-- skip if not in input structure
        ielem=tb_ielem(l)
c-- file element
        fid=lcase(trim(elem_data(ielem)%sym))
      do irho=tb_nrho,1,-1
c-- density value
        iirho=irho+3
        tb_rho(tb_nrho-irho+1)=10d0**(-iirho)
c-- file density id
        if(iirho/10==0) then
          write(fnum,'("0"i1)') iirho
        else
          write(fnum,'(i2)') iirho
        endif
c-- file name
        if(fid(2:2) == ' ') then
          fname=trim(word1//fid(1:1)//word2//fnum//word3)
        else
          fname=trim(word1//fid//word2//fnum//word3)
        endif
        open(4,file='Table/'//adjustl(fname),status='old',
     &     action='read',iostat=istat)
c-- require all possible data (for now)
        if(istat/=0) then
          write(6,*) 'file: ',fname
          stop 'read_tbxs: missing file'
        endif
        do itemp=1,tb_ntemp
c-- temperature value
          read(4,*) sdmy, dmy, tb_temp(itemp)
          read(4,*)
c-- all data at temp-rho(-elem) point
          read(4,*,iostat=ierr) tb_raw(:,:,itemp,tb_nrho-irho+1,l)
c
          if(ierr/=0) stop 'read_tbxs format err: body'
        enddo
c-- ensure no residual file data
        read(4,*,iostat=ierr) sdmy
        if(ierr/=-1) then
          write(6,*) 'sdmy: ',sdmy
          write(6,*) 'file: ',fname
          stop 'read_tbxs: body too long'
        endif
        close(4)
      enddo
      enddo
c
c-- check if there exists
      if (lemiss) then
c
c-- create emission subset mask
        allocate(tb_ielem_em_mask(tb_nelem, tb_nrho))
        tb_ielem_em_mask = .false.
c
c-- get number of elements with emission opacity tables
        do l=1,tb_nelem
        do irho=1,tb_nrho
c-- log10(density value)
          iirho=irho+3
c-- file element
          ielem = tb_ielem(l)
          fid=lcase(trim(elem_data(ielem)%sym))
c-- file density id
          if(iirho/10==0) then
            write(fnum,'("0"i1)') iirho
          else
            write(fnum,'(i2)') iirho
          endif
c-- construct file name and inquire existence
          if(fid(2:2) == ' ') then
            fname_em=trim(word4//word1//fid(1:1)//word2//fnum//word3)
          else
            fname_em=trim(word4//word1//fid//word2//fnum//word3)
          endif
          inquire(file='Table/'//adjustl(fname_em),exist=lexist_em)
c-- store file's existance and increment number of emission tables
          tb_ielem_em_mask(l,irho) = lexist_em
        enddo !irho
c-- check if any density had an emission opacity
          if(any(tb_ielem_em_mask(l,:))) then
            tb_nelem_em=tb_nelem_em+1
          endif
        enddo !l
c
c-- build element emission index to absorption index map
        allocate(tb_ielem_em_map(tb_nelem_em))
        ll = 1
        do l = 1, tb_nelem
          if(any(tb_ielem_em_mask(l,:))) then
            tb_ielem_em_map(ll) = l
            ll = ll + 1
          endif
        enddo
c
c-- element with max number of density values sets nrho_em
        lrho_mask = .false.
        do irho = 1, tb_nrho
          lrho_mask(irho) = any(tb_ielem_em_mask(:,irho))
        enddo
        tb_nrho_em = count(lrho_mask)
        write(*,*) 'tb_nelem_em = ', tb_nelem_em
        write(*,*) 'tb_nrho_em = ', tb_nrho_em
c-- build density emission index to full density index map
        allocate(tb_irho_em_map(tb_nrho_em))
        irho_em = 1
        do irho = 1, tb_nrho
          if(lrho_mask(irho)) then
            tb_irho_em_map(irho_em) = irho
            irho_em = irho_em + 1
          endif
        enddo
c
c-- raw emission table (TODO: do not assume same rho,T vals)
c-- TODO: allocate by number of available density points
        allocate(tb_em_raw(ncol,ngr,tb_ntemp,tb_nrho_em,tb_nelem_em))
        tb_em_raw = 0d0
        do ll=1,tb_nelem_em
        do irho_em=1,tb_nrho_em
          l = tb_ielem_em_map(ll)
          irho = tb_irho_em_map(irho_em)
c-- ignore elements without emission table
          if (.not.tb_ielem_em_mask(l,irho)) cycle
c-- log10(density value)
          iirho=irho+3
c-- file element
          ielem = tb_ielem(l)
          fid=lcase(trim(elem_data(ielem)%sym))
c-- file density id
          if(iirho/10==0) then
            write(fnum,'("0"i1)') iirho
          else
            write(fnum,'(i2)') iirho
          endif
c-- construct file name and inquire existence
          if(fid(2:2) == ' ') then
            fname_em=trim(word4//word1//fid(1:1)//word2//fnum//word3)
          else
            fname_em=trim(word4//word1//fid//word2//fnum//word3)
          endif
c-- read file
          open(4,file='Table/'//adjustl(fname_em),status='old',
     &         action='read',iostat=istat)
c-- require all possible data (for now)
          if(istat/=0) then
            write(6,*) 'file: ',fname_em
            stop 'read_tbxs: missing emission file'
          endif
          do itemp=1,tb_ntemp
c-- temperature value
            read(4,*) sdmy, dmy, temp_em(itemp)
            read(4,*)
            write(6,*) sdmy, dmy, temp_em(itemp)
c-- all data at temp-rho(-elem) point (TODO: fix rho index)
            read(4,*,iostat=ierr) tb_em_raw(:,:,itemp,irho_em,ll)
c
            if(ierr/=0) stop 'read_tbxs format err: body'
          enddo
c-- ensure no residual file data
          read(4,*,iostat=ierr) sdmy
          if(ierr/=-1) then
            write(6,*) 'sdmy: ',sdmy
            write(6,*) 'file: ',fname_em
            stop 'read_tbxs: body too long'
          endif
          close(4)
        enddo !irho_em
        enddo !l
c
c-- check emission table temperatures are the same
        if (any(abs(temp_em - tb_temp) > 1e-6 * tb_temp)) then
          write(*,*) 'tb_temp = ', tb_temp
          write(*,*) 'temp_em = ', temp_em
          stop 'read_tbxs: temp_em /= tb_temp'
        endif
c
      write(*,*) 'raw error = ',
     &     (sum(tb_em_raw(:,1,1,1),
     &     tb_em_raw(:,1,1,1) < huge(tb_em_raw(6,:,1,1,1)))
     &     - sum(tb_raw(:,1,1,1),
     &     tb_em_raw(:,1,1,1) < huge(tb_em_raw(6,:,1,1,1)))) /
     &     sum(tb_raw(:,1,1,1),
     &     tb_em_raw(:,1,1,1) < huge(tb_em_raw(6,:,1,1,1)))
      endif
c
      end subroutine read_tbxs
c
c
c
      subroutine coarsen_tbxs(lopac,lemiss,ng,wl)
c     ------------------------------------------------------------------
      use miscmod, only:binsrch
      use physconstmod
      implicit none
      logical, intent(in) :: lopac(4)
      logical, intent(in) :: lemiss
      integer, intent(in) :: ng
      real*8, intent(in) :: wl(ng+1)
************************************************************************
* Coursen Chris Fontes opacity tables.
************************************************************************
      integer :: ig,igr,itemp,irho,ielem,n
      integer :: ig1,ig2
      real*4 :: cap,cap1,cap2,capl,capr
      real*8 :: help1,help2,help3,help4,wll,wlr
      real*8 :: plck, f, uerg3, emu
      real*8,allocatable :: evinvarr(:)
c
c-- coarse energy array
      allocate(evinvarr(ng+1))
c-- final opacity tables
      allocate(tb_sig(tb_ntemp,tb_nrho,tb_nelem))
      allocate(tb_cap(ng,tb_ntemp,tb_nrho,tb_nelem))
c-- initialize
      evinvarr=pc_ev*wl/(pc_h*pc_c)
      tb_sig=0d0
      tb_cap=0.
c-- quick exit
      if(all(lopac)) return
c
      do ielem=1,tb_nelem
      do irho=1,tb_nrho
      do itemp=1,tb_ntemp
c-- scattering opacity
        if(.not.lopac(4)) tb_sig(itemp,irho,ielem)=
     &     sum(tb_raw(7,:,itemp,irho,ielem))/dble(ngr)
c-- absorption opacity
      do igr=ngr-1,1,-1
c-- binary search group containing igr point
        wll=1d0/tb_raw(1,igr+1,itemp,irho,ielem)
        wlr=1d0/tb_raw(1,igr,itemp,irho,ielem)
        ig1=binsrch(wll,evinvarr,ng+1,.true.)
        ig1=max(ig1,1)
        ig2=binsrch(wlr,evinvarr,ng+1,.true.)
        ig2=min(ig2,ng)
c-- exclude points out of wl domain
        if(ig2==0 .or. ig1==ng+1) cycle
        do ig=ig1,ig2
c-- interpolate trapezoid
          help1=max(wll,evinvarr(ig))
          help2=min(wlr,evinvarr(ig+1))
          help3=wlr*(help1-wll)/(help1*(wlr-wll)) ! high energy, low wavelength basis
          help4=wll*(wlr-help2)/(help2*(wlr-wll)) ! low energy, high wavelength basis
c-- add opacity contributions
          cap=0.
c-- bound-bound
          if(.not.lopac(1)) then
            capl=sngl(tb_raw(4,igr+1,itemp,irho,ielem))
            capr=sngl(tb_raw(4,igr,itemp,irho,ielem))
            cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
            cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
            cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
          endif
c-- bound-free
          if(.not.lopac(2)) then
            capl=sngl(tb_raw(5,igr+1,itemp,irho,ielem))
            capr=sngl(tb_raw(5,igr,itemp,irho,ielem))
            cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
            cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
            cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
          endif
c-- free-free
          if(.not.lopac(3)) then
            capl=sngl(tb_raw(6,igr+1,itemp,irho,ielem))
            capr=sngl(tb_raw(6,igr,itemp,irho,ielem))
            cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
            cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
            cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
          endif
c-- tally absorption opacity in group
          tb_cap(ig,itemp,irho,ielem) =
     &       tb_cap(ig,itemp,irho,ielem)+cap
        enddo
      enddo
c-- average opacity
      tb_cap(:,itemp,irho,ielem) =
     &   tb_cap(:,itemp,irho,ielem) *
     &   sngl(evinvarr(2:)*evinvarr(:ng) /
     &   (evinvarr(2:)-evinvarr(:ng)))
      enddo
      enddo
      enddo
c
c-- deallocate large raw table
      deallocate(tb_raw)
c
c-- emission opacity
      if (lemiss) then
c
c-- coarsen emission opacity groups
        allocate(tb_em_cap(ng,tb_ntemp,tb_nrho_em,tb_nelem_em))
c
        do ielem=1,tb_nelem_em
        do irho=1,tb_nrho_em
        do itemp=1,tb_ntemp
        do igr=ngr-1,1,-1
c-- binary search group containing igr point
          wll=1d0/tb_em_raw(1,igr+1,itemp,irho,ielem)
          wlr=1d0/tb_em_raw(1,igr,itemp,irho,ielem)
          ig1=binsrch(wll,evinvarr,ng+1,.true.)
          ig1=max(ig1,1)
          ig2=binsrch(wlr,evinvarr,ng+1,.true.)
          ig2=min(ig2,ng)
c-- exclude points out of wl domain
          if(ig2==0 .or. ig1==ng+1) cycle
          do ig=ig1,ig2
c-- interpolate trapezoid
            help1=max(wll,evinvarr(ig))
            help2=min(wlr,evinvarr(ig+1))
            help3=wlr*(help1-wll)/(help1*(wlr-wll)) ! high energy, low wavelength basis
            help4=wll*(wlr-help2)/(help2*(wlr-wll)) ! low energy, high wavelength basis
c-- add emission opacity contributions
            cap=0.
c$$$c-- bound-bound
c$$$            if(.not.lopac(1)) then
c$$$              capl=sngl(tb_em_raw(3,igr+1,itemp,irho,ielem))
c$$$              capr=sngl(tb_em_raw(3,igr,itemp,irho,ielem))
c$$$              cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
c$$$              cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
c$$$              cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
c$$$            endif
c$$$c-- bound-free
c$$$          if(.not.lopac(2)) then
c$$$            capl=sngl(tb_em_raw(4,igr+1,itemp,irho,ielem))
c$$$            capr=sngl(tb_em_raw(4,igr,itemp,irho,ielem))
c$$$            cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
c$$$            cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
c$$$            cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
c$$$          endif
c$$$c-- free-free
c$$$          if(.not.lopac(3)) then
c$$$            capl=sngl(tb_em_raw(5,igr+1,itemp,irho,ielem))
c$$$            capr=sngl(tb_em_raw(5,igr,itemp,irho,ielem))
c$$$            cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
c$$$            cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
c$$$            cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
c$$$          endif
          capl=min(sngl(tb_em_raw(6,igr+1,itemp,irho,ielem)), 1e34)
          capr=min(sngl(tb_em_raw(6,igr,itemp,irho,ielem)), 1e34)
          cap1=capl*sngl(1d0-help3)+capr*sngl(help3)
          cap2=capl*sngl(help4)+capr*sngl(1d0-help4)
          cap=cap+.5*sngl((help2-help1)/(help1*help2))*(cap1+cap2)
c$$$c-- convert to emission opacity
c$$$          f = 1d0/evinvarr(ig) ![eV]
c$$$          uerg3 = (pc_ev*f) * (pc_ev*f) * (pc_ev*f) ![erg^3]
c$$$          emu = exp(-f/tb_temp(itemp)) ![]
c$$$          plck = 9.611817d58*uerg3*emu/(1d0-emu)
c$$$          if (plck /= plck) then
c$$$            write(*,*)
c$$$            write(*,*) 'plck is NaN ...'
c$$$            write(*,*) 'f, tb_temp(itemp), uerg3, emu = ',
c$$$     &           f, tb_temp(itemp), uerg3, emu
c$$$          endif
c$$$          cap = sngl(dble(cap) /
c$$$     &         (tb_rho(tb_nrho-tb_irho_em_map(irho)+1)*plck)) ![cm^2/g]
c$$$          if (cap /= cap) then
c$$$            cap = 0.0
c$$$          elseif (cap > huge(cap)) then
c$$$            cap = 0.0
c$$$          endif
c-- tally emission opacity in group
          tb_em_cap(ig,itemp,irho,ielem) =
     &       tb_em_cap(ig,itemp,irho,ielem)+cap
        enddo
      enddo
c-- TODO: remove this Planck factor if contributions above become emission opacity
      ! erg3(i) =  hnu(i)^3       !(in ergs^3)
      ! u(i) = hnu(i) / kT
      ! emu(i) = exp(-u(i))

c-- average emission opacity
      tb_em_cap(:,itemp,irho,ielem) =
     &   tb_em_cap(:,itemp,irho,ielem) *
     &   sngl(evinvarr(2:)*evinvarr(:ng) /
     &   (evinvarr(2:)-evinvarr(:ng)))
      enddo
      enddo
      enddo
c
c-- deallocate large raw emission table
        deallocate(tb_em_raw)
      endif
c-- deallocate coarse energy array
      deallocate(evinvarr)
c
c-- convert temperature from eV to K
      tb_temp=tb_temp*pc_ev/pc_kb
c
c-- print alloc size (keep this updated)
c---------------------------------------
      n=tb_nelem*tb_nrho*tb_ntemp*(8+ng*4)/1024 !kB
      if (lemiss) n = n + tb_nelem_em*1*tb_ntemp*(ng*4)/1024 !kB
      write(6,*) 'ALLOC tbxs     :',n,"kB",n/1024,"MB",n/1024**2,"GB"
c
      write(*,*)
      write(*,*) 'cap = ', tb_cap(:,1,1,1)
      write(*,*)
      write(*,*) 'em_cap = ', tb_em_cap(:,1,1,1)
      write(*,*)
      write(*,*) 'error = ',
     &     (sum(tb_em_cap(:,1,1,1),
     &     tb_em_cap(:,1,1,1) < huge(tb_em_cap(:,1,1,1)))
     &     - sum(tb_cap(:,1,1,1),
     &     tb_em_cap(:,1,1,1) < huge(tb_em_cap(:,1,1,1)))) /
     &     sum(tb_cap(:,1,1,1),
     &     tb_em_cap(:,1,1,1) < huge(tb_em_cap(:,1,1,1)))
c
      stop 'testing stop ...'
c
      end subroutine coarsen_tbxs
c
      end module tbxsmod
