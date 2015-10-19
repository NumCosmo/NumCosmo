! code for evaluation of CMB likelihood in pixel space
! Based on BFLIKE by L. Colombo <colombo@usc.edu>
! Algorithm described in Tegmark & de Oliveira-Costa  arXiv:astro-ph/0012120
!*************************************************************************

module bflike_QU
  use healpix_types
  use fitstools_smw
  
  implicit none
  
  public init_pix_like, get_pix_loglike, clean_pix_like
  private
  
  integer,parameter,public :: myTT=1,myEE=2,myBB=3,myTE=4,myTB=5,myEB=6
  integer,save ::ntemp,npol,ntp,ntot,lmax,lsw,npixtot
  
  real(dp),allocatable,dimension(:,:),save ::NCVM,clnorm
  real(dp),allocatable,dimension(:,:) :: S ,Spol
  real(dp),allocatable,dimension(:),save:: dt,czs
  real(sp),allocatable,dimension(:),save :: c1s,c2s,s1s,s2s 
  real(sp),allocatable,dimension(:,:),save ::pls
  real(dp),allocatable,dimension(:):: plm,f1,f2,auxdt
  real(dp),allocatable,dimension(:,:):: cls
  
contains
  
  subroutine init_pix_like(clik_bflike_dir)
    implicit none
    character(LEN=*),intent(in)::clik_bflike_dir
    integer(i8b) :: npix
    integer(i4b) :: ordering,nside,nmaps,nsidemask
    integer(i4b)::i,l,uw,ub,j,unit,ct,lswitch
    real(dp) :: theta,phi,fct,fct2,chngconv,ell,nullval,calibration,dummy
    real(dp),allocatable,dimension(:,:) :: map,NCVMfull,beam
    real(dp),allocatable,dimension(:,:) :: mask,clfid
    real(dp),allocatable,dimension(:) :: alltheta,allphi,pl
    real(dp),allocatable,dimension(:,:) :: Tvec,Pvec
    integer(i4b),allocatable,dimension(:) :: pixact
    CHARACTER(len=80),DIMENSION(74) :: header1,header2,headercl
    character(len=20),dimension(3) :: units
    logical ::anynull
    real(dp) :: modv
    real(dp) :: tt,qq,uu,tq,tu,qu
    real(dp) :: a1,a2,c1c2,s1s2,c1s2,s1c2
    integer :: iele,nele
    
    character(len=256) :: mapTQU,maskTQU,ncvmfile,anglesfile,clfiducial,beamfile
    
    namelist /inputs/ mapTQU,maskTQU,ncvmfile,anglesfile,clfiducial,beamfile,lmax,lswitch,calibration
    
    write(*,'(a)') 'BFLIKE: Reading parameter file'
    !setting defaults
    mapTQU='mapTQU.fits'    
    maskTQU='maskTQU.fits'    
    ncvmfile='covmat.dat'
    anglesfile='pix2ang_nest_ns16.txt'
    clfiducial='fiducial.dat'
    beamfile='beamfile.fits'
    lmax=64
    lswitch=33
    calibration=1e6
    
    !read   
    open(11,file=trim(clik_bflike_dir)//'/params_bflike.ini',status='old')
    read(11,inputs)
    close(11)
    lsw=lswitch
    write(*,'(a)') 'BFLIKE: mapTQU: '//trim(mapTQU)
    write(*,'(a)') 'BFLIKE: maskTQU: '//trim(maskTQU)
    write(*,'(a)') 'BFLIKE: ncvmfile: '//trim(ncvmfile)
    write(*,'(a)') 'BFLIKE: anglesfile: '//trim(anglesfile)
    write(*,'(a)') 'BFLIKE: clfiducial: '//trim(clfiducial)
    write(*,'(a)') 'BFLIKE: beamfile: '//trim(beamfile)
    write(*,'(a,i6)') 'BFLIKE: lmax: ',lmax
    write(*,'(a,i6)') 'BFLIKE: lswitch: ',lsw
    write(*,'(a,e13.3)') 'BFLIKE: calibration: ',calibration
    
    write(*,'(a)') 'BFLIKE: Done'
    
    !data and NCVM
    npix=getsize_fits(trim(clik_bflike_dir)//'/'//trim(maskTQU),nside=nsidemask,ordering=ordering,nmaps=nmaps)
    if (nmaps .lt. 2) stop 'BFLIKE: Wrong number of fields in maskTQU'
    allocate(mask(0:12*nsidemask*nsidemask-1,1:nmaps))
    CALL read_dbintab(trim(clik_bflike_dir)//'/'//trim(maskTQU),mask,12*nsidemask*nsidemask,nmaps,nullval,anynull,units)
    if (ordering .eq. 1)write(*,'(a)') 'BFLIKE: WARNING!! map in ring'
    write(0,'(a)') 'BFLIKE: Reading Mask' 
    ntemp=0
    npol=0
    do i=0,12*nsidemask*nsidemask-1
       if (abs(mask(i,1)) .gt. 0.5) ntemp=ntemp+1
       if (abs(mask(i,2)) .gt. 0.5) npol=npol+1
    end do
    
    ntp = ntemp + npol
    ntot = ntemp +npol*2
    
    allocate(dt(ntot))
    allocate(pixact(ntot))
    
    npix=getsize_fits(trim(clik_bflike_dir)//'/'//trim(mapTQU),nside=nside,ordering=ordering,nmaps=nmaps)
    if (nside .ne. nsidemask) stop 'BFLIKE: Map and mask have different nside'
    if (nmaps .ne. 3) stop 'BFLIKE: Wrong number of fields in mapTQU'
    allocate(map(0:12*nside*nside-1,1:3))
    CALL read_dbintab(trim(clik_bflike_dir)//'/'//trim(mapTQU),map,12*nside*nside,3,nullval,anynull,units)
    if (ordering .eq. 1) write(*,'(a)') 'BFLIKE: WARNING!! map in ring'
    write(*,'(a)') 'BFLIKE: Reading Map'
    
    map=map*calibration
    npixtot=12*nside*nside
    
    ct=0
    do j=1,3
       do i=0,npixtot-1
          if (abs(mask(i,j)) .gt. 0.5) then
             ct=ct+1
             dt(ct)=map(i,j)
             pixact(ct)=(j-1)*npixtot+i+1 
          end if
       end do
    end do
    deallocate(map,mask)
    
    allocate(S(ntot,ntot))
    allocate(Spol(npol*2,npol*2))
    allocate(NCVMfull(3*npixtot,3*npixtot))
    allocate(NCVM(ntot,ntot))
    
    write(*,'(a)') 'BFLIKE: Reading covmat'
    ! this fails on old intel fortran !!!!
    !open(newunit=unit,file=trim(clik_bflike_dir)//'/'//trim(ncvmfile),form='binary',status='old')
    unit = 30
    open(unit=unit,file=trim(clik_bflike_dir)//'/'//trim(ncvmfile),access='stream',status='old')
    read(unit) NCVMfull
    close(unit)
    
    do i = 1,ntot
       do j=1,ntot
          NCVM(i,j)=NCVMfull(pixact(i),pixact(j))*calibration*calibration 
       end do
    end do
    
    deallocate(NCVMfull)
    
    ! geometry
    allocate(Tvec(ntemp,3))
    allocate(Pvec(npol,3))
    
    !open(newunit=unit,file=trim(clik_bflike_dir)//'/'//trim(anglesfile),form='formatted',status='old')
    open(unit=unit,file=trim(clik_bflike_dir)//'/'//trim(anglesfile),form='formatted',status='old')
    allocate(alltheta(0:npixtot-1),allphi(0:npixtot-1))
    do i=0,npixtot-1
       read(unit,*) dummy,alltheta(i),allphi(i)
    enddo
    close(unit)
    
    do i=1,ntemp
       theta=alltheta(pixact(i)-1)
       phi=allphi(pixact(i)-1) 
       Tvec(i,1) = sin(theta)*cos(phi)
       Tvec(i,2) = sin(theta)*sin(phi)
       Tvec(i,3) = cos(theta)
       modv = sqrt(sum(Tvec(i,:)**2))
       Tvec(i,:) = Tvec(i,:)/modv
    end do
    do i=1,npol
       theta=alltheta(pixact(i+ntemp)-npixtot-1)
       phi=allphi(pixact(i+ntemp)-npixtot-1)   
       Pvec(i,1) = sin(theta)*cos(phi)
       Pvec(i,2) = sin(theta)*sin(phi)
       Pvec(i,3) = cos(theta)
       modv = sqrt(sum(Pvec(i,:)**2))
       Pvec(i,:) = Pvec(i,:)/modv
    end do
    
    deallocate(pixact,alltheta,allphi)
    allocate(auxdt(ntot))
    
    !beam and window functions and various normalizations
    
    allocate(clnorm(2:lmax,6))
    allocate(beam(0:4*nside,1:3))
    allocate(clfid(2:lmax,4))
    
    call fits2cl(trim(clik_bflike_dir)//'/'//trim(beamfile),beam,lmax,3,headercl)
    
    fct2 = 1.d0/24.d0
    fct = sqrt(fct2)
    
    clnorm(2,myTT) = beam(2,1)*beam(2,1)/2.4d0
    clnorm(2,myTE) = beam(2,1)*beam(2,2)*fct/2.4d0
    clnorm(2,myTB) = beam(2,1)*beam(2,3)*fct/2.4d0
    clnorm(2,myEE) = beam(2,2)*beam(2,2)*fct2/2.4d0
    clnorm(2,myBB) = beam(2,3)*beam(2,3)*fct2/2.4d0
    clnorm(2,myEB) = beam(2,2)*beam(2,3)*fct2/2.4d0
    
    do l=3,lmax
       ell = real(l,kind=dp)
       
       fct2 = 1.d0/((ell+2.d0)*(ell+1.d0)*ell*(ell-1.d0))
       fct = sqrt(fct2)
       chngconv = (2.d0*ell +1.d0)/2.d0/(ell*(ell+1.d0))
       
       clnorm(l,myTT) = beam(l,1)*beam(l,1)*chngconv
       clnorm(l,myTE) = beam(l,1)*beam(l,2)*fct*chngconv
       clnorm(l,myTB) = beam(l,1)*beam(l,3)*fct*chngconv
       clnorm(l,myEE) = beam(l,2)*beam(l,2)*fct2*chngconv
       clnorm(l,myBB) = beam(l,3)*beam(l,3)*fct2*chngconv
       clnorm(l,myEB) = beam(l,2)*beam(l,3)*fct2*chngconv
    end do
    
    deallocate(beam)
    
    write(*,'(a)')'BFLIKE: Reading Fiducial'
    
    !open(newunit=unit,file=trim(clik_bflike_dir)//'/'//trim(clfiducial),status='old')
    open(unit=unit,file=trim(clik_bflike_dir)//'/'//trim(clfiducial),status='old')
    do l=2,lmax
       read(unit,*)ell,clfid(l,:)
    enddo
    close(unit)
    
    clfid=clfid*clnorm(:,1:4)
    
    allocate(pl(1:lmax))
    allocate(plm(1:lmax))
    allocate(f1(2:lmax))
    allocate(f2(2:lmax))
    
    S=0
    !TT diag
    tt = sum(clfid(lsw:lmax,myTT))
    tt = tt + 100.d0*clfid(2,myTT)*2.d0
    do i=1,ntemp
       S(i,i) = tt
    end do
    
    nele=int(ntp*(ntp-1)/2)
    allocate(czs(nele),c1s(nele),s1s(nele),c2s(nele),s2s(nele))
    
    iele=0
    do i=1,ntemp
       ! TT offdiag
       do j=i+1,ntemp
          iele=iele+1
          czs(iele) = sum(Tvec(j,:)*Tvec(i,:))
          pl(1) = czs(iele)
          pl(2) = 1.5d0*czs(iele)*czs(iele) -.5d0
          do l = 3,lmax
             pl(l) =(czs(iele)*(2*l -1)*pl(l-1) -(l-1)*pl(l-2))/l
          enddo
          
          tt = sum(clfid(lsw:lmax,myTT)*pl(lsw:lmax))
          tt = tt + 100.d0*clfid(2,myTT)*(1.d0 + pl(1))
          S(j,i) = tt
       enddo
       ! QT and UT
       do j=ntemp+1,ntp
          iele=iele+1
          czs(iele) = sum(Pvec(j-ntemp,:)*Tvec(i,:))
          plm(1) = 0.d0
          plm(2) = 3.d0*(1.d0 -czs(iele)*czs(iele))
          do l = 3,lmax
             plm(l) =(czs(iele)*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
          enddo
          tq = -sum(clfid(lsw:lmax,myTE)*plm(lsw:lmax))
          tu=0.d0
          !          tu = -sum(clfid(lsw:lmax,myTB)*plm(lsw:lmax))
          
          call get_rotation_angle(PVec(j-ntemp,:),TVec(i,:),a1,a2)
          
          c1s(iele)=cos(a1)
          c2s(iele)=cos(a2)
          s1s(iele)=-sin(a1)
          s2s(iele)=-sin(a2)
          
          S(j,i) = tq*c1s(iele) +tu*s1s(iele)
          S(j+npol,i)= -tq*s1s(iele) +tu*c1s(iele)
       enddo
    end do
    
    !QQ and UU diag
    plm(1) = 0.d0
    plm(2) = 3.d0
    f1(2)  = 12.d0
    f2(2)  = -12.d0
    do l = 3,lmax
       plm(l) =((2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l-2)
       f1(l) =2.d0*(-(l-4.d0)*plm(l)+(l+2)*plm(l-1))
       f2(l) = 4.d0*(-(l-1)*plm(l) +(l+2)*plm(l-1))
    enddo
    
    qq = sum(clfid(lsw:lmax,myEE)*f1(lsw:lmax) - clfid(lsw:lmax,myBB)*f2(lsw:lmax))
    uu = sum(clfid(lsw:lmax,myBB)*f1(lsw:lmax) - clfid(lsw:lmax,myEE)*f2(lsw:lmax))
    qu=0.d0
    !    qu = sum((f1(lsw:lmax) +f2(lsw:lmax))*clfid(lsw:lmax,myEB))
    
    do i=ntemp+1,ntp
       S(i,i)      =qq
       S(i+npol,i) =qu
       S(i+npol,i+npol) = uu
    end do
 
    !QQ and UU offdiag and QU
    do i=ntemp+1,ntp
       do j=i+1,ntp
          iele=iele+1
          czs(iele) = sum(Pvec(j-ntemp,:)*Pvec(i-ntemp,:))
          plm(1) = 0.d0
          plm(2) = 3.d0
          f1(2)  = 6.d0*(1d0+czs(iele)*czs(iele))
          f2(2)  = -12.d0*czs(iele)
          do l = 3,lmax
             plm(l) =(czs(iele)*(2*l -1)*plm(l-1) -(l+1)*plm(l-2))/(l -2)
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -czs(iele)*czs(iele)))*plm(l)+ &
                  (2*l+4)*czs(iele)*plm(l-1)
             f2(l) = 4.d0*(-(l-1)*czs(iele)*plm(l) +(l+2)*plm(l-1))
          enddo
          qq = sum(clfid(lsw:lmax,myEE)*f1(lsw:lmax) -clfid(lsw:lmax,myBB)*f2(lsw:lmax))
          uu = sum(clfid(lsw:lmax,myBB)*f1(lsw:lmax) -clfid(lsw:lmax,myEE)*f2(lsw:lmax))
          qu=0.d0
          !          qu = sum((f1(lsw:lmax) +f2(lsw:lmax))*clfid(lsw:lmax,myEB))
          
          call get_rotation_angle(PVec(j-ntemp,:),PVec(i-ntemp,:),a1,a2)

          c1s(iele)=cos(a1)
          c2s(iele)=cos(a2)
          s1s(iele)=-sin(a1)
          s2s(iele)=-sin(a2)
          
          c1c2 = c1s(iele)*c2s(iele)
          s1s2 = s1s(iele)*s2s(iele)
          c1s2 = c1s(iele)*s2s(iele)
          s1c2 = s1s(iele)*c2s(iele)
          S(j,i) = qq*c1c2 +uu*s1s2  +qu*(c1s2+s1c2)
          S(j+npol,i+npol) = qq*s1s2 +uu*c1c2 -qu*(c1s2+s1c2)
          S(j+npol,i) = -qq*s1c2 +uu*c1s2 +qu*(c1c2 -s1s2)
          S(i+npol,j) = -qq*c1s2 +uu*s1c2 +qu*(c1c2 -s1s2)
       end do
    end do
    
    deallocate(Pvec,Tvec,clfid,pl)
    
    NCVM = S+NCVM
    
    allocate(pls(1:lsw-1,nele))
    iele=0
    do i=1,ntemp
       ! TT offdiag
       do j=i+1,ntemp
          iele=iele+1
          pls(1,iele) = czs(iele)
          pls(2,iele) = 1.5d0*czs(iele)*czs(iele) -.5d0
          do l = 3,lsw-1
             pls(l,iele) =(czs(iele)*(2*l -1)*pls(l-1,iele) -(l-1)*pls(l-2,iele))/l
          enddo
       enddo
       ! QT and UT
       do j=ntemp+1,ntp
          iele=iele+1
          pls(1,iele) = 0.d0
          pls(2,iele) = 3.d0*(1.d0 -czs(iele)*czs(iele))
          do l = 3,lsw-1
             pls(l,iele) =(czs(iele)*(2*l -1)*pls(l-1,iele) -(l+1)*pls(l-2,iele))/(l-2)
          enddo
       enddo
    end do
    
    !QQ and UU offdiag and QU
    do i=ntemp+1,ntp
       do j=i+1,ntp
          iele=iele+1
          pls(1,iele) = 0.d0
          pls(2,iele) = 3.d0
          do l = 3,lsw-1
             pls(l,iele) =(czs(iele)*(2*l -1)*pls(l-1,iele) -(l+1)*pls(l-2,iele))/(l -2)
          enddo
       end do
    end do
    
    allocate(cls(2:lsw-1,6))
    
    
  end subroutine init_pix_like
  
  
  subroutine clean_pix_like
    implicit none

    if(allocated(S)) deallocate(S)
    if(allocated(NCVM)) deallocate(NCVM)
    if(allocated(clnorm)) deallocate(clnorm)
    if(allocated(c1s)) deallocate(c1s)
    if(allocated(c2s)) deallocate(c2s)
    if(allocated(s1s)) deallocate(s1s)
    if(allocated(s2s)) deallocate(s2s)
    
    if(allocated(pls)) deallocate(pls)
    if(allocated(plm)) deallocate(plm)
    if(allocated(f1)) deallocate(f1)
    if(allocated(f2)) deallocate(f2)
    if(allocated(cls)) deallocate(cls)
    if(allocated(czs)) deallocate(czs)
    
  end subroutine clean_pix_like
  
  
  subroutine get_pix_loglike(clsin,alike)
    !input clsin are l(l+1)C_l/2pi
    
    implicit none
    
    real(dp),intent(in)::clsin(2:,:)
    real(dp) :: logdet,argexp
    real(dp),intent(out) :: alike
    
    integer :: i,j,l,info,iele
    real(dp) :: tt,qq,uu,tq,tu,qu
    real(dp) :: c1c2,s1s2,c1s2,s1c2
    
    cls=0.d0
    cls = clsin(2:lsw-1,:)*clnorm(2:lsw-1,:)
    
    !  signal matrix
    S=0.d0
    
    !TT diag
    tt = sum(cls(2:lsw-1,myTT))
    do i=1,ntemp
       S(i,i) = tt
    end do
    
    iele=0
    do i=1,ntemp
       ! TT offdiag
       do j=i+1,ntemp
          iele=iele+1 
          S(j,i) = sum(cls(2:lsw-1,myTT)*pls(2:lsw-1,iele))
       enddo
       ! QT and UT
       do j=ntemp+1,ntp
          iele=iele+1
          tq = -sum(cls(2:lsw-1,myTE)*pls(2:lsw-1,iele))
          tu = -sum(cls(2:lsw-1,myTB)*pls(2:lsw-1,iele))
          S(j,i) = tq*c1s(iele) +tu*s1s(iele)
          S(j+npol,i)= -tq*s1s(iele) +tu*c1s(iele)
       enddo
    end do
    
    !QQ and UU diag
    plm(1) = 0.d0
    plm(2) = 3.d0
    f1(2)  = 12.d0
    f2(2)  = -12.d0
    do l = 3,lsw-1
       plm(l) =((2*l -1)*plm(l-1) -(l+1)* plm(l-2))/(l-2)
       f1(l) =2.d0*(-(l-4.d0)*plm(l)+(l+2)*plm(l-1))
       f2(l) = 4.d0*(-(l-1)*plm(l) +(l+2)*plm(l-1))
    enddo
    
    qq = sum(cls(2:lsw-1,myEE)*f1(2:lsw-1) - cls(2:lsw-1,myBB)*f2(2:lsw-1))
    uu = sum(cls(2:lsw-1,myBB)*f1(2:lsw-1) - cls(2:lsw-1,myEE)*f2(2:lsw-1))
    qu = sum((f1(2:lsw-1) +f2(2:lsw-1))*cls(2:lsw-1,myEB))
    
    do i=ntemp+1,ntp
       S(i,i)      =qq
       S(i+npol,i) =qu
       S(i+npol,i+npol) = uu
    end do
    
    !QQ and UU offdiag and QU
    
    do i=ntemp+1,ntp
       do j=i+1,ntp
          iele=iele+1
          f1(2)  = 6.d0*(1d0+czs(iele)*czs(iele))
          f2(2)  = -12.d0*czs(iele)
          do l = 3,lsw-1
             f1(l) =-(2*l-8 +l*(l-1)*(1.d0 -czs(iele)*czs(iele)))*pls(l,iele)+ &
                  (2*l+4)*czs(iele)*pls(l-1,iele)
             f2(l) = 4.d0*(-(l-1)*czs(iele)*pls(l,iele) +(l+2)*pls(l-1,iele))
          enddo
          qq = sum(cls(2:lsw-1,myEE)*f1(2:lsw-1) -cls(2:lsw-1,myBB)*f2(2:lsw-1))
          uu = sum(cls(2:lsw-1,myBB)*f1(2:lsw-1) -cls(2:lsw-1,myEE)*f2(2:lsw-1))
          qu = sum((f1(2:lsw-1) +f2(2:lsw-1))*cls(2:lsw-1,myEB))
          
          c1c2 = c1s(iele)*c2s(iele)
          s1s2 = s1s(iele)*s2s(iele)
          c1s2 = c1s(iele)*s2s(iele)
          s1c2 = s1s(iele)*c2s(iele)
          S(j,i) = qq*c1c2 +uu*s1s2  +qu*(c1s2+s1c2)
          S(j+npol,i+npol) = qq*s1s2 +uu*c1c2 -qu*(c1s2+s1c2)
          S(j+npol,i) = -qq*s1c2 +uu*c1s2 +qu*(c1c2 -s1s2)
          S(i+npol,j) = -qq*c1s2 +uu*s1c2 +qu*(c1c2 -s1s2)
       end do
    end do
    
    S = S+NCVM 
    
    !    cholesky decomposition:copy dt in auxdt
    auxdt =dt


    !    cholesky:  C^-1*X
!    call dposv('L',ntot,1,S,ntot,auxdt,ntot,info)
    !    cholesky: X*(C^-1*X)
    !write(0,*) 'Info: ',info
    
!    if (info.eq.0) then
!       argexp = sum(dt*auxdt)
       
!       logdet = 0.d0
!       do j=1,ntot
!          logdet =logdet +log(S(j,j))
!       enddo
       
!       logdet = 2.d0*logdet
!    else
!       argexp = 1.d37
!       logdet = 1.d37
!    endif
    
!    alike = -.5d0*(argexp+logdet)
    Spol = S(ntemp+1: ntot,ntemp+1: ntot)
    call dposv('L',npol*2,1,Spol,npol*2,auxdt(ntemp+1: ntot),npol*2,info)
    !    cholesky: X*(C^-1*X)
    write(0,*) 'Info: ',info

    if (info.eq.0) then
       argexp = sum(dt(ntemp+1: ntot)*auxdt(ntemp+1: ntot))

       logdet = 0.d0
       do j=1,2*npol
          logdet =logdet +log(Spol(j,j))
       enddo

       logdet = 2.d0*logdet
    else
       argexp = 1.d37
       logdet = 1.d37
    endif

    alike = -.5d0*(argexp+logdet)
   
  end subroutine get_pix_loglike
  
  subroutine get_rotation_angle(r1,r2,a12,a21)
    !computes TWICE the rotation angle
    
    implicit none
    real(dp),intent(in),dimension(3) :: r1,r2
    real(dp),intent(out) :: a12,a21
    
    real(dp),parameter ::eps = 3.141592653589793d0/180.d0/3600.d0/100.d0
    real(dp),parameter,dimension(3):: zz =(/0,0,1/),epsilon =(/eps,0.d0,0.d0/)
    
    real(dp),dimension(3) :: r12,r1star,r2star
    real(dp) :: modv
    
    call ext_prod(r1,r2,r12)
    modv = sqrt(sum(r12*r12))
    if(modv.lt.1.d-8) then !same or opposite pixels
       a12 = 0.d0
       a21 = 0.d0
       return
       
    end if
    r12 = r12/modv
    
    call ext_prod(zz,r1,r1star)
    r1star(3) = 0.d0
    modv = sqrt(sum(r1star*r1star))
    if(modv.lt.1.d-8) then   !r1 is at a pole            
       r1star = r1star+epsilon
       modv = sqrt(sum(r1star*r1star))
    end if
    r1star = r1star/modv
    
    call ext_prod(zz,r2,r2star)
    r2star(3) = 0.d0
    modv = sqrt(sum(r2star*r2star))
    if(modv.lt.1.d-8) then   !r2 is at a pole            
       r2star = r2star+epsilon
       modv = sqrt(sum(r2star*r2star))
    end if
    r2star = r2star/modv
    
    modv = sum(r12*r1star)
    modv = min(1.d0,modv)
    modv = max(-1.d0,modv)
    if(sum(r12*zz).gt.0.d0) then
       a12 = 2.d0*acos(modv)
    else
       a12 = -2.d0*acos(modv)
    end if
    
    r12 = -r12   !r21 = - r12
    modv = sum(r12*r2star)
    modv = min(1.d0,modv)
    modv = max(-1.d0,modv)
    if(sum(r12*zz).gt.0.d0) then
       a21 = 2.d0*acos(modv)
    else
       a21 = -2.d0*acos(modv)
    end if
    
  end subroutine get_rotation_angle
  
  
  subroutine ext_prod(v1,v2,v3)
    implicit none
    real(dp),intent(in) :: v1(3),v2(3)
    real(dp),intent(out) :: v3(3)
    
    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

  end subroutine ext_prod
  
end module bflike_QU
