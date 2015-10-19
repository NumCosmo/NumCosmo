! ===========================================================================
module planck_teeebb_lowl

! This code is taken from the WMAP likelihood
! This code calculates the likelihood function of Q/U 
! polarization for EE signal plus noise. 
!
! ===========================================================================

  implicit none
  private
  public :: teeebb_lowl_like_setup, teeebb_lowl_likelihood, teeebb_pixlike_dof

  integer :: nsmax,nlmax,np,mp
#ifdef OPTIMIZE
  REAL(8), allocatable, dimension(:) :: xxx, yyy
  REAL(8), allocatable, dimension(:,:) :: ninvplninv2
  REAL, allocatable, dimension(:,:,:), target :: ninvplninv3
#else
  REAL, allocatable, dimension(:,:,:,:), Target :: ninvplninv
#endif
  INTEGER, allocatable, dimension(:) :: ngood
  REAL(8), allocatable, dimension(:,:) :: Dp0
  REAL(8), allocatable, dimension(:) :: m_r3,w_r3,p_r3,f_r3,zzz
  COMPLEX, allocatable, dimension(:,:), Target :: alm_tt,NinvY
  REAL, allocatable, dimension(:) :: wl
contains

!===========================================
subroutine teeebb_lowl_like_setup
!===========================================

  use planck_options, only : Planck_data_dir, get_free_lun, use_wmap_pol
  implicit none

  character(LEN=2) :: rlz
  integer,parameter :: ires = 3
  REAL, dimension(:),allocatable :: T,N,Mask_R3
REAL, dimension(:,:), allocatable :: inmap
  CHARACTER(len=3) :: da(10),sband
  INTEGER :: ReadStatus
  CHARACTER(len=256) :: qfile,ufile,maskfile,filename(0:9),eebbdir

  INTEGER :: dum,nsmax,nlmax,l,ip,jp,m,lun,stat,i,j,k
  REAL, allocatable, dimension(:,:), Target :: NinvQUr3
  logical ::yes

  Complex (Kind=4), Dimension(:,:),   Pointer :: Cptr2
  Real (Kind=4),    Dimension(:,:),   Pointer :: Rptr2
  Real (Kind=4),    Dimension(:,:,:), Pointer :: Rptr3

  INCLUDE 'read_archive_map.fh'
  Include 'read_fits.fh'

!------------------------------
! Set initial parameters
!-----------------------------

  da=(/'K1','Ka','Q1','Q2','V1','V2','W1','W2','W3','W4'/)
  da(2)='Ka1'


  nsmax = 2**ires
  nlmax = 3*nsmax-1 ! use Nyquist sampling rather than 2*nsmax
  np = 12*nsmax**2

 if(use_wmap_pol) then
   eebbdir = trim(Planck_data_dir)//'lowlP/wmap9/'
   filename(0)=trim(eebbdir)//'masked_ee_ninvplninv_qu_r3_corrected_9yr.KaQV.fits'
   filename(1)=trim(eebbdir)//'masked_bb_ninvplninv_qu_r3_corrected_9yr.KaQV.fits'
   filename(2)=trim(eebbdir)//'masked_ninv_qu_r3_corrected_9yr.KaQV.fits'
   filename(3)=trim(eebbdir)//'wt_r3_9yr.KaQV.map_q'
   filename(4)=trim(eebbdir)//'wt_r3_9yr.KaQV.map_u'
   filename(6)=trim(eebbdir)//'masked_ninvy_e_qu_r3_corrected_9yr.KaQV.fits'
else
  eebbdir = trim(Planck_data_dir)//'lowlP/planck/'
   filename(0)=trim(eebbdir)//'masked_ee_ninvplninv_qu_r3.fits'
   filename(1)=trim(eebbdir)//'masked_bb_ninvplninv_qu_r3.fits'
   filename(2)=trim(eebbdir)//'masked_ninv_qu_r3_KaQV.fits'
   filename(3)=trim(eebbdir)//'wt_r3_KaQV_map_q.fits'
   filename(4)=trim(eebbdir)//'wt_r3_KaQV_map_u.fits'
   filename(6)=trim(eebbdir)//'masked_ninvy_e_qu_r3.fits'
endif

  filename(5)=trim(Planck_data_dir)//'lowlP/alms_dx9_cmb_single_sample_nopixwin.dat'
!  filename(5)=trim(Planck_data_dir)//'lowlP/alm_tt_commander.dat'
!  filename(5)=trim(Planck_data_dir)//'lowlP/alm_tt_fs_r9_ilc_nopixwin_7yr.dat'
  filename(9)=trim(Planck_data_dir)//'healpix_data/pixel_window_n0008.txt'
 
!------------------------------
! read in res3 mask
!------------------------------

  allocate(mask_r3(0:np-1))
  allocate(T(0:np-1))
  allocate(N(0:np-1))
  allocate(inmap(0:np-1,1))

  maskfile=trim(Planck_data_dir)//'lowlP/mask_r3.fits'
  inquire(file= maskfile,exist = yes)
  if(.not.yes)then
     write(*,*)"teeebb maskfile not found",trim(maskfile)
  endif

  CALL READ_ARCHIVE_MAP(maskfile,T,Mask_R3,dum,ReadStatus)
  if (ReadStatus/=0) then
     print*,'unable to read in mask',trim(maskfile)
     stop
  endif

  ALLOCATE(ngood(0:np-1))
  mp = 0
  do ip=0,np-1
     if (Mask_R3(ip)/=0) then
        ngood(mp) = ip

        mp = mp + 1

     endif
  enddo

  !------------------------------
  ! read in N^{-1}P_l N^{-1} at res3
  !------------------------------

#ifdef OPTIMIZE

  ALLOCATE( ninvplninv2(mp*(2*mp+1),2*(nlmax-1)), stat=stat )
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for ninvplninv2'
    Stop
  End If

  ALLOCATE(ninvplninv3(0:nlmax,0:2*np-1,0:2*np-1), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for ninvplninv3'
    Stop
  End If
  Rptr3 => ninvplninv3(:,:,:)

  Call Read_FITS_Real_3D (filename(0), Rptr3, stat)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(0))
    Stop
  End If
	k = 1
! --- original
!	do i = 0,2*mp-1
!	do j = i,2*mp-1
! --- RF	
        do j = 0,2*mp-1
           do i = 0,j
! ---
		ip = ngood(mod(i,mp)) + np*(i/mp)
		jp = ngood(mod(j,mp)) + np*(j/mp)

		ninvplninv2(k,1:nlmax-1) = ninvplninv3(2:nlmax,ip,jp)
		k = k + 1

	end do
	end do

  Call Read_FITS_Real_3D (filename(1), Rptr3, stat)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(1))
    Stop
  End If

	k = 1
! --- original
!	do i = 0,2*mp-1
!	do j = i,2*mp-1
! --- RF
        do j = 0,2*mp-1
           do i = 0,j
! ---
		ip = ngood(mod(i,mp)) + np*(i/mp)
		jp = ngood(mod(j,mp)) + np*(j/mp)

		ninvplninv2(k,nlmax:2*(nlmax-1)) = ninvplninv3(2:nlmax,ip,jp)
		k = k + 1

	end do
	end do

	deallocate( ninvplninv3 )

	allocate( xxx(2*(nlmax-1)) )
	allocate( yyy(mp*(2*mp+1)) )

#else

  ALLOCATE(ninvplninv(0:nlmax,0:2*np-1,0:2*np-1,2), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for ninvplninv'
    Stop
  End If

  Rptr3 => ninvplninv(:,:,:,1)
  Call Read_FITS_Real_3D (filename(0), Rptr3, stat)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(0))
    Stop
  End If

  Rptr3 => ninvplninv(:,:,:,2)
  Call Read_FITS_Real_3D (filename(1), Rptr3, stat)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(1))
    Stop
  End If

#endif

  !------------------------------
  ! read in N^{-1} at res3
  !------------------------------

  ALLOCATE(NinvQUr3(0:1535,0:1535), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for NinvQUr3'
    Stop
  End If
	
  Rptr2 => NinvQUr3
  Call Read_FITS_Real_2D (filename(2), Rptr2, stat)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(2))
    Stop
  End If

  !------------------------------
  ! read in maps at res3
  !------------------------------

  ALLOCATE(w_r3(0:2*mp-1),m_r3(0:2*mp-1),p_r3(0:2*mp-1), stat=stat)
  allocate( zzz(0:2*mp-1) )
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for w_r3, m_r3, or p_r3'
    Stop
  End If
  qfile=filename(3)
  CALL READ_ARCHIVE_MAP(qfile,T,N,np,ReadStatus)

  do ip=0,mp-1
     w_r3(ip) = T(ngood(ip))*Mask_R3(ngood(ip))
  enddo

  ufile=filename(4)
  CALL READ_ARCHIVE_MAP(ufile,T,N,np,ReadStatus)

  do ip=0,mp-1
     w_r3(mp+ip) = T(ngood(ip))*Mask_R3(ngood(ip))
  enddo

  ALLOCATE(Dp0(0:2*mp-1,0:2*mp-1), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for Dp0'
    Stop
  End If
! --- original
!  do ip=0,mp-1
!     do jp=0,mp-1
! --- RF
  do jp=0,mp-1
     do ip=0,mp-1
! ---
        Dp0(ip,jp) = NinvQUr3(ngood(ip),ngood(jp))
        Dp0(ip,mp+jp) = NinvQUr3(ngood(ip),np+ngood(jp))
        Dp0(mp+ip,jp) = NinvQUr3(np+ngood(ip),ngood(jp))
        Dp0(mp+ip,mp+jp) = NinvQUr3(np+ngood(ip),np+ngood(jp))
     enddo
  enddo

  !------------------------------
  ! read in alm_tt 
  !------------------------------
  ALLOCATE(alm_tt(0:nlmax,0:nlmax), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for alm_tt'
    Stop
  End If
  call get_free_lun( lun )
  open(lun,file=filename(5),status='old')
  do l=0,nlmax
     do m=0,l
        read(lun,*)alm_tt(l,m)
     enddo
  enddo
  close(lun)

  !------------------------------
  ! read in Ninv Y
  !------------------------------
  ALLOCATE(NinvY(0:1535,300), stat=stat)
  If (stat .NE. 0) Then
    Print *, 'Memory allocation error for NinvY'
    Stop
  End If

  Cptr2 => NinvY
  Call Read_FITS_Complex_2D_LM (filename(6), Cptr2, stat, IndFmt=2)
  If (stat .NE. 0) Then
    Print *, 'Error ', stat, ' while reading ', trim(filename(6))
    Stop
  End If

!
! read in pixel window
!
  ALLOCATE(wl(0:nlmax))
  call get_free_lun( lun )
  OPEN(lun,file=filename(9),status='old')
  do l = 0,nlmax
     read(lun,*) wl(l)
  end do
  CLOSE(lun)

  deallocate(T,N,Mask_R3)
  DEALLOCATE(NinvQUr3)

end subroutine teeebb_lowl_like_setup


  function teeebb_pixlike_dof()
	integer :: teeebb_pixlike_dof
	teeebb_pixlike_dof = 2*mp
  end function

!===========================================
subroutine teeebb_lowl_likelihood(nlmaxin,clttin,cltein,cleein,clbbin,chisq_r3,lndet)
!=======================================

  use planck_options, only : teeebb_pixlike_lndet_offset

  IMPLICIT NONE

  REAL(8), dimension(2:*),intent(in) :: clttin,cltein,cleein,clbbin
  real(8),intent(out) ::lndet,chisq_r3 ! like
  integer, intent(in) ::nlmaxin
#ifndef OPTIMIZE
  REAL(8), allocatable, dimension(:,:) :: NinvSNinv
  real(8),allocatable,dimension(:,:) ::CQU
#endif
  REAL(8), external :: DDOT
  REAL(8), allocatable, dimension(:,:) :: Dp 
  real(8),allocatable,dimension(:) ::clee,clbb,cltt,clte
  INTEGER :: ip,jp,nlmax,info,l,m,i,k
  integer :: t_start,t_end,crate,cmax
  REAL :: Omega_pix

  nlmax = 23
  Omega_pix = 3.14159/(3.*8.**2.)
  if(nlmaxin.ne.nlmax)then
     write(*,*)"need nlmax",nlmax,", in teeebb likelihood, currently",nlmaxin
     stop
  endif

  allocate( Cltt(0:nlmax),Clte(0:nlmax),Clee(0:nlmax),Clbb(0:nlmax) )
  Cltt(0:1)=0.
  Clte(0:1)=0.
  Clee(0:1)=0.
  Clbb(0:1)=0.
  do l=2,nlmax
     Cltt(l) = Clttin(l)/dble(l*(l+1))*(2.*3.14159)*1.d-6*wl(l)**2.
     Clte(l) = Cltein(l)/dble(l*(l+1))*(2.*3.14159)*1.d-6*wl(l)**2.
     Clee(l) = Cleein(l)/dble(l*(l+1))*(2.*3.14159)*1.d-6*wl(l)**2.
     Clbb(l) = Clbbin(l)/dble(l*(l+1))*(2.*3.14159)*1.d-6*wl(l)**2.
  enddo

  !------------------------------
  ! compute [N^{-1}S N^{-1} + N^{-1}]^{-1}
  !------------------------------
  ALLOCATE(Dp(0:2*mp-1,0:2*mp-1))

  lndet = 0d0

#ifdef OPTIMIZE
		!! MRN

	do l = 2,nlmax
		xxx(l-1) = clee(l)-clte(l)**2./cltt(l)
		xxx(l+nlmax-2) = clbb(l)
	end do

	call DGEMV( 'N', mp*(2*mp+1), 2*(nlmax-1), 1.d0, ninvplninv2, &
		mp*(2*mp+1), xxx, 1, 0.d0, yyy, 1 )

	k = 1
! --- original ---
!	do ip = 0,2*mp-1
!	do jp = ip,2*mp-1
! --- RF
	do jp = 0,2*mp-1
	do ip = 0,jp
! ---
		Dp(ip,jp) = Dp0(ip,jp) &
                        + yyy(k)
		k = k + 1
	end do
	end do

#else

!
! fill only the upper triangular part of N^{-1}SN^{-1}.
! add the foreground error term.
!
  ALLOCATE(NinvSNinv(0:2*mp-1,0:2*mp-1))
  ALLOCATE(CQU(0:2*mp-1,0:2*mp-1))
  NinvSNinv = 0.
! --- original
!  do ip=0,mp-1   
!     do jp=ip,mp-1
! --- RF
  do jp=0,mp-1
     do ip=0,jp   
! ---
        NinvSNinv(ip,jp) = &
             +DDOT(nlmax-1,clee(2:)-clte(2:)**2./cltt(2:),1,DBLE(ninvplninv(2:,ngood(ip),ngood(jp),1)),1) &
             +DDOT(nlmax-1,clbb(2:),1,DBLE(ninvplninv(2:,ngood(ip),ngood(jp),2)),1)
        
        NinvSNinv(ip,mp+jp) = &
             +DDOT(nlmax-1,clee(2:)-clte(2:)**2./cltt(2:),1,DBLE(ninvplninv(2:,ngood(ip),np+ngood(jp),1)),1)&
             +DDOT(nlmax-1,clbb(2:),1,DBLE(ninvplninv(2:,ngood(ip),np+ngood(jp),2)),1)

        NinvSNinv(mp+ip,mp+jp) = &
             +DDOT(nlmax-1,clee(2:)-clte(2:)**2./cltt(2:),1,DBLE(ninvplninv(2:,np+ngood(ip),np+ngood(jp),1)),1) &
             +DDOT(nlmax-1,clbb(2:),1,DBLE(ninvplninv(2:,np+ngood(ip),np+ngood(jp),2)),1)
     enddo
! --- original
!   do jp=0,ip-1   
! --- RF 
  end do   

  do jp=0,mp-1
     do ip=jp+1,mp-1
!--- 
        NinvSNinv(ip,mp+jp) = NinvSNinv(ip,mp+jp) &
             +DDOT(nlmax-1,clee(2:)-clte(2:)**2./cltt(2:),1,DBLE(ninvplninv(2:,ngood(ip),np+ngood(jp),1)),1)&
             +DDOT(nlmax-1,clbb(2:),1,DBLE(ninvplninv(2:,ngood(ip),np+ngood(jp),2)),1)
     enddo
    
  enddo

  Dp = Dp0 + NinvSNinv 

  DEALLOCATE(NinvSNinv)

#endif

  CALL DPOTRF('U',2*mp,Dp,2*mp,info) 
 
  IF(info.NE.0)then
        write(*,*) 'error' !	call wmap_likelihood_error( 'teeebb: bad dpotrf', info )
	chisq_r3 = 1d30
	lndet = 0d0
	return
  endif

  do ip=0,2*mp-1
!     write(*,*) ip, lndet
     lndet = lndet + 2.*log(real(Dp(ip,ip)))
  enddo

#ifndef OPTIMIZE
   CALL DPOTRI('U',2*mp,Dp,2*mp,info)
#endif

  IF(info.NE.0)then
!	call wmap_likelihood_error( 'teeebb: bad dpotri', info )
	write(*,*) 'Bad'
	chisq_r3 = 1d30
	lndet = 0d0
	return
  endif

#ifndef OPTIMIZE
  CQU(:,:) = Dp 
  DEALLOCATE(Dp)
#endif


  !-----------------------------------
  ! calculate the predicted QU at res3
  !-----------------------------------
  p_r3 = 0.

  do ip=0,mp-1
     i = 3
     do l=2,23
        i = i+1
        p_r3(ip)    = p_r3(ip)    &
		+ clte(l)/cltt(l)*wl(l)*alm_tt(l,0)*NinvY(ngood(ip),i)
        p_r3(ip+mp) = p_r3(ip+mp) &
		+ clte(l)/cltt(l)*wl(l)*alm_tt(l,0)*NinvY(ngood(ip)+np,i)
        do m=1,l
           i = i+1
           p_r3(ip) = p_r3(ip)&
		+ clte(l)/cltt(l)*wl(l)*( alm_tt(l,m)*NinvY(ngood(ip),i)&
		                  +conjg(alm_tt(l,m)*NinvY(ngood(ip),i)) )
           p_r3(ip+mp) = p_r3(ip+mp)&
                + clte(l)/cltt(l)*wl(l)*( alm_tt(l,m)*NinvY(ngood(ip)+np,i)&
		                  +conjg(alm_tt(l,m)*NinvY(ngood(ip)+np,i)) )
        enddo
     enddo 
  enddo

#ifdef OPTIMIZE

	m_r3 = w_r3-p_r3
	call DPOTRS( 'U', 2*mp, 1, Dp, 2*mp, m_r3, 2*mp, info )
	!chisq_r3 = DDOT(2*mp,m_r3,1,w_r3-p_r3,1) 
	zzz = w_r3-p_r3
	chisq_r3 = sum( m_r3*zzz )
#else

  CALL DSYMV('U',2*mp,1.d0,CQU(:,:),2*mp,w_r3-p_r3,1,0.d0,m_r3,1)
  chisq_r3 = DDOT(2*mp,m_r3,1,w_r3-p_r3,1) 

  DEALLOCATE(CQU)

#endif

  chisq_r3 = chisq_r3/2d0
  lndet = (lndet - teeebb_pixlike_lndet_offset)/2d0
  
END subroutine TEEEBB_LOWL_LIKELIHOOD

end module Planck_teeebb_lowl
