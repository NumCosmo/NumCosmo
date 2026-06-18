module temp_like

  implicit none

  real*8, dimension(:), allocatable :: X
  real*8,  dimension(:,:), allocatable :: c_inv
  integer :: Nspec,nX,num_ells,nXfromdiff
  integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
  real*8 :: sz_143_temp(0:5000)
  real*8 :: ksz_temp(0:5000), tszxcib_temp(0:5000)
  integer :: lmax_sz
  real*8, dimension(:,:), allocatable :: beam_cov_inv,beam_cov_full
  real*8, dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
  integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim
  real*8,dimension(:,:),allocatable :: Cl_fg
  real(8),dimension(:,:),allocatable::gcal
  integer:: num_non_beam=14
  real(8), dimension(:),  allocatable::  X_beam_corr_model, Y
  logical:: has_dust, has_calib_prior
  logical :: needinit=.true.

  real(8), parameter :: sz_bandpass100_nom143 = 2.022d0
  real(8), parameter :: cib_bandpass143_nom143 = 1.134d0
  real(8), parameter :: sz_bandpass143_nom143 = 0.95d0
  real(8), parameter :: cib_bandpass217_nom217 = 1.33d0


  integer, allocatable :: marge_indices(:),marge_indices_reverse(:)
  integer, allocatable :: keep_indices(:),keep_indices_reverse(:)
  real(8), allocatable :: beam_conditional_mean(:,:)
  logical, allocatable :: want_marge(:)
  integer marge_num, keep_num
  real(8) :: beam_factor = 1 !new beam modes are already scaled (2.7_campc)
  logical :: make_cov_marged = .false.

contains



!!  subroutine like_init(like_file, sz143_file, tszxcib_file, ksz_file, beam_file)
!!
!!    real*8, dimension(:), allocatable :: X
!!    real*8,  dimension(:,:), allocatable :: c_inv
!!    integer :: Nspec,nX,num_ells,nXfromdiff
!!    integer, dimension(:), allocatable  :: lminX, lmaxX, np, npt
!!    real*8 :: sz_143_temp(0:5000)
!!    real*8 :: ksz_temp(0:5000), tszxcib_temp(0:5000)
!!    integer :: lmax_sz
!!    real*8, dimension(:,:), allocatable :: beam_cov_inv
!!    real*8, dimension(:,:,:), allocatable :: beam_modes ! mode#, l, spec#
!!    integer :: num_modes_per_beam,beam_lmax,beam_Nspec,cov_dim
!!  
!!    integer :: i, j, l,dummy
!!
!!    character*100 like_file, sz143_file, ksz_file, tszxcib_file, beam_file
!!
!!    ! cl_ksz_148_tbo.dat file is in D_l, format l D_l, from l=2 to 10000
!!    ! tsz_x_cib_template.txt is is (l D_l), from l=2 to 9999, normalized to unity 
!!    !    at l=3000 
!!
!!    if(needinit .eqv. .false.) then
!!       return
!!    endif
!!
!!    open(48, file=like_file, form='unformatted', status='unknown')
!!
!!    read(48) Nspec,nX
!!    allocate(lminX(Nspec))
!!    allocate(lmaxX(Nspec))
!!    allocate(np(Nspec))
!!    allocate(npt(Nspec))
!!    allocate(X(nX))
!!    allocate(c_inv(nX,nX))
!!
!!    read(48) (lminX(i), lmaxX(i), np(i), npt(i), i = 1, Nspec)
!!    read(48) (X(i), i=1, nX)
!!    read(48) 
!!    read(48) ((c_inv(i, j), j = 1, nX), i = 1,  nX)
!!    close(48)
!!
!!    !  open(48, file=sz100_file, form='unformatted', status='unknown')
!!    !  read(48) lmax_sz
!!    !  read(48) (sz_100_temp(l), l = 0, lmax_sz)
!!    !  close(48)
!!
!!    lmax_sz=5000
!!
!!    open(48, file=sz143_file, form='formatted', status='unknown')
!!    do i=2,lmax_sz
!!       read(48,*) dummy,sz_143_temp(i)
!!    enddo
!!    close(48)
!!
!!    open(48, file=ksz_file, form='formatted',status='unknown')
!!    do i=2,lmax_sz
!!       read(48,*) dummy,ksz_temp(i)
!!    enddo
!!    close(48)
!!
!!    open(48, file=tszxcib_file,form='formatted',status='unknown')
!!    do i=2,lmax_sz
!!       read(48,*) dummy,tszxcib_temp(i)
!!    enddo
!!    close(48)
!!
!!    open(48, file=beam_file, form='unformatted', status='unknown')
!!    read(48) beam_Nspec,num_modes_per_beam,beam_lmax
!!    if(beam_Nspec.ne.Nspec) stop 'Problem: beam_Nspec != Nspec'
!!    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
!!    cov_dim=beam_Nspec*num_modes_per_beam
!!    allocate(beam_cov_inv(cov_dim,cov_dim))
!!    read(48) (((beam_modes(i,l,j),j=1,Nspec),l=0,beam_lmax),i=1,num_modes_per_beam)
!!    read(48) ! skipping beam_cov
!!    read(48) ((beam_cov_inv(i,j),j=1,cov_dim),i=1,cov_dim)
!!    close(48)
!!    call like_init_frommem(Nspec, nX,lminX,lmaxX,np,npt, c_inv,X,lmax_sz,sz_143_temp,ksz_temp,tszxcib_temp,beam_Nspec,num_modes_per_beam,beam_lmax,cov_dim,beam_cov_inv,beam_modes,0,1)
!!
!!    deallocate(lminX)
!!    deallocate(lmaxX)
!!    deallocate(np)
!!    deallocate(npt)
!!    deallocate(X)
!!    deallocate(c_inv)
!!    deallocate(beam_cov_inv)
!!    deallocate(beam_modes)
!!
!!    return
!!  end subroutine like_init

  subroutine like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes,dust_flag,calib_flag,marge_flags,marge_mode)
    integer,intent(in)::iNspec,inX,inum_modes_per_beam,ibeam_lmax,ibeam_Nspec,icov_dim,ilmax_sz
    integer,dimension(:)::ilminX,ilmaxX,inp,inpt
    real*8, dimension(:) :: iX
    real*8,  dimension(:,:) ::ic_inv
    real*8,intent(in) :: iksz_temp(:), itszxcib_temp(:)
    real*8,intent(in) :: isz_143_temp(:)
    real*8, dimension(:,:) :: ibeam_cov_inv
    real*8, dimension(:,:,:) :: ibeam_modes ! mode#, l, spec#
    real(8)::renorm    
    integer,intent(in)::dust_flag,calib_flag
    logical,intent(in),dimension(:)::marge_flags
    real(8),intent(in),dimension(:,:)::marge_mode
    integer::i,j,k

    Nspec = iNspec
    nX = inX

    allocate(lminX(Nspec))
    allocate(lmaxX(Nspec))
    allocate(np(Nspec))
    allocate(npt(Nspec))
    allocate(X(nX))
    allocate(c_inv(nX,nX))

    lminX = ilminX
    lmaxX = ilmaxX
    X = iX
    c_inv = ic_inv
    np = inp
    npt = inpt
    
    lmax_sz = ilmax_sz
    
    ksz_temp = iksz_temp
    renorm=1.d0/ksz_temp(3000)
    ksz_temp=ksz_temp*renorm
    
    sz_143_temp = isz_143_temp
    renorm=1.d0/sz_143_temp(3000)
    sz_143_temp=sz_143_temp*renorm

    tszxcib_temp = itszxcib_temp
    renorm=1.d0/tszxcib_temp(3000)
    tszxcib_temp=tszxcib_temp*renorm

    beam_Nspec = ibeam_Nspec
    num_modes_per_beam = inum_modes_per_beam
    beam_lmax = ibeam_lmax
    cov_dim = icov_dim

    allocate(beam_modes(num_modes_per_beam,0:beam_lmax,beam_Nspec))
    allocate(beam_cov_inv(cov_dim,cov_dim))

    beam_modes = ibeam_modes
    beam_cov_inv = ibeam_cov_inv

    allocate(Cl_fg(Nspec,0:maxval(lmaxX)))
    allocate(gcal(1:Nspec,0:beam_lmax))
    allocate(X_beam_corr_model(1:nX))
    allocate(Y(1:nX))

    num_non_beam=14
    if (dust_flag .eq. 1) then
      num_non_beam = 15
      has_dust = .true.
    endif
    has_calib_prior = .true.
    if (calib_flag .eq. 0) then
      has_calib_prior = .false.
    endif

    allocate(want_marge(cov_dim))
    do i=1,cov_dim
      want_marge(i) =  marge_flags(i)
    enddo
    
    marge_num=count(want_marge)
    keep_num=cov_dim-marge_num

    
    allocate(marge_indices(marge_num))
    allocate(marge_indices_reverse(cov_dim))
    allocate(keep_indices(keep_num))
    allocate(keep_indices_reverse(cov_dim))

    j=0
    k=0
    do i=1,cov_dim
        if (want_marge(i)) then
            j=j+1
            marge_indices(j) = i
        else
            k=k+1
            keep_indices(k)=i
        end if
        marge_indices_reverse(i)=j
        keep_indices_reverse(i)=k
    end do

    
    if (marge_num.ne.0) then
      allocate(beam_conditional_mean(marge_num, keep_num))
      do i=1, marge_num
        beam_conditional_mean(i,:) = marge_mode(i,:)
      end do
      if (keep_num>0) then
        allocate(beam_cov_full(cov_dim,cov_dim))
        beam_cov_full = beam_cov_inv
        call Matrix_inverse_internal(beam_cov_full)
        deallocate(beam_cov_inv)
        allocate(beam_cov_inv(keep_num,keep_num))
        beam_cov_inv = beam_cov_full(keep_indices,keep_indices)
        call Matrix_inverse_internal(beam_cov_inv)
        deallocate(beam_cov_full)
      end if

    endif

    needinit=.false.

  end subroutine like_init_frommem


  subroutine compute_fg(Cl_fg,freq_params)
    real(8), intent(in)  :: freq_params(:)
    real(8), intent(inout)  :: Cl_fg(1:,0:)
    integer :: i, j, l, ii,jj
    real(8):: A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, xi, A_ksz, ncib,A_dust
    real(8):: zell, zCIB
    real(8):: Lfac
    integer::lmax

    
    if(needinit .eqv. .true.) then
       stop "initialize first"
    endif
    if (size(freq_params) < num_non_beam) stop 'CAMspec: not enough parameters'
    A_ps_100=freq_params(1)
    A_ps_143 = freq_params(2)
    A_ps_217 = freq_params(3)
    A_cib_143=freq_params(4)
    A_cib_217=freq_params(5) 
    A_sz =freq_params(6)  !143
    r_ps = freq_params(7)
    r_cib=freq_params(8)
    ncib = freq_params(9)
    xi=freq_params(13)
    A_ksz=freq_params(14)
    

    if(has_dust) then
      A_dust = freq_params(15)
    endif

    lmax = ubound(Cl_fg,2)
    !100x100
    
    Cl_fg(:,0) = 0
    do l = 1, lmax
      zell = real(l,8)
      Lfac = real(l*(l+1),8)
      
      !100x100
      Cl_fg(1,l) = A_ps_100*1.d-6/9.d0 + &
        A_ksz*ksz_temp(l)/Lfac+ &
        A_sz*sz_bandpass100_nom143*sz_143_temp(l)/Lfac
      
      ! 143x143
      zCIB = cib_bandpass143_nom143*A_cib_143*(zell/3000)**(ncib)/Lfac
      Cl_fg(2,l) = A_ps_143*1.d-6/9.d0 + zCIB + &
        A_ksz*ksz_temp(l)/Lfac + A_sz*sz_bandpass143_nom143*sz_143_temp(l)/Lfac  &
        -2.0*sqrt(cib_bandpass143_nom143*A_cib_143*sz_bandpass143_nom143*A_sz)*xi*tszxcib_temp(l)/Lfac
      
      ! 217x217
      zCIB = cib_bandpass217_nom217*A_cib_217*(zell/3000)**(ncib)/Lfac
      Cl_fg(3,l) = A_ps_217*1.d-6/9.d0 + zCIB &
        + A_ksz*ksz_temp(l)/Lfac
      
      ! 143X217
      zCIB = dsqrt(cib_bandpass143_nom143*A_cib_143*cib_bandpass217_nom217*A_cib_217)*(zell/3000)**(ncib) &
        /Lfac
      Cl_fg(4,l) = r_ps*dsqrt(A_ps_143*A_ps_217)*1.d-6/9.d0 + r_cib*zCIB &
        +A_ksz*ksz_temp(l)/Lfac  &
        -sqrt(cib_bandpass217_nom217*A_cib_217*sz_bandpass143_nom143*A_sz)*xi*tszxcib_temp(l)/Lfac
    end do    
    if(has_dust) then
      A_dust = freq_params(15)
      call add_dust(Cl_fg(1,:), Cl_fg(2,:), Cl_fg(3,:), Cl_fg(4,:), 1,lmax,A_dust)
    endif
  end subroutine
      
  subroutine beams_and_cal_prior(logprior, freq_params)
    real(8),intent(out)::logprior
    real(8),intent(in)::freq_params(:)
    real(8)::cal0,cal1,cal2
    integer::i,j,k,l,ii,jj
    real(8) beam_coeffs(Nspec,num_modes_per_beam)
    real(8) beam_params(cov_dim)

    if(needinit .eqv. .true.) then
       stop "initialize first"
    endif

    if (size(freq_params) < num_non_beam +  keep_num) stop 'CAMspec: not enough parameters'

    cal0=freq_params(10)
    cal1=freq_params(11) 
    cal2=freq_params(12)
    
    call fill_beam_params(freq_params,beam_params)

    do ii=1,beam_Nspec
      do jj=1,num_modes_per_beam
          beam_coeffs(ii,jj)=beam_params(jj+num_modes_per_beam*(ii-1))
      enddo
    enddo

    logprior = 0

    ! add prior on the beam coefs
    if (keep_num>0) then
      do i=1,beam_Nspec
        do j=1,beam_Nspec
          do k=1,num_modes_per_beam
            jj = k+num_modes_per_beam*(i-1)
            if (.not. want_marge(jj)) then
              do l=1,num_modes_per_beam
                ii = l+num_modes_per_beam*(j-1)
                if (.not. want_marge(ii)) then
                  logprior = logprior + beam_coeffs(j,l)*beam_cov_inv(keep_indices_reverse(ii),keep_indices_reverse(jj))*beam_coeffs(i,k)
                endif
              enddo
            endif
          enddo
        enddo
      enddo
    endif

  

    ! add prior on the calibration coefs

    if (has_calib_prior) then
      logprior = logprior + ((cal2/cal1-0.9966d0)/0.0015d0)**2  + ((cal0/cal1-1.0006d0)/0.0004d0)**2
    endif

  end subroutine

   subroutine plik_dust_template(dust_template__l, lmin, lmax)

    integer,                                intent(in)    :: lmin
    integer,                                intent(in)    :: lmax
    real(8),      dimension(lmin:),         intent(out)   :: dust_template__l

    integer                                               :: i
    real(8)                                               :: l_pivot
    real(8)                                               :: gal_index
    real(8),                                parameter     :: twopi = 6.283185307179586476925286766559005768394

    l_pivot   = 500.0
    gal_index = -2.6

    do i = lmin,lmax
       dust_template__l(i) = (dble(i)/l_pivot)**gal_index / twopi
    enddo


  end subroutine plik_dust_template



  subroutine add_dust(cl100, cl143, cl217, cl143x217, lmin,lmax, amplitude)

    integer,                                intent(in)    :: lmin
    integer,                                intent(in)    :: lmax
    real(8),                                intent(in)    :: amplitude
    real(8),      dimension(lmin:),      intent(inout) :: cl100
    real(8),      dimension(lmin:),      intent(inout) :: cl143
    real(8),      dimension(lmin:),      intent(inout) :: cl217
    real(8),      dimension(lmin:),  intent(inout) :: cl143x217

    real(8)                                               :: cc_100
    real(8)                                               :: cc_143
    real(8)                                               :: cc_217
    real(8)                                               :: rescaling_143_to_100
    real(8)                                               :: rescaling_143_to_217
    real(8),      dimension(:),             allocatable   :: dust_template__l

! frequency averaged color corrections
    cc_100 = 1.090
    cc_143 = 1.023
    cc_217 = 1.130

! amplitude rescalings between nominal frequencies
! include mask 1 -> mask 3 conversion factor for 100 ghz
    rescaling_143_to_100 = 0.466 * 1.66
    rescaling_143_to_217 = 3.170

    allocate(dust_template__l(lmin:lmax))
    dust_template__l = 0.0

    call plik_dust_template(dust_template__l, lmin, lmax)

    cl100 = cl100 + amplitude * cc_100**2 * rescaling_143_to_100**2&
         & * dust_template__l(lmin:lmax)

    cl143 = cl143 + amplitude * cc_143**2&
         & * dust_template__l(lmin:lmax)

    cl217 = cl217 + amplitude * cc_217**2 * rescaling_143_to_217**2&
         & * dust_template__l(lmin:lmax)

    cl143x217 = cl143x217 + amplitude * cc_143 * cc_217 * rescaling_143_to_217&
         & * dust_template__l(lmin:lmax)

    deallocate(dust_template__l)


  end subroutine add_dust

  subroutine fill_beam_params(freq_params,beam_params)
    real(8),intent(in),dimension(:)::freq_params
    real(8),intent(out),dimension(:)::beam_params
    integer::i

    do i=1,keep_num
      beam_params(keep_indices(i)) = freq_params(num_non_beam+i)
    enddo

    if (marge_num>0) then
      beam_params(marge_indices) = matmul(beam_conditional_mean, freq_params(num_non_beam+1:))
    endif

  end subroutine  

  subroutine compute_beams_and_cal(gcal,freq_params)
    real(8),dimension(1:Nspec,0:beam_lmax)::gcal
    real(8), intent(in)  :: freq_params(:)
    real(8)::cal0, cal1, cal2
    real(8),dimension(Nspec)::cal
    real(8)::corrected_beam
    integer::i,ispec,l
    real(8),dimension(cov_dim)::beam_params

    if(needinit .eqv. .true.) then
       stop "initialize first"
    endif

    if (size(freq_params) < num_non_beam +  keep_num) stop 'CAMspec: not enough parameters'
    cal0=freq_params(10)
    cal1=freq_params(11) 
    cal2=freq_params(12)
    cal(1) = cal0
    cal(2) = cal1
    cal(3) = cal2
    cal(4) = dsqrt(cal1*cal2)

    call fill_beam_params(freq_params,beam_params)
    do ispec=1,Nspec
      do l = 0, beam_lmax
        corrected_beam=1.d0
        do i=1,num_modes_per_beam
          corrected_beam=corrected_beam + beam_params(i+num_modes_per_beam*(ispec-1)) * beam_modes(i,l,ispec) * beam_factor
        enddo
        gcal(ispec,l) = corrected_beam/cal(ispec)
      enddo
    enddo

  end subroutine
  
  subroutine calc_like(zlike,  cell_cmb, freq_params)
    real(8), intent(in)  :: freq_params(:)
    real(8), dimension(0:) :: cell_cmb
    integer :: i,j
    real(8) zlike,ztemp
    integer :: ofs,ispec 
    real(8) logprior
    
    if(needinit .eqv. .true.) then
       stop "initialize first"
    endif

    
    if(Nspec.ne.4) then
      print*, 'Nspec inconsistent with foreground corrections in calc_like.'
      stop
    end if

    ! compute forgrounds
    call compute_fg(Cl_fg,freq_params)

    !compute beam and cal corrections
    call compute_beams_and_cal(gcal,freq_params)
    
    !   100 foreground
    ofs = 1
    do ispec=1,Nspec
      !print *,nX,ispec,npt(ispec),ofs,lmaxX(ispec),lminX(ispec),npt(ispec)-ofs+lmaxX(ispec)-lminX(ispec)+1
      X_beam_corr_model(npt(ispec):npt(ispec)-ofs+lmaxX(ispec)-lminX(ispec)+1) = (cell_cmb(lminX(ispec):lmaxX(ispec)) + Cl_fg(ispec,lminX(ispec):lmaxX(ispec))) * gcal(ispec,lminX(ispec):lmaxX(ispec))
      ofs = 1
    enddo

    !X_beam_corr_model = X
    Y = X - X_beam_corr_model
    
    zlike = 0
    !$OMP parallel do private(j,i,ztemp) reduction(+:zlike) schedule(static,16)
    do  j = 1, nX
      ztemp= dot_product(Y(j+1:nX), c_inv(j+1:nX, j))
      zlike=zlike+ (ztemp*2 +c_inv(j, j)*Y(j))*Y(j)
    end do
    !    zlike = 0
    !    do  j = 1, nX
    !       ztemp= 0
    !       do  i = 1, nX
    !          ztemp = ztemp + Y(i)*c_inv(i, j)
    !       end do
    !       zlike=zlike+ztemp*Y(j)
    !    end do
    !   zlike = CAMSpec_Quad(c_inv, Y)

    ! compute prior on beams and calibration
    call beams_and_cal_prior(logprior,freq_params)

    !PRINT *,zlike/2,logprior/2,(zlike+logprior)/2

    zlike = zlike + logprior

    end subroutine calc_like

    subroutine Matrix_inverse_internal(M)
      !Inverse of symmetric matrix
      real(8), intent(inout):: M(:,:)
      integer i,j,n
      integer info

      n=Size(M,DIM=1)
      call dpotrf ('L', n, M, n, info)
      if (info/=0)  stop 'Matrix_Cholesky: not positive definite '
      call dpotri ('L', n, M, n, info)
      if (info/=0) stop 'Matrix_inverse: error '
      do i=1,n
          do j=1,i-1
              M(j,i) = M(i,j)
          end do
      end do
    end subroutine Matrix_inverse_internal

end module temp_like
