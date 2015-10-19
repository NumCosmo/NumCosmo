module comm_lowl_mod_dist
  implicit none

  ! *********************************************************************
  ! *  comm_lowl_mod -- An F90 module for computing a Gaussian          *
  ! *                   low-l likelihood by brute force                 *
  ! *                                                                   *
  ! *           H. K. Eriksen, E. Gjerløw (University of Oslo)          *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *      L. P. L. Colombo, J. B. Jewell and I. K. Wehus (JPL)         *   
  ! *                                                                   *
  ! *                                                                   *
  ! *  Please cite the following papers when using this code:           *
  ! *                                                                   *
  ! *  - Gjerløw et al. 2014, ApJ, in preparation                       *
  ! *                                                                   *
  ! *  History:                                                         *
  ! *      April 25th, 2014 -- First fully functional version           *
  ! *                                                                   *
  ! *********************************************************************

  ! ==========================================================================
  ! User routines:
  !
  !  ** subroutine comm_lowl_initialize_object(paramfile, handle)
  !
  !    Input parameters:
  !      paramfile   (character=*) :: Low-l parameter file
  ! 
  !    Output parameters:
  !      handle (i4b)  :: this code supports multiple objects (optional; default=1),
  !                       and the handle identifies specific data sets
  !
  ! 
  !  ** function comm_lowl_compute_lnL(cls, ierr, handle)
  !
  !    Input parameters:
  !      cls     (dp)  :: array containing at least cls(2:lmax) in units of l(l+1)/2pi
  !      ierr    (i4b) :: error flag
  !      handle  (i4b) :: handle selecting which data set to use (optional; default = 1)
  !    
  !  
  !  ** subroutine comm_lowl_deallocate_object(handle)
  ! 
  !    Input parameters:
  !      handle  (i4b) :: handle of object to deallocate (optional; default=remove all)
  !
  ! ==============================================================================

  integer,      parameter, private :: i4b = selected_int_kind(9)
  integer,      parameter, private :: sp  = selected_real_kind(5,30)
  integer,      parameter, private :: dp  = selected_real_kind(12,200)
  integer,      parameter, private :: lgt = kind(.true.)
  real(dp),     parameter, private :: pi = 3.141592653589793238462643383279502d0

  integer(i4b), parameter, private :: MAX_N_LOWL = 5

  type comm_lowl_data
     logical(lgt) :: initialized=.false.
     integer(i4b) :: n, n_h, n_d, nmaps, llow, lhigh, lmax
     real(dp)     :: loglike_weight, cond_threshold
     real(dp)     :: lnL_recent, chisq_recent, red_chisq_recent
     real(dp)     :: cond_number_recent
     real(dp), allocatable, dimension(:,:)   :: d
     real(dp), allocatable, dimension(:,:)   :: cl_fid
     real(dp), allocatable, dimension(:,:)   :: N_cov, P_harm
     real(dp), allocatable, dimension(:,:,:) :: w
     real(dp), allocatable, dimension(:,:)   :: beam
  end type comm_lowl_data

  type(comm_lowl_data), dimension(MAX_N_LOWL), private :: comm_lowl
  
  private

  public :: comm_lowl_initialize_object
  public :: comm_lowl_deallocate_object
  public :: comm_lowl_compute_lnL
  public :: comm_lowl_get_recent_value
  public :: read_fiducial_spectrum

contains

  ! Initialization routines
  subroutine comm_lowl_initialize_object(paramfile, handle)
    implicit none

    character(len=*),  intent(in)  :: paramfile
    integer(i4b),      intent(out), optional :: handle

    integer(i4b)       :: firstchain, lastchain, firstsample, lastsample, thinstep
    integer(i4b)       :: b, i, j, k, l, m, n, q, c1, c2, f, ind, d, n_h, n_d
    integer(i4b)       :: unit, numsamples, numchains, lmax_chain, id, bin(2), n_p, n_g, nmode
    integer(i4b)       :: col, p, nmaps, nspec, nsamp, lmax_cl, lmin_bin, lmax_bin, ncomp
    logical(lgt)       :: pattern(3,3), polarization, exist
    real(dp)           :: lnL, t1, t2
    character(len=512) :: line, s1, s2, datafile
    character(len=128) :: sigmafile, clfile
    real(dp),     allocatable, dimension(:,:)     :: cls
    integer(i4b), allocatable, dimension(:)       :: i2p
    real(dp),     allocatable, dimension(:,:,:,:) :: sigma
    real(dp),     allocatable, dimension(:,:,:)   :: sigma_1D
    real(dp),     allocatable, dimension(:,:)     :: S, sqrtS_P

    id = 1
    unit = comm_getlun()
    if (present(handle)) then
       do while (comm_lowl(id)%initialized .and. id < MAX_N_LOWL)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > MAX_N_LOWL) stop 'Error -- planck_br_mod: lowl id out of range'

    ! Initialize all distributions
    comm_lowl(id)%initialized = .true.
    call comm_get_parameter(paramfile, 'DATAFILE',         par_string=datafile)
    call comm_get_parameter(paramfile, 'FIDUCIAL_CL_FILE', par_string=clfile)
    call comm_get_parameter(paramfile, 'LOGLIKE_WEIGHT',   par_dp=comm_lowl(id)%loglike_weight)
    call comm_get_parameter(paramfile, 'LMIN',             par_int=comm_lowl(id)%llow)
    call comm_get_parameter(paramfile, 'LMAX',             par_int=comm_lowl(id)%lhigh)
    call comm_get_parameter(paramfile, 'CONDITION_NUMBER_THRESHOLD', &
         & par_dp=comm_lowl(id)%cond_threshold)

    ! Read data file
    inquire(file=trim(datafile), exist=exist)
    if (.not. exist) then
       write(*,*) 'Error -- low-l datafile = ', trim(datafile), ' does not exist'
       stop
    end if

    call read_lowl_datafile(datafile, comm_lowl(id)%n, comm_lowl(id)%n_d, comm_lowl(id)%n_h, &
         & comm_lowl(id)%lmax, comm_lowl(id)%nmaps, comm_lowl(id)%d, comm_lowl(id)%N_cov, &
         & comm_lowl(id)%beam, comm_lowl(id)%P_harm, comm_lowl(id)%w)
    ncomp               = (comm_lowl(id)%lmax+1)**2
    nmode = comm_lowl(id)%n
    nmaps = comm_lowl(id)%nmaps
    n_h   = comm_lowl(id)%n_h

    if (comm_lowl(id)%lhigh > comm_lowl(id)%lmax) then
       write(*,*) 'comm_lowl_mod -- Error: Requested LMAX greater than maximum '
       write(*,*) '                        multipole in datafile'
       stop
    end if

    ! Read in fiducial spectrum, to be conditioned upon outside range of interest
    call read_fiducial_spectrum(clfile, comm_lowl(id)%cl_fid)

    ! Add basis vectors with fixed C_l's directly into noise covariance
    allocate(sqrtS_P(n_h,nmode), S(nmaps,nmaps))
    sqrtS_P = 0.d0
    ind     = 5
    do l = 2, comm_lowl(id)%lmax
       if (l >= comm_lowl(id)%llow .and. l <= comm_lowl(id)%lhigh) then
          ind = ind + 2*l+1
          cycle
       end if
       call cl2s(comm_lowl(id)%cl_fid(l,:), S) ! Fix to fiducial
       S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
       call cholesky_decompose_with_mask_dp(S)
       do m = -l, l
          call dgemm('N', 'N', nmaps, nmode, nmaps, 1.d0, S, nmaps, &
               & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
               & nmaps, 0.d0, sqrtS_P(ind:n_h:ncomp,:), nmaps)
          ind = ind+1
       end do
    end do
    call dsyrk('L','T',nmode,n_h,1.d0,sqrtS_P,n_h,1.d0,comm_lowl(id)%N_cov,nmode)
    deallocate(sqrtS_P, S)

    ! Symmetrize matrix
    do i = 1, nmode
       do j = i+1, nmode
          comm_lowl(id)%N_cov(i,j) = comm_lowl(id)%N_cov(j,i)
       end do
    end do

  end subroutine comm_lowl_initialize_object

  subroutine comm_lowl_deallocate_object(handle)
    implicit none
    
    integer(i4b), optional :: handle

    integer(i4b) :: i, j, k, id, id_min, id_max

    id_min = 1; id_max = MAX_N_LOWL; 
    if (present(handle)) then 
       id_min = handle
       id_max = handle
    end if

    do id = id_min, id_max
       if (comm_lowl(id)%initialized) then
          deallocate(comm_lowl(id)%d, comm_lowl(id)%cl_fid, comm_lowl(id)%N_cov)
          deallocate(comm_lowl(id)%P_harm, comm_lowl(id)%w, comm_lowl(id)%beam)
       end if
      comm_lowl(id)%initialized = .false.
    end do

  end subroutine comm_lowl_deallocate_object


  ! Base computation routine
  function comm_lowl_compute_lnL(cls, sqrt_S, chisq, red_chisq, handle, enforce_pos_def, &
       & lnL_multi, cond_number, ierr)
    implicit none

    real(dp),     dimension(0:,1:), intent(in),  optional :: cls
    integer(i4b),                   intent(out), optional :: ierr
    integer(i4b),                   intent(in),  optional :: handle
    logical(lgt),                   intent(in),  optional :: enforce_pos_def
    real(dp),     dimension(1:,1:), intent(in),  optional :: sqrt_S
    real(dp),                       intent(out), optional :: chisq, red_chisq, cond_number
    real(dp),     dimension(1:),    intent(out), optional :: lnL_multi
    real(dp)                                              :: comm_lowl_compute_lnL

    integer(i4b) :: i, j, l, m, n, id, k, stat, n_d
    logical(lgt) :: posdef, enf_pos_def
    real(dp)     :: chi2, logdet, t1, t2, L_max, L_min, cond
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: C, invC_d, V, map

    if (present(ierr)) ierr = 0
    id = 1; if (present(handle)) id = handle
    enf_pos_def = .true.; if (present(enforce_pos_def)) enf_pos_def = enforce_pos_def

    ! Check that likelihood structure is initialized
    if (.not. comm_lowl(id)%initialized) then
       write(*,*) 'Error -- comm_lowl_mod: Requested handle ', id, ' is not initialized'
       stop
    end if

    ! Compute covariance matrix
    n   = comm_lowl(id)%n
    n_d = comm_lowl(id)%n_d
    allocate(C(n,n))
    if (present(sqrt_S)) then
       call comm_lowl_getC(enf_pos_def, C, stat, sqrt_S=sqrt_S)
    else
       call comm_lowl_getC(enf_pos_def, C, stat, cls=cls)
    end if

    if (stat /= 0) then
       if (present(ierr)) ierr = 1
       comm_lowl_compute_lnL = -1.d30
       deallocate(C)
       return
    end if

    ! Cholesky decompose matrix
    call dpotrf('L', n, C, n, stat)
    if (stat /= 0) then
       if (present(ierr)) ierr = ierr + 1
       comm_lowl_compute_lnL = -1.d30
       deallocate(C)
       return
    end if

    ! Compute log-determinant
    logdet =  0.d0
    L_max  = -1.d30
    L_min  =  1.d30
    do i = 1, n
       logdet = logdet + 2.d0 * log(C(i,i))
       L_max  = max(L_max, C(i,i))
       L_min  = min(L_min, C(i,i))
    end do
    cond = (L_max/L_min)**2
    if (present(cond_number)) cond_number = cond

    ! Compute chi-square term
    allocate(invC_d(n,n_d))
    invC_d = comm_lowl(id)%d
    call dpotrs('L', n, n_d, C, n, invC_d, n, stat)

    ! Return log-like value
    if (stat == 0 .and. cond < comm_lowl(id)%cond_threshold) then
       do i = 1, n_d
          chi2 = sum(comm_lowl(id)%d(:,i)*invC_d(:,i))
          if (present(lnL_multi)) lnL_multi(i) = -0.5d0 * (chi2 + logdet)
          if (i == 1) then
             if (present(chisq))     chisq     = chi2
             if (present(red_chisq)) red_chisq = chi2 / n
             comm_lowl_compute_lnL = -0.5d0 * (chi2 + logdet)
          end if
       end do
    else
       comm_lowl_compute_lnL = -1.d30
       if (present(ierr))      ierr      = 1
       if (present(lnL_multi)) lnL_multi = -1.d30
    end if

    ! Update data structures for quick look-up 
    comm_lowl(id)%lnL_recent         = comm_lowl_compute_lnL
    comm_lowl(id)%chisq_recent       = sum(comm_lowl(id)%d(:,1)*invC_d(:,1))
    comm_lowl(id)%red_chisq_recent   = comm_lowl(id)%chisq_recent / n
    comm_lowl(id)%cond_number_recent = cond

    deallocate(C, invC_d)

  end function comm_lowl_compute_lnL

  subroutine comm_lowl_get_recent_value(handle, cond_number, chisq, red_chisq, lnL)
    implicit none

    integer(i4b), intent(in),  optional :: handle
    real(dp),     intent(out), optional :: cond_number, chisq, red_chisq, lnL

    integer(i4b) :: id

    id = 1; if (present(handle)) id = handle
    if (present(cond_number)) cond_number = comm_lowl(id)%cond_number_recent
    if (present(chisq))       chisq       = comm_lowl(id)%chisq_recent
    if (present(red_chisq))   red_chisq   = comm_lowl(id)%red_chisq_recent
    if (present(lnL))         lnL         = comm_lowl(id)%lnL_recent

  end subroutine comm_lowl_get_recent_value


  subroutine comm_lowl_getC(enforce_pos_def, C, ierr, handle, cls, sqrt_S)
    implicit none

    real(dp),              dimension(0:,1:), intent(in), optional :: cls
    real(dp),              dimension(1:,1:), intent(in), optional :: sqrt_S
    integer(i4b),                            intent(in), optional :: handle
    logical(lgt),                            intent(in)           :: enforce_pos_def
    real(dp), allocatable, dimension(:,:),   intent(out)          :: C
    integer(i4b),                            intent(out)          :: ierr

    real(dp)     :: t1, t2, var
    integer(i4b) :: i, j, l, m, ind, ind1, ind2, id, stat
    integer(i4b) :: n_h, n, nmaps, lmax, ncomp, ind_min, ind_max, nmode
    real(dp), allocatable, dimension(:,:) :: S, sqrtS_P, V, P_sub, map

    id = 1; if (present(handle)) id = handle
    n_h     = comm_lowl(id)%n_h
    n       = comm_lowl(id)%n
    nmaps   = comm_lowl(id)%nmaps
    lmax    = comm_lowl(id)%lmax
    ncomp   = comm_lowl(id)%n_h / comm_lowl(id)%nmaps
    ind_min = comm_lowl(id)%llow**2+1
    ind_max = (comm_lowl(id)%lhigh+1)**2
    nmode   = (ind_max-ind_min+1)*nmaps
    ierr    = 0

    allocate(C(n,n), sqrtS_P(nmode,n), S(nmaps,nmaps))

    if (present(sqrt_S)) then
       ! Use user-supplied harmonic space covariance

       ! Extract vectors to be multiplied with sqrt_S, and multiply with beam
       allocate(P_sub(nmode,n))
       ind2 = 1
       do i = 1, nmaps
          ind1 = (i-1)*ncomp + comm_lowl(id)%llow**2 + 1
          do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
             do m = -l, l
                P_sub(ind2,:) = comm_lowl(id)%beam(l,i)*comm_lowl(id)%P_harm(ind1,:)
                ind2 = ind2+1
                ind1 = ind1+1
             end do
          end do
       end do

       ! Multiply basis vectors with sqrt(S), ie., Cholesky factor of S
       call dgemm('T', 'N', nmode, n, nmode, 1.d0, sqrt_S, nmode, &
            & P_sub, nmode, 0.d0, sqrtS_P, nmode)
       call dsyrk('L','T',n,nmode,1.d0,sqrtS_P,nmode,0.d0,C,n)
       do i = 1, n
          do j = i+1, n
             C(i,j) = C(j,i)
          end do
       end do
       deallocate(P_sub)

    else if (enforce_pos_def) then
       ! Use user-supplied power spectrum for S

       ! Check that spectrum is positive definite
       if (.not. comm_cls_posdef(cls(2:lmax,:))) then
          ierr = 1
          C    = 0.d0
          deallocate(sqrtS_P, S)
          return
       end if

       ! Compute signal covariance matrix
       sqrtS_P        = 0.d0
       ind            = comm_lowl(id)%llow**2+1
       ind2           = 1
       do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
          call cl2s(cls(l,:), S)      
          S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
          call cholesky_decompose_with_mask_dp(S,ierr=stat)
          if (stat /= 0) then
             ierr = 1       
             C    = 0.d0
             deallocate(sqrtS_P, S)
             return
          end if
          do m = -l, l
             call dgemm('N', 'N', nmaps, n, nmaps, 1.d0, S, nmaps, &
                  & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
                  & nmaps, 0.d0, sqrtS_P(ind2:ind2+nmaps-1,:), nmaps)
             ind  = ind  + 1
             ind2 = ind2 + nmaps
          end do
       end do
       call dsyrk('L','T',n,nmode,1.d0,sqrtS_P,nmode,0.d0,C,n)
       do i = 1, n
          do j = i+1, n
             C(i,j) = C(j,i)
          end do
       end do
       
    else 
       ! Use user-supplied power spectrum for S; do not enforce positive definiteness

       ! Compute signal covariance matrix
       allocate(P_sub(nmode,n))
       sqrtS_P        = 0.d0
       ind            = comm_lowl(id)%llow**2+1
       ind2           = 1
       do l = comm_lowl(id)%llow, comm_lowl(id)%lhigh
          call cl2s(cls(l,:), S)      
          S = S * comm_lowl(id)%w(l,:,:) / (l*(l+1)/(2.d0*pi)) ! w = (b_l*p_l)^2
          do m = -l, l
             call dgemm('N', 'N', nmaps, n, nmaps, 1.d0, S, nmaps, &
                  & comm_lowl(id)%P_harm(ind:n_h:ncomp,:), &
                  & nmaps, 0.d0, sqrtS_P(ind2:ind2+nmaps-1,:), nmaps)
             P_sub(ind2:ind2+nmaps-1,:) = comm_lowl(id)%P_harm(ind:n_h:ncomp,:)
             ind  = ind  + 1
             ind2 = ind2 + nmaps
          end do
       end do
       call dgemm('T', 'N', n, n, nmode, 1.d0, P_sub, nmode, &
            & sqrtS_P, nmode, 0.d0, C, n)

       deallocate(P_sub)
    end if

    ! Add noise
    C = C + comm_lowl(id)%N_cov

    ! Validate against correlation function; only works for pixel space basis
    if (.false.) then
       var = 0.d0
       do l = 2, lmax
          var = var + (2*l+1)/(4.d0*pi) * cls(l,1) * &
               & 2.d0*pi/real(l*(l+1),dp) * comm_lowl(id)%beam(l,1)**2
       end do
       write(*,*) C(1,1), var
       stop
    end if

    deallocate(sqrtS_P, S)

  end subroutine comm_lowl_getC


  ! ==========================================================================
  ! Utility routines
  ! ==========================================================================
  ! *****************************************************************************************
  !
  ! Routine for reading one parameter from a parameter file
  !     Example usage:   call comm_get_parameter(21, "parfile.txt", "NSIDE", par_int=nside)
  !
  ! *****************************************************************************************

  subroutine comm_get_parameter(parfile, parname, par_int, par_char, &
       & par_string, par_sp, par_dp, par_lgt)
    implicit none

    character(len=*),  intent(in)  :: parfile, parname
    integer(i4b),      intent(out), optional :: par_int
    logical(lgt),      intent(out), optional :: par_lgt
    character(len=1),  intent(out), optional :: par_char
    character(len=*),  intent(out), optional :: par_string
    real(sp),          intent(out), optional :: par_sp
    real(dp),          intent(out), optional :: par_dp


    integer(i4b)        :: i, unit
    character(len=128)  :: string, variable, value
    character(len=1)    :: equals

    unit = comm_getlun()
    open(unit, file=trim(parfile))

    do while (.true.)
       read(unit,*,end=1) string

       if (string(1:1)=='#') cycle

       backspace(unit)
       read(unit,*) variable, equals, value

       if (trim(variable) == trim(parname)) then

          if (present(par_int)) then
             read(value,*) par_int
          else if (present(par_char)) then
             read(value,*) par_char
          else if (present(par_string)) then
             read(value,'(a)') par_string
          else if (present(par_sp)) then
             read(value,*) par_sp
          else if (present(par_dp)) then
             read(value,*) par_dp
          else if (present(par_lgt)) then
             read(value,*) par_lgt
          end if

          close(unit)
          return

       end if

    end do

1   write(*,*) 'GET_PARAMETER:    Critical error -- parameter not found'
    write(*,*) 'GET_PARAMETER:       Parameter file = ', trim(parfile)
    write(*,*) 'GET_PARAMETER:       Parameter name = ', trim(parname)

    close(unit)
    stop

  end subroutine comm_get_parameter

  function comm_getlun()
    implicit none
    integer(i4b) :: comm_getlun
    logical(lgt) :: exists, isopen
    comm_getlun = 9
    do
       comm_getlun = comm_getlun+1
       inquire(unit=comm_getlun,exist=exists)
       if(exists) then
          inquire(unit=comm_getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function comm_getlun


  subroutine read_lowl_datafile(filename, n, n_d, n_h, lmax, nmaps, d, C, beam, P_harm, window)
    implicit none

    character(len=*),                        intent(in)  :: filename
    integer(i4b),                            intent(out) :: n, n_d, n_h, lmax, nmaps
    real(dp), allocatable, dimension(:,:),   intent(out) :: d, C, P_harm
    real(dp), allocatable, dimension(:,:),   intent(out) :: beam
    real(dp), allocatable, dimension(:,:,:), intent(out) :: window
    
    integer(i4b)         :: i, j, l, m, unit
    integer(i4b)         :: status, blocksize, readwrite, hdutype
    integer(i4b)         :: nelements, fpixel, group, bitpix
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline

    unit  = comm_getlun()
    m     = len(trim(filename))

    if (filename(m-2:m) == 'unf') then

       open(unit,file=trim(filename),form='unformatted')
       read(unit) lmax, nmaps
       read(unit) n_d, n_h, n
       allocate(window(0:lmax,nmaps,nmaps), beam(0:lmax,nmaps))
       allocate(d(n,n_d), C(n,n), P_harm(n_h,n))
       read(unit) d
       read(unit) C
       read(unit) beam
       read(unit) P_harm
       close(unit)

    else if (filename(m-3:m) == 'fits') then

       status    = 0
       readwrite = 0
       nullval   = 0.d0
       group     = 1
       fpixel    = 1

       ! Open the data file
       call ftopen(unit, trim(filename), readwrite, blocksize, status)

       ! Read data vector and keywords
       call ftmahd(unit, 1, hdutype, status)
       call ftgkyj(unit, 'LMAX',  lmax, comment,status)
       call ftgkyj(unit, 'NMAPS', nmaps, comment,status)
       call ftgkyj(unit, 'N_D',   n_d,  comment,status)
       call ftgkyj(unit, 'N',     n,    comment,status)
       call ftgkyj(unit, 'N_H',   n_h,  comment,status)
       allocate(window(0:lmax,nmaps,nmaps), beam(0:lmax,nmaps))
       allocate(d(n,n_d), C(n,n), P_harm(n_h,n))
       call ftgpvd(unit, group, fpixel, size(d), nullval, d, anyf, status)

       ! Read covariance matrix
       call ftmahd(unit, 2, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(C), nullval, C, anyf, status)

       ! Read beam
       call ftmahd(unit, 3, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(beam), nullval, beam, anyf, status)

       ! Read projection operator
       call ftmahd(unit, 4, hdutype, status)
       call ftgpvd(unit, group, fpixel, size(P_harm), nullval, P_harm, anyf, status)

       ! Close file
       call ftclos(unit, status)
       
       if (status /= 0) then
          write(*,*) 'Error in FITS operations -- status = ', status
          call ftgmsg(errorline)
          do while (trim(errorline) /= '')
             write(*,*) errorline
             call ftgmsg(errorline)
          end do
       end if

    else
       write(*,*) 'Error -- unknown filetype = ', trim(filename)
       stop
    end if

    do l = 0, lmax
       do i = 1, nmaps
          do j = 1, nmaps
             window(l,i,j) = beam(l,i)*beam(l,j)
          end do
       end do
    end do
    
  end subroutine read_lowl_datafile

  subroutine read_fiducial_spectrum(clfile, cls)
    implicit none

    character(len=*),                              intent(in)  :: clfile
    real(dp),         allocatable, dimension(:,:), intent(out) :: cls

    integer(i4b)        :: unit, l, lmax
    real(dp)            :: cl_in(4)
    character(len=2048) :: line

    ! Read fiducial spectrum file
    unit = comm_getlun()
    open(unit,file=trim(clfile))
    lmax = -1
    do while (.true.)
       read(unit,'(a)',end=89) line
       line = trim(line)
       if (line(1:1) == '#') cycle
       read(line,*) l
       lmax = max(lmax, l)
    end do
89  close(unit)
    if (lmax > -1) then
       allocate(cls(0:lmax,6))
       cls = 0.d0
       open(unit,file=trim(clfile))
       do while (.true.)
          read(unit,'(a)',end=90) line
          line = trim(line)
          if (line(1:1) == '#') cycle
          read(line,*) l, cl_in
          cls(l,1) = cl_in(1) ! Assume (TT, EE, BB, TE) ordering for now
          cls(l,2) = cl_in(4)          
          cls(l,4) = cl_in(2)
          cls(l,6) = cl_in(3)
       end do
90     close(unit)
    else
       write(*,*) 'Error -- comm_mod: No valid entries in clfile = ', trim(clfile)
       stop
    end if

  end subroutine read_fiducial_spectrum

  subroutine cholesky_decompose_with_mask_dp(M, ierr)
    implicit none

    real(dp),     dimension(1:,1:), intent(inout) :: M
    integer(i4b),                   intent(out), optional :: ierr

    integer(i4b)     :: i, j, q, n, info
    integer(i4b), allocatable, dimension(:)   :: indmap
    real(dp),     allocatable, dimension(:,:) :: L

    if (present(ierr)) ierr = 0

    q     = size(M,1)
    allocate(indmap(q))
    n = 0
    do i = 1, q
       if (M(i,i) > 0.d0) then
          n = n+1
          indmap(n) = i
       else if (M(i,i) < 0.d0) then
          M    = 0.d0
          ierr = 1
          deallocate(indmap)
          return
       end if
    end do
    if (n == 0) then
       M = 0.d0
       deallocate(indmap)
       return
    end if

    allocate(L(n,n))
    L = M(indmap(1:n),indmap(1:n))

    call dpotrf('U', n, L, n, info)
    if (info /= 0) then
       if (present(ierr)) then
          ierr = 1
       else
          write(*,*) 'DPOTRF: Cholesky factorization failed. Info = ', info
          stop
       end if
    end if

    ! Nullify lower triangle and missing rows/columns
    do i = 1, n
       L(i,1:i-1) = 0.d0
    end do

    ! Copy back to main matrix
    M                          = 0.d0
    M(indmap(1:n),indmap(1:n)) = L

    deallocate(indmap, L)

  end subroutine cholesky_decompose_with_mask_dp

  function comm_cls_posdef(cls)
    implicit none

    real(dp),     dimension(:,:), intent(in) :: cls
    logical(lgt)                             :: comm_cls_posdef

    integer(i4b) :: i, nmaps
    logical(lgt) :: ok
    real(dp), allocatable, dimension(:)   :: W
    real(dp), allocatable, dimension(:,:) :: S

    nmaps = size(cls,2)
    comm_cls_posdef = .true.
    do i = 1, size(cls,1)
       if (nmaps == 1) then
          ! Only temperature
          ok = cls(i,1) >= 0.d0
       else if (cls(i,3) == 0.d0 .and. cls(i,5) == 0.d0) then
          ! TB = EB = 0 => 2x2 + 1
          ok = cls(i,1)*cls(i,4)-cls(i,2)**2 >= 0.d0 .and. cls(i,6) >= 0.d0
       else
          ! All six spectra are included
          allocate(S(nmaps,nmaps), W(nmaps))
          call cl2S(cls(i,:), S)
          call comm_get_eigenvalues(S, W)
          ok = all(W > 0.d0)
          deallocate(S, W)
       end if
       if (.not. ok) then
          comm_cls_posdef = .false.
          return
       end if
    end do

  end function comm_cls_posdef

  subroutine cl2s(cl, S)
    implicit none

    real(dp), dimension(1:),    intent(in)  :: cl
    real(dp), dimension(1:,1:), intent(out) :: S

    integer(i4b) :: i, j, k, nmaps

    nmaps = size(S,1)

    S = 0.d0
    k = 1
    do i = 1, nmaps
       do j = i, nmaps
          S(i,j) = cl(k)
          S(j,i) = cl(k)
          k      = k+1
       end do
    end do

  end subroutine cl2s

  subroutine comm_get_eigenvalues(A, eigenvals)
    
    real(dp), dimension(1:,1:), intent(in)            :: A
    real(dp), dimension(1:),    intent(out)           :: eigenvals

    integer(i4b)     :: n, liwork, lwork, lda, info
    character(len=1) :: job, uplo
    real(dp),     allocatable, dimension(:,:) :: A_copy
    real(dp),     allocatable, dimension(:)   :: W, work
    integer(i4b), allocatable, dimension(:)   :: iwork    

    n      = size(eigenvals)

    if (n == 1) then

       eigenvals(1) = A(1,1)

    else if (n == 2) then

       eigenvals(1) = 0.5d0*(A(1,1)+A(2,2) + sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))
       eigenvals(2) = 0.5d0*(A(1,1)+A(2,2) - sqrt(4.d0*A(1,2)*A(2,1)+(A(1,1)-A(2,2))**2))

    else

       job    = 'n'
       uplo   = 'l'
       lda    = n
       liwork = 1
       lwork  = 2*n+1

       ! Perform eigenvalue decomposition
       allocate(work(lwork))
       allocate(iwork(liwork))
       allocate(A_copy(n,n))
       A_copy = A
       call dsyevd(job, uplo, n, A_copy, lda, eigenvals, work, lwork, iwork, liwork, info)
       if (info /= 0) write(*,*) 'get_eigenvalues -- dsyevd info = ', info
       
       deallocate(work)
       deallocate(iwork)
       deallocate(A_copy)

    end if

  end subroutine comm_get_eigenvalues

end module comm_lowl_mod_dist
