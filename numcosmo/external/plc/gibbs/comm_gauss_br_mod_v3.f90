module comm_gauss_br_mod_v3
  implicit none

  ! *********************************************************************
  ! *      comm_gauss_br_mod -- An F90 module for computing the         *
  ! *           Gaussianized Blackwell-Rao estimator                    *
  ! *                                                                   *
  ! *          H. K. Eriksen (Oslo) and I. K. Wehus (JPL)               *
  ! *                                                                   *
  ! *  Please cite the following papers when using this code:           *
  ! *                                                                   *
  ! *  - H. K. Eriksen et al. 2008, ApJ, 676, 10  (Commander)           *
  ! *  - M. Chu et al. 2005, Phys. Rev. D, 71, 103002 (Blackwell-Rao)   *
  ! *  - Ø. Rudjord et al. 2009, 2009, ApJ, 692, 1669 (Gauss BR)        *
  ! *                                                                   *
  ! *********************************************************************

  ! ================================================================================================
  ! User routines:
  !
  !  ** subroutine comm_gauss_br_initialize_object(gaussfile, lmin, lmax, max_delta_l, handle)
  !
  !    Input parameters:
  !      gaussfile   (character=*) :: FITS file containing Gaussianized BR data objects
  !      lmin        (i4b)         :: Lowest multipole to include in BR estimator
  !      lmax        (i4b)         :: Highest multipole to include in BR estimator
  !      delta_l     (i4b)         :: Width of banded covariance matrix
  ! 
  !    Output parameters:
  !      handle      (i4b)         :: this code supports multiple BR objects (optional; default=1), and
  !                                   the handle identifies specific BR data sets
  !
  ! 
  !  ** function comm_gauss_br_compute_lnL(cls, handle)
  !
  !    Input parameters:
  !      cls         (dp)          :: array containing at least cls(2:lmax) in units of l(l+1)/2pi
  !      handle      (i4b)         :: handle selecting which BR estimator to use (optional; default = 1)
  !    
  !  
  !  ** subroutine comm_gauss_br_deallocate_object(handle)
  ! 
  !    Input parameters:
  !      handle      (i4b)         :: handle of object to deallocate (optional; default=remove all)
  !
  ! ================================================================================================

  integer, parameter, private :: i4b = selected_int_kind(9)
  integer, parameter, private :: sp  = selected_real_kind(5,30)
  integer, parameter, private :: dp  = selected_real_kind(12,200)
  integer, parameter, private :: lgt = kind(.true.)

  type comm_gauss_br_data
     integer(i4b) :: lmin, lmax, nbin
     real(dp)     :: offset
     real(dp), allocatable, dimension(:)     :: mu, mu_sigma
     real(dp), allocatable, dimension(:,:)   :: cov, prior
     real(dp), allocatable, dimension(:,:,:) :: cl2x
  end type comm_gauss_br_data

  integer(i4b),             parameter,         private :: N_gauss = 100
  type(comm_gauss_br_data), dimension(N_gauss)         :: comm_gauss_br

contains

  ! Initialization routines
  subroutine comm_gauss_br_initialize_object(gaussfile, lmin, lmax, delta_l, handle)
    implicit none
    
    integer(i4b),      intent(in)  :: lmin, lmax, delta_l
    character(len=*),  intent(in)  :: gaussfile
    integer(i4b),      intent(out), optional :: handle
    
    integer(i4b)       :: i, j, l, pos(1), id, info
    real(dp)           :: cl
    
    id = 1
    if (present(handle)) then
       do while (allocated(comm_gauss_br(id)%mu) .and. id < N_gauss)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > N_gauss) stop 'Error -- comm_gauss_br_mod: id out of range'
    
    ! Initialize Gaussianized BR object
    if (allocated(comm_gauss_br(id)%mu))       deallocate(comm_gauss_br(id)%mu)
    if (allocated(comm_gauss_br(id)%mu_sigma)) deallocate(comm_gauss_br(id)%mu_sigma)
    if (allocated(comm_gauss_br(id)%cov))      deallocate(comm_gauss_br(id)%cov)
    if (allocated(comm_gauss_br(id)%cl2x))     deallocate(comm_gauss_br(id)%cl2x)
    call read_gauss_BR_datafile_int(gaussfile, lmin, lmax, comm_gauss_br(id)%mu, &
         & comm_gauss_br(id)%cov, comm_gauss_br(id)%cl2x, comm_gauss_br(id)%mu_sigma)
    comm_gauss_br(id)%lmin        = lmin
    comm_gauss_br(id)%lmax        = lmax
    comm_gauss_br(id)%nbin        = size(comm_gauss_br(id)%cl2x,1)

    ! Bandlimit covariance matrix
    do i = lmin, lmax
       do j = lmin, lmax
          if (abs(i-j) > delta_l) comm_gauss_br(id)%cov(i,j) = 0.d0
       end do
    end do
    
    ! Set up priors
    allocate(comm_gauss_br(id)%prior(lmin:lmax,2))
    do l = lmin, lmax
       j = 1
       do while (abs(comm_gauss_br(id)%cl2x(j,l,2)+5.d0) < 1.d-4)
          j = j+1
       end do
       comm_gauss_br(id)%prior(l,1) = comm_gauss_br(id)%cl2x(j+2,l,1)
    
       j = comm_gauss_br(id)%nbin
       do while (abs(comm_gauss_br(id)%cl2x(j,l,2)-5.d0) < 1.d-4)
          j = j-1
       end do
       comm_gauss_br(id)%prior(l,2) = comm_gauss_br(id)%cl2x(j-2,l,1)
    end do

    ! Invert covariance matrix
    call DPOTRF('L', lmax-lmin+1, comm_gauss_br(id)%cov, lmax-lmin+1, info)
    if (info /= 0) then
       write(*,*) 'DGETRF: Factorization of covariance matrix failed. Info = ', info
       stop
    else
       call DPOTRI('L', lmax-lmin+1, comm_gauss_br(id)%cov, lmax-lmin+1, info)
    end if

    do i = lmin, lmax
       do j = lmin+1, lmax
          comm_gauss_br(id)%cov(i,j) = comm_gauss_br(id)%cov(j,i)
       end do
    end do

    ! Compute offset for "chisq"-like normalization, defined by lnL(mean(sigma_l)) = 0
    comm_gauss_br(id)%offset = 0.d0
    comm_gauss_br(id)%offset = comm_gauss_br_compute_lnL(comm_gauss_br(id)%mu_sigma(2:lmax), handle)
    
  end subroutine comm_gauss_br_initialize_object

  subroutine comm_gauss_br_deallocate_object(handle)
    implicit none
    
    integer(i4b), optional :: handle

    integer(i4b) :: i, id_min, id_max

    id_min = 1; id_max = N_gauss; 
    if (present(handle)) then 
       id_min = handle
       id_max = handle
    end if

    do i = id_min, id_max
       if (allocated(comm_gauss_br(i)%mu))    deallocate(comm_gauss_br(i)%mu)
       if (allocated(comm_gauss_br(i)%cov))   deallocate(comm_gauss_br(i)%cov)
       if (allocated(comm_gauss_br(i)%cl2x))  deallocate(comm_gauss_br(i)%cl2x)
       if (allocated(comm_gauss_br(i)%prior)) deallocate(comm_gauss_br(i)%prior)
       comm_gauss_br(i)%lmin        = -1
       comm_gauss_br(i)%lmax        = -1
    end do

  end subroutine comm_gauss_br_deallocate_object


  ! Base computation routine
  function comm_gauss_br_compute_lnL(cls, handle)
    implicit none

    real(dp),     dimension(2:), intent(in)           :: cls
    real(dp)                                          :: comm_gauss_br_compute_lnL
    integer(i4b),                intent(in), optional :: handle

    integer(i4b) :: i, j, l, id, lmin, lmax, nbin
    real(dp)     :: lnL, dxdCl
    real(dp), allocatable, dimension(:) :: x

    id = 1; if (present(handle)) id = handle    
    if (.not. allocated(comm_gauss_br(id)%mu)) then
       write(*,*) 'Error -- comm_gauss_br_mod: Requested handle ', id, ' is not allocated'
       stop
    end if
    lmin = comm_gauss_br(id)%lmin
    lmax = comm_gauss_br(id)%lmax

    ! Check that spectrum is within the allowed range
    nbin = size(comm_gauss_br(id)%cl2x(:,l,1))
    do l = lmin, lmax
       if (cls(l) < comm_gauss_br(id)%prior(l,1) .or. &
            & cls(l) > comm_gauss_br(id)%prior(l,2)) then
          comm_gauss_br_compute_lnL = -1.d30
          return
       end if
    end do

    ! Convert C_l's to Gaussianized variables
    lmin = comm_gauss_br(id)%lmin
    lmax = comm_gauss_br(id)%lmax
    allocate(x(lmin:lmax))
    do l = lmin, lmax
       x(l) = splint_gauss_br(comm_gauss_br(id)%cl2x(:,l,1), comm_gauss_br(id)%cl2x(:,l,2), &
            & comm_gauss_br(id)%cl2x(:,l,3), cls(l))
    end do

    ! Compute the chi-square term of the Gaussianized Blackwell-Rao estimator
    lnL = -0.5d0 * sum((x-comm_gauss_br(id)%mu) * matmul(comm_gauss_br(id)%cov, &
         & (x-comm_gauss_br(id)%mu)))
    
    ! Add the Jacobian terms from the transformation
    do l = lmin, lmax
       dxdCl = splint_deriv_gauss_br(comm_gauss_br(id)%cl2x(:,l,1), &
            & comm_gauss_br(id)%cl2x(:,l,2), comm_gauss_br(id)%cl2x(:,l,3), cls(l))
       if (dxdCl <= 0.d0) then
          lnL = -1.d30
          exit
       else
          lnL = lnL + log(dxdCl)
       end if
    end do

    comm_gauss_br_compute_lnL = lnL - comm_gauss_br(id)%offset

  end function comm_gauss_br_compute_lnL

  function comm_gauss_br_getlun()
    implicit none
    integer(i4b) :: comm_gauss_br_getlun
    logical(lgt) :: exists, isopen
    comm_gauss_br_getlun = 9
    do
       comm_gauss_br_getlun = comm_gauss_br_getlun+1
       inquire(unit=comm_gauss_br_getlun,exist=exists)
       if(exists) then
          inquire(unit=comm_gauss_br_getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function comm_gauss_br_getlun
  

  subroutine read_gauss_BR_datafile_int(filename, lmin, lmax, mu, cov, cl2x, mu_sigma)
    implicit none

    character(len=*),                        intent(in)  :: filename
    integer(i4b),                            intent(in)  :: lmin, lmax
    real(dp), allocatable, dimension(:),     intent(out) :: mu, mu_sigma
    real(dp), allocatable, dimension(:,:),   intent(out) :: cov
    real(dp), allocatable, dimension(:,:,:), intent(out) :: cl2x
    
    integer(i4b)         :: i, j, l, m, unit
    integer(i4b)         :: status, blocksize, readwrite, hdutype, nbin
    integer(i4b)         :: nelements, fpixel, group, bitpix, lmax_in, lmin_in
    logical(lgt)         :: simple, extend, anyf, exist
    real(dp)             :: nullval
    character(len=80)    :: comment, errorline
    real(dp), allocatable, dimension(:)     :: mu_in, mu_sigma_in
    real(dp), allocatable, dimension(:,:)   :: cov_in
    real(dp), allocatable, dimension(:,:,:) :: cl2x_in

    unit  = comm_gauss_br_getlun()
    m     = len(trim(filename))

    status    = 0
    readwrite = 0
    nullval   = 0.d0
    group     = 1
    fpixel    = 1

    ! Open the data file
    call ftopen(unit, trim(filename), readwrite, blocksize, status)

    ! Read data vector and keywords
    call ftmahd(unit, 1, hdutype, status)
    call ftgkyj(unit, 'LMIN',  lmin_in, comment,status)
    call ftgkyj(unit, 'LMAX',  lmax_in, comment,status)
    call ftgkyj(unit, 'NBIN',  nbin,    comment,status)
    allocate(cl2x_in(nbin,lmin_in:lmax_in,3), mu_in(lmin_in:lmax_in), mu_sigma_in(lmin_in:lmax_in))
    allocate(cov_in(lmin_in:lmax_in,lmin_in:lmax_in))
    call ftgpvd(unit, group, fpixel, size(cl2x_in), nullval, cl2x_in, anyf, status)

    ! Read mean vector
    call ftmahd(unit, 2, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(mu_in), nullval, mu_in, anyf, status)

    ! Read covariance matrix
    call ftmahd(unit, 3, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(cov_in), nullval, cov_in, anyf, status)

    ! Read mean sigma vector
    call ftmahd(unit, 4, hdutype, status)
    call ftgpvd(unit, group, fpixel, size(mu_sigma_in), nullval, mu_sigma_in, anyf, status)

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

    if (lmax > lmax_in) then
       write(*,*) 'Error: Requested lmax greater than that supported by ', trim(filename)
       stop
    end if

    if (lmin < lmin_in) then
       write(*,*) 'Error: Requested lmin smaller than that supported by ', trim(filename)
       stop
    end if

    ! Copy relevant data into output data structures
    allocate(cl2x(nbin,lmin:lmax,3), mu(lmin:lmax), cov(lmin:lmax,lmin:lmax), mu_sigma(2:lmax))
    cl2x     = cl2x_in(:,lmin:lmax,:)
    mu       = mu_in(2:lmax)
    cov      = cov_in(lmin:lmax,lmin:lmax)
    mu_sigma = mu_sigma_in(2:lmax)

    deallocate(cl2x_in, mu_in, cov_in)
    
  end subroutine read_gauss_BR_datafile_int


  function splint_gauss_br(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_gauss_br

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n   = size(xa)
    klo = max(min(locate_gauss_br(xa,x),n-1),1)
    khi = klo+1
    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    splint_gauss_br = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo) + (b**3-b)*y2a(khi))*(h**2)/6.d0
    
  end function splint_gauss_br

  function splint_deriv_gauss_br(xa, ya, y2a, x)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: xa, ya, y2a
    real(dp)                            :: splint_deriv_gauss_br

    integer(i4b) :: khi, klo, n
    real(dp)     :: a, b, h

    n   = size(xa)
    klo = max(min(locate_gauss_br(xa,x),n-1),1)
    khi = klo+1
    h   = xa(khi) - xa(klo)
    a   = (xa(khi) - x) / h
    b   = (x - xa(klo)) / h
    splint_deriv_gauss_br = (ya(khi) - ya(klo)) / h - (3.d0 * a**2 - 1.d0) / 6.d0 * h * y2a(klo) + &
         & (3.d0 * b**2 - 1.d0) / 6.d0 * h * y2a(khi)

  end function splint_deriv_gauss_br

  function locate_gauss_br(xx, x)
    implicit none

    real(dp),               intent(in) :: x
    real(dp), dimension(:), intent(in) :: xx
    integer(i4b)                       :: locate_gauss_br

    integer(i4b) :: n, jl, jm, ju
    logical(lgt) :: ascnd

    n     = size(xx)
    ascnd = (xx(n) >= xx(1))
    jl    = 0
    ju    = n+1

    do 
       if (ju-jl <= 1) exit
       jm = (ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl = jm
       else
          ju = jm
       end if
    end do

    if (x == xx(1)) then
       locate_gauss_br = 1
    else if (x == xx(n)) then
       locate_gauss_br = n-1
    else
       locate_gauss_br = jl
    end if

  end function locate_gauss_br


end module comm_gauss_br_mod_v3
