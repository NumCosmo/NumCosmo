module comm_br_mod
  implicit none

  ! *********************************************************************
  ! *  comm_br_mod -- An F90 module for computing the Blackwell-Rao     *
  ! *                 estimator given signal samples from the posterior *
  ! *                                                                   *
  ! *  H. K. Eriksen, E. Gjerløw, K. Mikkelsen (University of Oslo)     *
  ! *                                                                   *   
  ! *                               and                                 *
  ! *                                                                   *   
  ! *  I. K. Wehus, D. Pietrobon, K. M. Gorski, L. Colombo, G. Rocha,   *
  ! *  J. Bartlett, C. R. Lawrence, J. B. Jewell, S. Hildebrandt (JPL)  *   
  ! *                                                                   *
  ! *                                                                   *
  ! *  Please cite the following papers when using this code:           *
  ! *                                                                   *
  ! *  - H. K. Eriksen et al. 2008, ApJ, 676, 10  (Commander)           *
  ! *  - M. Chu et al. 2005, Phys. Rev. D, 71, 103002 (Blackwell-Rao)   *
  ! *                                                                   *
  ! *********************************************************************

  ! ================================================================================================
  ! User routines:
  !
  !  ** subroutine comm_br_initialize_object(sigmafile, clfile, lmin, lmax, firstchain, lastchain, &
  !           & firstsample, lastsample, thinstep, handle)
  !
  !    Input parameters:
  !      sigmafile   (character=*) :: FITS file containing sigma_l samples in units of l(l+1)/2pi 
  !      clfile      (character=*) :: ASCII file containing (l, C_l) in units of l(l+1)/2pi;
  !                                   used only for initialization (ie., offset estimation) to avoid
  !                                   numerical errors during exponentiation
  !      lmin        (i4b)         :: Lowest multipole to include in BR estimator
  !      lmax        (i4b)         :: Highest multipole to include in BR estimator
  !      firstchain  (i4b)         :: First chain to include in BR estimator
  !      lastchain   (i4b)         :: Last chain to include in BR estimator
  !      firstsample (i4b)         :: First sample to include in BR estimator; allow for burn-in if necessary
  !      lastsample  (i4b)         :: Last sample to include in BR estimator
  !      thinstep    (i4b)         :: Sample thinning factor -- typically 1
  ! 
  !    Output parameters:
  !      handle      (i4b)         :: this code supports multiple BR objects (optional; default=1), and
  !                                   the handle identifies specific BR data sets
  !
  ! 
  !  ** function comm_br_compute_lnL(cls, handle)
  !
  !    Input parameters:
  !      cls         (dp)          :: array containing at least cls(2:lmax) in units of l(l+1)/2pi
  !      handle      (i4b)         :: handle selecting which BR estimator to use (optional; default = 1)
  !    
  !  
  !  ** subroutine comm_br_deallocate_object(handle)
  ! 
  !    Input parameters:
  !      handle      (i4b)         :: handle of object to deallocate (optional; default=remove all)
  !
  ! ================================================================================================

  integer, parameter, private :: i4b = selected_int_kind(9)
  integer, parameter, private :: sp  = selected_real_kind(5,30)
  integer, parameter, private :: dp  = selected_real_kind(12,200)
  integer, parameter, private :: lgt = kind(.true.)

  type comm_br_data
     integer(i4b) :: lmin, lmax, numsamples, numchain
     integer(i4b) :: firstchain, lastchain, firstsample, lastsample, thinstep
     real(dp)     :: offset
     real(dp), allocatable, dimension(:,:,:) :: sigmas
  end type comm_br_data

  integer(i4b),       parameter       :: N_br = 100
  type(comm_br_data), dimension(N_br) :: comm_br
  

contains

  ! Initialization routines
  subroutine comm_br_initialize_object(sigmafile, clfile, lmin, lmax, firstchain, lastchain, &
       & firstsample, lastsample, thinstep, handle)
    implicit none

    integer(i4b),      intent(in)  :: lmin, lmax, firstchain, lastchain, firstsample, lastsample, thinstep
    character(len=*),  intent(in)  :: sigmafile, clfile
    integer(i4b),      intent(out), optional :: handle

    real(dp)           :: cl
    integer(i4b)       :: i, j, l, q, unit, numsamples, numchains, lmax_chain, id
    character(len=512) :: line
    real(dp), allocatable, dimension(:) :: cls

    id = 1
    if (present(handle)) then
       do while (allocated(comm_br(id)%sigmas) .and. id < N_br)
          id = id+1
       end do
       handle = id
    end if
    if (id < 1 .or. id > N_br) stop 'Error -- planck_br_mod: BR id out of range'

    ! Initialize BR object
    unit = comm_br_getlun()
    if (allocated(comm_br(id)%sigmas)) deallocate(comm_br(id)%sigmas)
    call comm_br_read_chain(sigmafile, unit, lmax_chain, numchains, numsamples, comm_br(id)%sigmas)
    if (lmax_chain < lmax .or. numchains < lastchain .or. numsamples < lastsample) then
       write(*,*) 'Error -- comm_br_mod: Requested parameters inconsistent with sigma file content'
       write(*,*) ''
       write(*,*) '    sigmafile  = ', trim(sigmafile)
       write(*,*) '    lmax_chain = ', lmax_chain, ' vs requested lmax       = ', lmax
       write(*,*) '    numchains  = ', numchains,  ' vs requested lastchain  = ', lastchain
       write(*,*) '    numsamples = ', numsamples, ' vs requested lastsample = ', lastsample
       stop
    end if
    comm_br(id)%lmin        = lmin
    comm_br(id)%lmax        = lmax
    comm_br(id)%firstsample = firstsample
    comm_br(id)%lastsample  = lastsample
    comm_br(id)%firstchain  = firstchain
    comm_br(id)%lastchain   = lastchain
    comm_br(id)%thinstep    = thinstep

    do j = comm_br(id)%firstsample, comm_br(id)%lastsample
       do i = comm_br(id)%firstchain, comm_br(id)%lastchain
          do l = comm_br(id)%lmin, comm_br(id)%lmax
             if (comm_br(id)%sigmas(l,i,j) <= 0.0) then
                write(*,*) 'Error -- comm_br_mod: Negative sigma_l = ', comm_br(id)%sigmas(l,i,j)
                write(*,*) '   l      = ', l
                write(*,*) '   chain  = ', i
                write(*,*) '   sample = ', j
                stop
             endif
          enddo
       enddo
    enddo

    ! Compute offset relative to input C_l spectrum
    allocate(cls(2:lmax))
    open(unit,file=trim(clfile))
    do while (.true.)
       read(unit,'(a)',end=99) line
       line = trim(line)
       if (line(1:1) == '#') cycle
       read(line,*) l, cl
       if (l >= lmin .and. l <= lmax) cls(l) = cl
    end do
99  close(unit)
    
    comm_br(id)%offset = comm_br_compute_offset(cls(2:lmax), handle)

  end subroutine comm_br_initialize_object


  subroutine comm_br_deallocate_object(handle)
    implicit none
    
    integer(i4b), optional :: handle

    integer(i4b) :: i, id_min, id_max

    id_min = 1; id_max = N_br; 
    if (present(handle)) then 
       id_min = handle
       id_max = handle
    end if

    do i = id_min, id_max
       if (allocated(comm_br(i)%sigmas)) deallocate(comm_br(i)%sigmas)
       comm_br(i)%lmin        = -1
       comm_br(i)%lmax        = -1
       comm_br(i)%firstchain  = -1
       comm_br(i)%lastchain   = -1
       comm_br(i)%firstsample = -1
       comm_br(i)%lastsample  = -1
       comm_br(i)%thinstep    =  0
    end do

  end subroutine comm_br_deallocate_object



  ! Base computation routine
  function comm_br_compute_lnL(cls, handle)
    implicit none

    real(dp),     dimension(2:), intent(in)           :: cls
    real(dp)                                          :: comm_br_compute_lnL
    integer(i4b),                intent(in), optional :: handle

    integer(i4b) :: i, j, l, id
    real(dp)     :: subtotal, x, lnL

    id = 1; if (present(handle)) id = handle    
    if (.not. allocated(comm_br(id)%sigmas)) then
       write(*,*) 'Error -- comm_br_mod: Requested BR handle ', id, ' is not allocated'
       stop
    end if

    ! Compute the Blackwell-Rao estimator
    lnL = 0.d0
    do j = comm_br(id)%firstsample, comm_br(id)%lastsample, comm_br(id)%thinstep
       do i = comm_br(id)%firstchain, comm_br(id)%lastchain
          subtotal = 0.d0
          do l = comm_br(id)%lmin, comm_br(id)%lmax
             x = comm_br(id)%sigmas(l,i,j)/cls(l)
             subtotal = subtotal + &
                  & 0.5d0 * real(2*l+1,dp) * (-x + log(x)) - log(comm_br(id)%sigmas(l,i,j))
          end do
          lnL = lnL + exp(subtotal-comm_br(id)%offset)
       end do
    end do

    if (lnL > 1e-20) then
       lnL = log(lnL)
    else
       lnL = log(1e-30)
    end if

    comm_br_compute_lnL = lnL

  end function comm_br_compute_lnL


  ! Utility routine for initializing the offset to be subtracted from each term
  ! to avoid overflow errors. Only called with the first power spectrum
  function comm_br_compute_offset(cls, handle)
    implicit none

    real(dp),     dimension(2:), intent(in)           :: cls
    real(dp)                                          :: comm_br_compute_offset
    integer(i4b),                intent(in), optional :: handle

    integer(i4b) :: i, j, l, id
    real(dp)     :: subtotal, x, offset

    id = 1; if (present(handle)) id = handle    

    ! Compute the Blackwell-Rao estimator
    offset = -1.6375d30
    do j = comm_br(id)%firstsample, comm_br(id)%lastsample, comm_br(id)%thinstep
       do i = comm_br(id)%firstchain, comm_br(id)%lastchain
          subtotal = 0.d0
          do l = comm_br(id)%lmin, comm_br(id)%lmax
             x = comm_br(id)%sigmas(l,i,j)/cls(l)
             subtotal = subtotal + &
                  & 0.5d0 * real(2*l+1,dp) * (-x + log(x)) - log(comm_br(id)%sigmas(l,i,j))
          end do
          offset = max(offset,subtotal)
       end do
    end do

    if (offset < -1.637e30) then
       write(*,*) 'Error -- comm_br_mod: Offset evaluation failed'
       write(*,*) '    handle      = ', id
       write(*,*) '    lmin        = ', comm_br(id)%lmin
       write(*,*) '    lmax        = ', comm_br(id)%lmax
       write(*,*) '    firstchain  = ', comm_br(id)%firstchain
       write(*,*) '    lastchain   = ', comm_br(id)%lastchain
       write(*,*) '    firstsample = ', comm_br(id)%firstsample
       write(*,*) '    lastsample  = ', comm_br(id)%lastsample
       write(*,*) '    thinstep    = ', comm_br(id)%thinstep
       write(*,*) ''
       write(*,*) '    cls         = ', cls(comm_br(id)%lmin:comm_br(id)%lmax)
       write(*,*) '    offset      = ', offset
       stop
    endif

    comm_br_compute_offset = offset

  end function comm_br_compute_offset

  function comm_br_getlun()
    implicit none
    integer(i4b) :: comm_br_getlun
    logical(lgt) :: exists, isopen
    comm_br_getlun = 9
    do
       comm_br_getlun = comm_br_getlun+1
       inquire(unit=comm_br_getlun,exist=exists)
       if(exists) then
          inquire(unit=comm_br_getlun,opened=isopen)
          if(.not. isopen) return
       end if
    end do
  end function comm_br_getlun
  
  ! Routine for reading the Gibbs sigma samples 
  subroutine comm_br_read_chain(filename, unit, lmax, numchains, numsamples, data)
    implicit none

    character(len=*),                            intent(in)  :: filename
    integer(i4b),                                intent(in)  :: unit
    integer(i4b),                                intent(out) :: lmax, numchains, numsamples
    real(dp),         allocatable, dimension(:,:,:)          :: data

    integer(i4b)         :: l, status, blocksize, readwrite, numspec, i, j, k
    integer(i4b)         :: fpixel, group, numargs
    logical(lgt)         :: anyf
    real(dp)             :: nullval
    character(len=80)    :: comment

    integer(i4b),          dimension(4)     :: naxes
    real(dp),     pointer, dimension(:,:,:,:) :: indata

    status = 0
    readwrite = 0
    nullval = 0.

    ! numargs = 1
    numargs = 0

    ! Open the result file
    call ftopen(unit,trim(filename),readwrite,blocksize,status)

    ! Read keywords
    call ftgkyj(unit,'LMAX',     lmax,       comment,status)
    call ftgkyj(unit,'NUMSAMP',  numsamples, comment,status)
    call ftgkyj(unit,'NUMCHAIN', numchains,  comment,status)
    call ftgkyj(unit,'NUMSPEC',  numspec,    comment,status)

    allocate(data(0:lmax,numchains,numsamples))
    nullify(indata)
    allocate(indata(0:lmax,1:1,1:numchains,1:numargs+numsamples))

    ! Read the binned power spectrum array
    group  = 1
    fpixel = 1
    call ftgpvd(unit,group,fpixel,size(indata),nullval,indata,anyf,status)
    call ftclos(unit,status)

    do k = numargs+1, numargs+numsamples
       do j = 1, numchains
          do i = 0, lmax
             data(i, j, k) = indata(i, 1, j, k)
          enddo
       enddo
    enddo
    deallocate(indata)

  end subroutine comm_br_read_chain

end module comm_br_mod