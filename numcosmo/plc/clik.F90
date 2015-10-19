module clik

  integer(kind=4), parameter :: PN_SIZE=256
  integer(kind=4), parameter :: MAX_NUMNAMES=1000
  
  type clik_object
    integer(kind=8)::ptr
  end type clik_object


contains
  
  subroutine clik_get_version(clikid,version)
    type(clik_object), intent(in) :: clikid
    character(len=*), intent(out) :: version
    
    call fortran_clik_get_version(clikid,version)  
  end subroutine

  subroutine clik_init(clikid,hdffilepath)
    
    ! Input
    type(clik_object), intent(out) :: clikid
    character(len=*), intent(in) :: hdffilepath
    ! Local
    integer(kind=4) :: fpathlen
    fpathlen = len_trim(hdffilepath)
    ! Call wrapping routine
    call fortran_clik_init(clikid%ptr,trim(hdffilepath),fpathlen)
    
  end subroutine clik_init
  
  subroutine clik_get_has_cl(clikid,has_cl)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4), dimension(6), intent(out) :: has_cl
    ! Call wrapping routine
    call fortran_clik_get_has_cl(clikid%ptr,has_cl)

  end subroutine clik_get_has_cl

  integer(kind=4) function clik_get_extra_parameter_names(clikid,names)

    ! Input
    type(clik_object), intent(in) :: clikid
    character(len=256), dimension(:), pointer :: names
    integer(kind=4) :: numnames

    ! Local
    character(len=PN_SIZE*MAX_NUMNAMES) :: buf_names

    ! Get number of extra parameters
    call fortran_clik_get_extra_parameter_number(clikid%ptr,numnames)
    if (numnames > 0) then
       ! Allocate character buffer
       allocate(names(numnames))
       ! Get the numnames, and copy into the names array
       call fortran_clik_get_extra_parameter_names(clikid%ptr,buf_names)
				do i=1,numnames
          names(i) = buf_names(1+(i-1)*PN_SIZE:(i)*PN_SIZE)
       enddo
    endif

    clik_get_extra_parameter_names=numnames

  end function clik_get_extra_parameter_names

  subroutine clik_get_lmax(clikid,lmax)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4), dimension(6), intent(out) :: lmax
    
    call fortran_get_lmax(clikid%ptr,lmax)

  end subroutine clik_get_lmax

  real(kind=8) function clik_compute(clikid,cl_and_pars)

    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:) :: cl_and_pars
    ! Local 
    real(kind=8) :: lkl
    call fortran_clik_compute(clikid%ptr,cl_and_pars,lkl)

    clik_compute = lkl
    return

  end function clik_compute
  
  real(kind=8) function clik_compute_with_error(clikid,cl_and_pars,ler)

    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:) :: cl_and_pars
    logical,intent(inout)::ler
    ! Local 
    real(kind=8) :: lkl
    integer(kind=4)::ier
    ier = 0
    call fortran_clik_compute_with_error(clikid%ptr,cl_and_pars,lkl,ier)
    if (ier.eq.0) then
      ler = .FALSE.
    else
      ler = .TRUE.
    endif

    clik_compute_with_error = lkl
    return

  end function clik_compute_with_error

  subroutine clik_cleanup(clikid)

    !Input
    type(clik_object), intent(in) :: clikid
    call fortran_clik_cleanup(clikid%ptr)

  end subroutine clik_cleanup
  
  subroutine clik_fill_cl(clikid,cl_and_pars,clTT,clEE,clBB,clTE,clTB,clEB,flag_cl_start_at_two,flag_is_llp1_over_2pi,extrapars)

    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:), pointer,intent(out) :: cl_and_pars
    real(kind=8), dimension(:), intent(in) :: clTT,clEE,clBB,clTE,clTB,clEB,extrapars
    logical,intent(in) ::  flag_cl_start_at_two,flag_is_llp1_over_2pi
    integer :: l, offl, curl,cli
    real(kind=8):: llp1_over_2pi,over_2pi
    integer(kind=4), dimension(6) :: lmax
    integer(kind=4) :: nextra, ntot, clmax

    over_2pi = 1./(6.283185307179586476925286766559005768394) 
    
    if (flag_cl_start_at_two) then
      offl = 2
    else
      offl = 0
    end if

    call clik_get_lmax(clikid,lmax)
    if (size(clTT)+offl<lmax(1)+1) then
      stop "not enough data in clTT"
    end if
    if (size(clEE)+offl<lmax(2)+1) then
      stop "not enough data in clTT"
    end if
    if (size(clBB)+offl<lmax(3)+1) then
      stop "not enough data in clBB"
    end if
    if (size(clTE)+offl<lmax(4)+1) then
      stop "not enough data in clTE"
    end if
    if (size(clTB)+offl<lmax(5)+1) then
      stop "not enough data in clTB"
    end if
    if (size(clEB)+offl<lmax(6)+1) then
      stop "not enough data in clEB"
    end if
    
    call fortran_clik_get_extra_parameter_number(clikid%ptr,nextra)
    
    if (size(extrapars)<nextra) then
      stop "not enough data in extrapars"
    end if
    
    ntot = 6 + nextra
    do cli=1,6
      ntot = ntot+lmax(cli)

    end do

    allocate(cl_and_pars(ntot))

    curl = 1
    
    !TT
    clmax = lmax(1)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTT(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if
    
    !EE
    clmax = lmax(2)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clEE(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !BB
    clmax = lmax(3)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clBB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !TE
    clmax = lmax(4)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTE(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !TB
    clmax = lmax(5)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clTB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    !EB
    clmax = lmax(6)
    if (clmax .ne. -1) THEN
      do l = 0,offl
        cl_and_pars(curl) = 0
        curl = curl + 1
      end do
      do l = offl,clmax
        llp1_over_2pi = 1.
        if (flag_is_llp1_over_2pi) THEN
          llp1_over_2pi = l*(l+1.) * over_2pi
        end if
        cl_and_pars(curl) = clEB(l+1-offl) / llp1_over_2pi 
        curl = curl + 1
      end do
    end if

    ! copy extra parameters
    cl_and_pars(curl:ntot) = extrapars(:nextra)

  end subroutine  clik_fill_cl  

!#ifdef CLIK_LENSING
  subroutine clik_lensing_init(clikid,hdffilepath)
    
    ! Input
    type(clik_object), intent(out) :: clikid
    character(len=*), intent(in) :: hdffilepath
    ! Local
    integer(kind=4) :: fpathlen
    fpathlen = len_trim(hdffilepath)
    ! Call wrapping routine
    call fortran_clik_lensing_init(clikid%ptr,trim(hdffilepath),fpathlen)
    
  end subroutine clik_lensing_init
  
  subroutine clik_lensing_get_lmax(clikid,lmax)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4), intent(out) :: lmax
    
    call fortran_clik_lensing_get_lmax(lmax,clikid%ptr)

  end subroutine clik_lensing_get_lmax

  subroutine clik_lensing_get_lmaxs(clikid,lmax)

    ! Input
    type(clik_object), intent(in) :: clikid
    integer(kind=4),dimension(7),intent(out) :: lmax
    
    call fortran_clik_lensing_get_lmaxs(clikid%ptr,lmax)

  end subroutine clik_lensing_get_lmaxs

  integer(kind=4) function clik_lensing_get_extra_parameter_names(clikid,names)

    ! Input
    type(clik_object), intent(in) :: clikid
    character(len=256), dimension(:), pointer :: names
    integer(kind=4) :: numnames

    ! Local
    character(len=PN_SIZE*MAX_NUMNAMES) :: buf_names

    ! Get number of extra parameters
    call fortran_clik_lensing_get_extra_parameter_number(clikid%ptr,numnames)
    if (numnames > 0) then
       ! Allocate character buffer
       allocate(names(numnames))
       ! Get the numnames, and copy into the names array
       call fortran_clik_lensing_get_extra_parameter_names(clikid%ptr,buf_names)
        do i=1,numnames
          names(i) = buf_names(1+(i-1)*PN_SIZE:(i)*PN_SIZE)
       enddo
    endif

    clik_lensing_get_extra_parameter_names=numnames

  end function clik_lensing_get_extra_parameter_names

  real(kind=8) function clik_lensing_compute(clikid,cl_and_pars)

    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8), dimension(:) :: cl_and_pars
    ! Local 
    real(kind=8) :: lkl
    call fortran_clik_lensing_compute(lkl,clikid%ptr,cl_and_pars)
    clik_lensing_compute = lkl
    return

  end function clik_lensing_compute
  
  subroutine clik_lensing_cleanup(clikid)

    !Input
    type(clik_object), intent(in) :: clikid
    call fortran_clik_lensing_cleanup(clikid%ptr)

  end subroutine clik_lensing_cleanup

  subroutine clik_lensing_cltt_fid(clikid,cltt)
    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8),dimension(:),allocatable,intent(out)::cltt
    integer,dimension(7)::lmax

    call clik_lensing_get_lmaxs(clikid,lmax)
    allocate(cltt(0:lmax(2)))
    call fortran_clik_lensing_cltt_fid(clikid,cltt)
  end subroutine clik_lensing_cltt_fid

  subroutine clik_lensing_clcmb_fid(clikid,cltt)
    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8),dimension(:),allocatable,intent(out)::cltt
    integer,dimension(7)::lmax
    integer::ntot,i

    call clik_lensing_get_lmaxs(clikid,lmax)
    ntot = lmax(2)
    do i=3,7
      ntot = ntot + lmax(i)+1
    enddo

    allocate(cltt(0:ntot-1))
    call fortran_clik_lensing_clcmb_fid(clikid,cltt)
  end subroutine clik_lensing_clcmb_fid

  subroutine clik_lensing_clpp_fid(clikid,cltt)
    !Input
    type(clik_object), intent(in) :: clikid
    real(kind=8),dimension(:),allocatable,intent(out)::cltt
    integer,dimension(7)::lmax

    call clik_lensing_get_lmaxs(clikid,lmax)
    allocate(cltt(0:lmax(1)))
    call fortran_clik_lensing_clpp_fid(clikid,cltt)
  end subroutine clik_lensing_clpp_fid

  subroutine clik_try_lensing(is_lensing,hdffilepath)
    
    ! Input
    logical,intent(out)::is_lensing
    character(len=*), intent(in) :: hdffilepath
    ! Local
    integer(kind=4) :: fpathlen
    integer(kind=4):: isl

    fpathlen = len_trim(hdffilepath)
    ! Call wrapping routine
    call fortran_clik_try_lensing(isl,trim(hdffilepath),fpathlen)
    
    is_lensing = .FALSE.
    if (isl==1) then
      is_lensing = .TRUE.
    endif

  end subroutine clik_try_lensing
!#endif

end module clik
