! ===========================================================================
MODULE PLANCK_OPTIONS

! This module contains the options in the likelihood code
!
! ===========================================================================

!---------------------------------------------------
! location of input data
! ---------------------------------------------------
	 
  character(len=*), parameter :: Planck_data_dir = 'data/'


!---------------------------------------------------
! likelihood terms from Planck
!---------------------------------------------------
  integer,parameter   :: num_pl   = 4    ! number of individual chi2 terms in likelihood

  integer,parameter   :: ttlowllike = 1    ! low tttt chisq flag
  integer,parameter   :: ttlowldet  = 2    ! low tttt determinant flag
  integer,parameter   :: lowllike   = 3    ! TE/EE/BB lowl chisq flag
  integer,parameter   :: lowldet    = 4    ! TE/EE/BB lowl determinant flag

!---------------------------------------------------
! l range to be used in the likelihood code
!---------------------------------------------------
  integer :: ttmin                = 2    ! must be l.ge.2 
  integer :: temin                = 2    ! must be l.ge.2 
   integer :: ttmax                = 32    ! must be l.ge.2
  integer :: temax                = 32    ! must be l.ge.2


!---------------------------------------------------
! various likelihood options
! change these to include/ exclude various likelihood aspects
!---------------------------------------------------
  logical :: use_lowl_pol         = .true. ! include TE,EE,BB pixel likelihood for l<24

  logical :: use_wmap_pol	  = .false. !Use WMAP for low-ell pol

  integer :: lowl_max             = 32     ! use low l TT code 2<l<lowl_max

 double precision, parameter :: teeebb_pixlike_lndet_offset &
                = 16078.083180d0


contains

  subroutine planck_print_options()
    print *, "-----------------------------------------------------"
    print *, "Planck_data_dir = ", trim(Planck_data_dir) 
    print *, ""
    print *, "ttmin = ", ttmin 
    print *, "temin = ", temin
    print *, ""
    print *, "use_lowl_pol =      ", use_lowl_pol         
    print *, ""
    print *, "lowl_tt_res = ", lowl_tt_res
    print *, "lowl_max =    ", lowl_max
    print *, ""
  end subroutine planck_print_options
  
  subroutine get_free_lun( lun )

        implicit none

        integer, intent(out) :: lun

        integer, save :: last_lun = 19
        logical :: used

        lun = last_lun
        do
                inquire( unit=lun, opened=used )
                if ( .not. used ) exit
                lun = lun + 1
        end do

        last_lun = lun

  end subroutine


END MODULE PLANCK_OPTIONS		
