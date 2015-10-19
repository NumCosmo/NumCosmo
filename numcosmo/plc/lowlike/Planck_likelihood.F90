! ===========================================================================
MODULE planck_likelihood

! This code is the central likelihood routine from which the other subroutines
! are called.
!
! Parameters are defined in associated Planck_options module
!
! This code is based on code used by the WMAP team
! ===========================================================================


  USE Planck_options
  use Planck_teeebb_lowl

  character(len=*), parameter, public :: planck_lowlike_version='v2.1'

  logical :: initialise_pass2=.true.


PRIVATE

	public :: planck_lowlike_init
	public :: planck_lowlike_compute

contains

! ===========================================================================
SUBROUTINE planck_lowlike_init 
! ===========================================================================

  IMPLICIT NONE

  INTEGER  :: l,ll
  integer  :: il,ill,i,j
  REAL(8)  :: dummy
  LOGICAL  :: good
  integer ::  lun
  real :: rtmp, rtmp2

  print *, 'Initializing Planck low-likelihood, version '//planck_lowlike_version

!-----------------------------------------------
! initialise low l codes
!-----------------------------------------------
  if(use_lowl_pol)then
!        write(*,*) use_lowl_pol, 'using low ell'
        call teeebb_lowl_like_setup
 endif

  initialise_pass2 = .false.

END SUBROUTINE

! ===========================================================================

SUBROUTINE planck_lowlike_compute(cltt,clte,clee,clbb,like)
! ===========================================================================

  IMPLICIT NONE
  REAL(8), intent(in) :: cltt(2:*), clte(2:*), clee(2:*), clbb(2:*)
 REAL(8),intent(out) :: like(num_pl)
  INTEGER :: il, ill
  REAL(8) :: dlnlike_tot, dlnlike 
  REAL(8) :: cltt_temp(2:lowl_max)
  real, allocatable,dimension(:) :: lowl_cl
  integer :: tt_hi_l_start, te_hi_l_start
  REAL(8) :: correlation_coefficient_cl,tol=1d-10


  if(initialise_pass2)then
     call planck_lowlike_init
  endif

  Like = 0.d0

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! allocate memory

  !--------------------------------------
  ! Are cltt, clte, and clee consistent?
  !--------------------------------------
  do il = 2,lowl_max
    if ( abs(clte(il)) > 0d0 ) then
      correlation_coefficient_cl = abs(clte(il))/sqrt(cltt(il)*clee(il))
      if( correlation_coefficient_cl-1d0 > tol ) then
      write(*,*) 'TE inconsistent'
         Like(lowllike)=1d30	
         return
      end if
    end if
  enddo

  cltt_temp(2:lowl_max)=cltt(2:lowl_max)

  !---------------------------------------------------------------------------
  ! low l TE/EE/BB likelihood
  !---------------------------------------------------------------------------
  if(use_lowl_pol)then
     call teeebb_lowl_likelihood(23,cltt_temp,clte,clee,clbb,Like(lowllike),Like(lowldet))
     te_hi_l_start = 24
  else
     te_hi_l_start = temin
  endif


10 continue

end SUBROUTINE

END MODULE

