MODULE LOWLIKE_EXTRA

	IMPLICIT NONE

	INTEGER:: BOK = 0
	INTEGER:: CLIK_LMAX,CLIK_LMIN
	real(8), dimension(:), allocatable :: cltt,clte,clee,clbb
	
END MODULE LOWLIKE_EXTRA



SUBROUTINE LOWLIKE_EXTRA_ONLY_ONE(MOK)
	USE LOWLIKE_EXTRA
	INTEGER,INTENT(OUT)::MOK
	MOK = BOK
	BOK = 1
END SUBROUTINE 	LOWLIKE_EXTRA_ONLY_ONE

SUBROUTINE LOWLIKE_EXTRA_FREE()
	USE LOWLIKE_EXTRA
	BOK =0
	deallocate(cltt)
	deallocate(clte)
	deallocate(clee)
	deallocate(clbb)
END SUBROUTINE 	LOWLIKE_EXTRA_FREE

SUBROUTINE LOWLIKE_EXTRA_LKL(LKL,CL)
	USE LOWLIKE_EXTRA
	use planck_likelihood
	use planck_options
  	use healpix_types

	REAL(8),INTENT(OUT)::LKL
	REAL(8) :: like(num_pl)
	REAL(8),INTENT(IN),DIMENSION(0:4*CLIK_LMAX+3)::CL
	INTEGER::i,cur

	!TT
	cur = 0
	cltt = 0
	clee = 0
	clbb = 0
	clte = 0

	DO i = ttmin,ttmax
		cltt(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1	

	!EE	
	DO i = temin,temax
		clee(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1

	!BB
	DO i = temin,temax
		clbb(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1
	
	!TE
	DO i = temin,temax
		clte(i)=CL(cur+i)*(i*(i+1.))/TWOPI
	END DO	
	cur = cur+clik_lmax+1
	
	CALL planck_lowlike_compute(cltt,clte,clee,clbb,like)
	!do i =1,num_pl
	!	write(*,*) i,like(i)
	!enddo	

	LKL = -sum(like(1:num_pl))
END SUBROUTINE 	LOWLIKE_EXTRA_LKL

SUBROUTINE LOWLIKE_EXTRA_PARAMETER_INIT(tt_min,tt_max,te_min,te_max,m_use_gibbs,m_use_lowl_pol,m_use_wmap_pol)
	USE LOWLIKE_EXTRA
	use planck_likelihood
	use planck_options
	INTEGER,INTENT(IN)::tt_min,tt_max,te_min,te_max,m_use_gibbs,m_use_lowl_pol
	
 	ttmin = tt_min
	ttmax = tt_max
	temin = te_min
	temax = te_max
	lowl_max = 32
	
	!if (ttmin>lowl_max) then
	!	use_lowl_TT = .false.
	!endif
	
	use_lowl_pol = .false.
	if (m_use_lowl_pol==1) then
		use_lowl_pol = .true.
	endif
	
	use_wmap_pol = .false.
	if (m_use_wmap_pol==1) then
		use_wmap_pol = .true.
	endif
	
	clik_lmax = tt_max
	if (te_max>clik_lmax) then
		clik_lmax = te_max
	endif
	clik_lmin = tt_min
	if (te_max<clik_lmin) then
		clik_lmin = te_min
	endif
	
	allocate( cltt(2:clik_lmax) )
	allocate( clte(2:clik_lmax) )
	allocate( clee(2:clik_lmax) )
	allocate( clbb(2:clik_lmax) )
	
	CALL planck_lowlike_init()

END SUBROUTINE 	LOWLIKE_EXTRA_PARAMETER_INIT
