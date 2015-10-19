MODULE ACTSPT_EXTRA

	IMPLICIT NONE

	INTEGER:: BOK = 0
	INTEGER:: CLIK_LMAX,CLIK_LMIN,big_clik_lmax
	real(8), dimension(:), allocatable :: cltt
	
END MODULE ACTSPT_EXTRA



SUBROUTINE ACTSPT_EXTRA_ONLY_ONE(MOK)
	USE ACTSPT_EXTRA
	INTEGER,INTENT(OUT)::MOK
	MOK = BOK
	BOK = 1
END SUBROUTINE 	ACTSPT_EXTRA_ONLY_ONE

SUBROUTINE ACTSPT_EXTRA_FREE()
	USE ACTSPT_EXTRA
	BOK =0
	deallocate(cltt)
END SUBROUTINE 	ACTSPT_EXTRA_FREE

SUBROUTINE ACTSPT_EXTRA_LKL(LKL,CL)
	use highell_options
	use highell_likelihood
	use ACTSPT_EXTRA
	
	REAL(8),INTENT(OUT)::LKL
	REAL(8),INTENT(IN),DIMENSION(0:CLIK_LMAX+1+21-CLIK_LMIN)::CL
	INTEGER::i,cur
	real(8)::amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,ncib,rps0,rps1,rps,rcib,cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3,like_tot,ags,age

	!TT
	cur = 0
	cltt = 0
	DO i = clik_lmin,clik_lmax
		cltt(i)=CL(cur)*(i*(i+1.))/2./PI
		!if (cur<10) then
		!	print *,i,cur,CL(cur),cltt(i)
		!endif
		cur = cur + 1
	END DO	

	amp_tsz  = CL(cur)
	cur = cur + 1
	amp_ksz  = CL(cur)
	cur = cur + 1
	xi       = CL(cur)
	cur = cur + 1
	aps148   = CL(cur)
	cur = cur + 1
	aps217   = CL(cur)
	cur = cur + 1
	aps95    = CL(cur)
	cur = cur + 1
	aps150   = CL(cur)
	cur = cur + 1
	aps220   = CL(cur)
	cur = cur + 1
	acib150  = CL(cur)
	cur = cur + 1
	acib220  = CL(cur)
	cur = cur + 1
	ncib  = CL(cur)
	cur = cur + 1
	rps0     = CL(cur)
	cur = cur + 1
	rps1     = CL(cur)
	cur = cur + 1
	rps      = CL(cur)
	cur = cur + 1
	rcib     = CL(cur)
	cur = cur + 1
	ags     = CL(cur)
	cur = cur + 1
	age     = CL(cur)
	cur = cur + 1
	cas1     = CL(cur)
	cur = cur + 1
	cas2     = CL(cur)
	cur = cur + 1
	cae1     = CL(cur)
	cur = cur + 1
	cae2    = CL(cur)
	cur = cur + 1
	cal_1    = CL(cur)
	cur = cur + 1
	cal_2    = CL(cur)
	cur = cur + 1
	cal_3    = CL(cur)
	cur = cur + 1
	
	CALL highell_likelihood_compute(cltt,amp_tsz,amp_ksz,xi,aps148,aps217,aps95,aps150,aps220,acib150,acib220,ncib,rps0,rps1,rps,rcib,ags,age,cas1,cas2,cae1,cae2,cal_1,cal_2,cal_3,like_tot)
	
	LKL = -like_tot

END SUBROUTINE 	ACTSPT_EXTRA_LKL



SUBROUTINE ACTSPT_EXTRA_PARAMETER_INIT(datadir,l_datadir,ilmin11,ilmin12,ilmin22,ilmax11,ilmax12,ilmax22,itt_lmax_mc,iuse_act_south  , iuse_act_equa    , iuse_spt_highell)
	use highell_options
	use highell_likelihood
	use ACTSPT_EXTRA

	INTEGER,INTENT(IN)::l_datadir,ilmin11,ilmin12,ilmin22,ilmax11,ilmax12,ilmax22,iuse_act_south  , iuse_act_equa    , iuse_spt_highell
	character(len=l_datadir)::datadir
	
	
	data_dir = TRIM(datadir)
	ACT_data_dir = TRIM(data_dir)//"/data_act/"
	SPT_data_dir = TRIM(data_dir)//"/data_spt/"

	lmin11 = ilmin11
	lmin12 = ilmin12
	lmin22 = ilmin22
	lmax11 = ilmax11
	lmax12 = ilmax12
	lmax22 = ilmax22

	tt_lmax_mc = itt_lmax_mc

	use_spt_lowell = .false.

	use_spt_highell = .false.
	if (iuse_spt_highell.EQ.1) then
		use_spt_highell = .true.
	endif

	use_act_equa = .false.
	if (iuse_act_equa.EQ.1) then
		use_act_equa = .true.
	endif

	use_act_south = .false.
	if (iuse_act_south.EQ.1) then
		use_act_south = .true.
	endif

	clik_lmax = tt_lmax_k
	if (lmax11>clik_lmax) then
		clik_lmax = lmax11
	endif
	if (lmax12>clik_lmax) then
		clik_lmax = lmax12
	endif
	if (lmax22>clik_lmax) then
		clik_lmax = lmax22
	endif
	big_clik_lmax = clik_lmax
	if (tt_lmax_mc<clik_lmax) then
		clik_lmax = tt_lmax_mc
	endif

	
	clik_lmin = 2
	if (lmin12<clik_lmin) then
		clik_lmin = lmin12
	endif
	if (lmin12<clik_lmin) then
		clik_lmin = lmin12
	endif
	if (lmin22<clik_lmin) then
		clik_lmin = lmin22
	endif

	
	allocate( cltt(2:big_clik_lmax) )
	
	!PRINT *,lmin11,lmin12,lmin22,lmax11,lmax12,lmax22,use_act_south,use_act_equa,use_spt_highell,data_dir,ACT_data_dir,SPT_data_dir
	call highell_likelihood_init
	
END SUBROUTINE 	ACTSPT_EXTRA_PARAMETER_INIT
