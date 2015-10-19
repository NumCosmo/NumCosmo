MODULE plik_cmbonly_EXTRA

	IMPLICIT NONE

	INTEGER:: BOK = 0
	INTEGER:: CLIK_LMAX,CLIK_LMIN
	real(8), dimension(:), allocatable :: cltt,clte,clee
	real(8), parameter :: PI    = 3.14159265358979323846264d0

END MODULE plik_cmbonly_EXTRA



SUBROUTINE plik_cmbonly_EXTRA_only_ONE(MOK)
	USE plik_cmbonly_EXTRA
	INTEGER,INTENT(OUT)::MOK
	MOK = BOK
	BOK = 1
END SUBROUTINE 	plik_cmbonly_EXTRA_only_ONE

SUBROUTINE plik_cmbonly_extra_FREE()
	USE plik_cmbonly_EXTRA
	BOK =0
	deallocate(cltt)
	deallocate(clte)
	deallocate(clee)
END SUBROUTINE 	plik_cmbonly_extra_FREE

SUBROUTINE plik_cmbonly_extra_LKL(LKL,CL)
	use plik_cmbonly_EXTRA
	use Plik_CMBonly

	REAL(8),INTENT(OUT)::LKL
	REAL(8),INTENT(IN),DIMENSION(0:(CLIK_LMAX+1-CLIK_LMIN)*3)::CL
	INTEGER::i,cur
	real(8)::like_tot,calPlanck

	!TT
	cur = 0
	cltt = 0

	if (use_tt) then
		DO i = clik_lmin,clik_lmax
			cltt(i)=CL(cur)*(i*(i+1.))/2./PI
			cur = cur + 1
		END DO	
	endif

	if (use_ee) then
		DO i = clik_lmin,clik_lmax
			clee(i)=CL(cur)*(i*(i+1.))/2./PI
			cur = cur + 1
		END DO	
	endif
	if (use_te) then
		DO i = clik_lmin,clik_lmax
			clte(i)=CL(cur)*(i*(i+1.))/2./PI
			cur = cur + 1
		END DO	
	endif
	calPlanck = Cl(cur)

	call calc_like_cmbonly(like_tot,cltt,clte,clee,calPlanck)	
	LKL = -like_tot

END SUBROUTINE 	plik_cmbonly_extra_LKL



SUBROUTINE plik_cmbonly_extra_INIT(datadir,l_datadir,iuse_tt, iuse_ee, iuse_te)
	use Plik_CMBonly
	use plik_cmbonly_extra

	INTEGER,INTENT(IN)::l_datadir
	character(len=l_datadir)::datadir
	
	data_dir = TRIM(datadir)
	!write(*,*),data_dir
	clik_lmin = plmin
	clik_lmax = plmax
	tt_lmax = plmax !dans ton cul

	allocate( cltt(2:clik_lmax) )
	allocate( clee(2:clik_lmax) )
	allocate( clte(2:clik_lmax) )
	
	!PRINT *,lmin11,lmin12,lmin22,lmax11,lmax12,lmax22,use_act_south,use_act_equa,use_spt_highell,data_dir,ACT_data_dir,SPT_data_dir
	use_tt = iuse_tt.NE.0
	use_te = iuse_te.NE.0
	use_ee = iuse_ee.NE.0
	
	call like_init_cmbonly
	
	
END SUBROUTINE 	plik_cmbonly_extra_INIT
