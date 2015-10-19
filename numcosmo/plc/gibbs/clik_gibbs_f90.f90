MODULE GIBBS_EXTRA

	IMPLICIT NONE

	INTEGER,dimension(100):: CLIK_LMAX,CLIK_LMIN,clik_approx_chi2
	real(8),dimension(2:1000) :: cltt
	real(8),parameter:: PI    = 3.141592653589793238462643383279502884197

END MODULE GIBBS_EXTRA




SUBROUTINE GIBBS_EXTRA_FREE(handle)
	use comm_br_mod
	use GIBBS_EXTRA
	INTEGER,intent(in)::handle
	
	call comm_br_deallocate_object(handle)
	
	clik_lmin(handle) = -1
	clik_lmax(handle) = -1
	
END SUBROUTINE 	GIBBS_EXTRA_FREE

SUBROUTINE GIBBS_EXTRA_LKL(LKL,handle,CL)
	use comm_br_mod
	use GIBBS_EXTRA
	
	REAL(8),INTENT(OUT)::LKL
	INTEGER,intent(in)::handle
	REAL(8),INTENT(IN),DIMENSION(0:CLIK_LMAX(handle)-CLIK_LMIN(handle))::CL
	INTEGER::i,cur

	!TT
	cur = 0
	cltt = 0
	DO i = clik_lmin(handle),clik_lmax(handle)
		cltt(i)=CL(cur)*(i*(i+1.))/2./PI
		cur = cur + 1
	END DO	

	LKL = comm_br_compute_lnL(cltt(2:clik_lmax(handle)),handle)		
	
	if (clik_approx_chi2(handle)==1) then
		LKL = LKL - 0.5d0 * (clik_lmax(handle)-clik_lmin(handle)+1)
	endif

END SUBROUTINE 	GIBBS_EXTRA_LKL



SUBROUTINE GIBBS_EXTRA_PARAMETER_INIT(handle,datadir,l_datadir,lmin,lmax,firstchain,lastchain,firstsample,lastsample,step,approx_chi2)
	use comm_br_mod
	use GIBBS_EXTRA
	
	INTEGER,INTENT(IN)::l_datadir,lmin,lmax,firstchain,lastchain,firstsample,lastsample,step,approx_chi2
	character(len=l_datadir)::datadir
	INTEGER,INTENT(OUT)::handle
	character(len=1000)::sigma_file,cl_file,data_dir

		
	data_dir = TRIM(datadir)
	sigma_file = TRIM(data_dir)//"/sigma.fits"
	cl_file = TRIM(data_dir)//"/cl.dat"

	handle = 0

	call comm_br_initialize_object(sigma_file, cl_file, lmin, lmax, firstchain, lastchain, &
       & firstsample, lastsample, step, handle)
	
	clik_lmin(handle) = lmin
	clik_lmax(handle) = lmax
	clik_approx_chi2(handle) = approx_chi2
		
END SUBROUTINE 	GIBBS_EXTRA_PARAMETER_INIT

SUBROUTINE GIBBS_GAUSS_EXTRA_FREE(handle)
	use comm_gauss_br_mod
	use GIBBS_EXTRA
	INTEGER,intent(in)::handle
	
	call comm_gauss_br_deallocate_object(handle)
	
	clik_lmin(handle) = -1
	clik_lmax(handle) = -1
	
END SUBROUTINE 	GIBBS_GAUSS_EXTRA_FREE

SUBROUTINE GIBBS_GAUSS_EXTRA_LKL(LKL,handle,CL)
	use comm_gauss_br_mod
	use GIBBS_EXTRA
	
	REAL(8),INTENT(OUT)::LKL
	INTEGER,intent(in)::handle
	REAL(8),INTENT(IN),DIMENSION(0:CLIK_LMAX(handle)-CLIK_LMIN(handle))::CL
	INTEGER::i,cur

	!TT
	cur = 0
	cltt = 0
	DO i = clik_lmin(handle),clik_lmax(handle)
		cltt(i)=CL(cur)*(i*(i+1.))/2./PI
		cur = cur + 1
	END DO	

	LKL = comm_gauss_br_compute_lnL(cltt(2:clik_lmax(handle)),handle)		
	
END SUBROUTINE 	GIBBS_GAUSS_EXTRA_LKL



SUBROUTINE GIBBS_GAUSS_EXTRA_PARAMETER_INIT(handle,lmin,lmax,delta_l)
	use comm_gauss_br_mod
	use GIBBS_EXTRA
	
	INTEGER,INTENT(IN)::lmin,lmax,delta_l
	INTEGER,INTENT(OUT)::handle
	character(len=1000)::sigma_file

		
	sigma_file = "sigma.fits"
	
	handle = 0
	
	call comm_gauss_br_initialize_object(sigma_file, lmin, lmax, delta_l, handle)
	
	clik_lmin(handle) = lmin
	clik_lmax(handle) = lmax
		
END SUBROUTINE 	GIBBS_GAUSS_EXTRA_PARAMETER_INIT

MODULE COMM_LOWL_EXTRA

	IMPLICIT NONE

	INTEGER,dimension(100):: CLIK_LMAX, CLIK_LMIN
	real(8),dimension(0:100,1:6) :: cl_
	real(8),parameter:: PI    = 3.141592653589793238462643383279502884197

END MODULE COMM_LOWL_EXTRA

SUBROUTINE COMM_LOWL_EXTRA_PARAMETER_INIT(handle,parfile,l_parfile,lmin,lmax)
use comm_lowl_mod_dist
use COMM_LOWL_EXTRA
character(len=l_parfile)::parfile
integer,intent(OUT)::handle
integer,intent(IN)::l_parfile,lmax

call comm_lowl_initialize_object(parfile,handle)
clik_lmax(handle) = lmax
clik_lmin(handle) = lmin

END SUBROUTINE COMM_LOWL_EXTRA_PARAMETER_INIT

SUBROUTINE COMM_LOWL_EXTRA_LKL(LKL,handle,CL)
use comm_lowl_mod_dist
use COMM_LOWL_EXTRA
REAL(8),INTENT(OUT)::LKL
INTEGER,intent(in)::handle
REAL(8),INTENT(IN),DIMENSION(0:(CLIK_LMAX(handle)-clik_lmin(handle)+1)*6-1)::CL
integer::zero,i

cl_ = 0

!TT
zero = 0
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,1) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

!EE
zero = zero + clik_lmax(handle)+1-clik_lmin(handle)
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,4) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

!BB
zero = zero + clik_lmax(handle)+1-clik_lmin(handle)
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,6) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

!TE
zero = zero + clik_lmax(handle)+1-clik_lmin(handle)
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,2) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

!TB
zero = zero + clik_lmax(handle)+1-clik_lmin(handle)
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,3) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

!EB
zero = zero + clik_lmax(handle)+1-clik_lmin(handle)
DO i = clik_lmin(handle),clik_lmax(handle)
	cl_(i,5) = CL(zero+i-clik_lmin(handle)) *i*(i+1.)/2./PI
END DO	

LKL = comm_lowl_compute_lnL(cls=cl_, handle=handle)


END SUBROUTINE COMM_LOWL_EXTRA_LKL

SUBROUTINE COMM_LOWL_EXTRA_FREE(handle)
	use comm_lowl_mod_dist
	use COMM_LOWL_EXTRA
	INTEGER,intent(in)::handle
	
	call comm_lowl_deallocate_object(handle)
	
	clik_lmax(handle) = -1
	
END SUBROUTINE 	COMM_LOWL_EXTRA_FREE
