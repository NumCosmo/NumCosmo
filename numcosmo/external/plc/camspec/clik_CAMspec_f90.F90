MODULE CAMSPEC_EXTRA

    IMPLICIT NONE

    INTEGER:: BOK = 0
    real(8), dimension(:), allocatable :: cltt,nuisance,clte,clee
    INTEGER::lmin,lmax,lminte,lmaxte,lminee,lmaxee,clik_lmin,clik_lmax
    logical,dimension(3)::l_has_cl
    INTEGER::xcase,npar,n_nuisance
    integer:: cam_version
END MODULE CAMSPEC_EXTRA



SUBROUTINE CAMSPEC_EXTRA_ONLY_ONE(MOK)
    USE CAMSPEC_EXTRA
    INTEGER,INTENT(OUT)::MOK
    MOK = BOK
    BOK = 1
END SUBROUTINE  CAMSPEC_EXTRA_ONLY_ONE

SUBROUTINE CAMSPEC_EXTRA_FREE()
    USE CAMSPEC_EXTRA
    DEALLOCATE(cltt)
    BOK =1
    ! forbid double instantiation of camspec. Need more work to cleanup correctly it seems...

END SUBROUTINE  CAMSPEC_EXTRA_FREE

SUBROUTINE CAMSPEC_EXTRA_INIT(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes,ihas_dust,ihas_calib,imarge_flag, imarge_mode,imarge_num, ikeep_num,bs_factor)
    USE CAMSPEC_EXTRA
    USE temp_like
    implicit none
    integer,intent(in)::iNspec,inX,inum_modes_per_beam,ibeam_lmax,ibeam_Nspec,icov_dim,ilmax_sz,ihas_dust,ihas_calib
    integer,dimension(1:iNspec)::ilminX,ilmaxX,inp,inpt
    real*8, dimension(1:iNx) :: iX
    real*8,  dimension(1:iNx,1:iNx) ::ic_inv
    real*8,intent(in) :: iksz_temp(0:ilmax_sz), itszxcib_temp(0:ilmax_sz)
    real*8,intent(in) :: isz_143_temp(0:ilmax_sz)
    real*8, dimension(1:icov_dim,1:icov_dim) :: ibeam_cov_inv
    real*8, dimension(1:inum_modes_per_beam,0:ibeam_lmax,1:ibeam_Nspec) :: ibeam_modes ! mode#, l, spec#
    integer, dimension(1:icov_dim)::imarge_flag
    logical,dimension(1:icov_dim)::marge_flag
    integer::imarge_num,ikeep_num
    real*8,dimension(1:imarge_num, 1:ikeep_num)::imarge_mode
    real*8::bs_factor

    integer::i
    
    do i=1,icov_dim
        marge_flag(i) = (imarge_flag(i).eq.1)
    enddo

    call like_init_frommem(iNspec, inX,ilminX,ilmaxX,inp,inpt, ic_inv,iX,ilmax_sz,isz_143_temp,iksz_temp,itszxcib_temp,ibeam_Nspec,inum_modes_per_beam,ibeam_lmax,icov_dim,ibeam_cov_inv,ibeam_modes,ihas_dust,ihas_calib,marge_flag,imarge_mode)
    
    lmax = ilmaxX(1)
    DO i=2,iNspec
        IF (lmax<ilmaxX(i)) lmax = ilmaxX(i)
    ENDDO
    lmin = ilminX(1)
    DO i=2,iNspec
        IF (lmin>ilminX(i)) lmin = ilminX(i)
    ENDDO
    
    xcase = 1
    !IF (mlmax2.ne.0) xcase = 1

    ALLOCATE(cltt(0:lmax+1))

    allocate(nuisance(num_non_beam+ihas_dust+beam_Nspec*num_modes_per_beam - marge_num))
    npar = lmax+1-lmin+num_non_beam+ beam_Nspec * num_modes_per_beam +ihas_dust - marge_num
    beam_factor = bs_factor
    cam_version =1

END SUBROUTINE CAMSPEC_EXTRA_INIT

#ifndef CAMSPEC_V1
SUBROUTINE   camspec_extra_init_v2(ipre_marged,like_file,l_like_file,sz143_file,l_sz143_file,tszxcib_file,l_tszxcib_file,ksz_file,l_ksz_file,beam_file,l_beam_file,data_vector,l_data_vector,icamspec_fiducial_foregrounds,l_camspec_fiducial_foregrounds,icamspec_fiducial_cl,l_camspec_fiducial_cl,lmins,lmaxs,spec_flag,icamspec_beam_mcmc_num,xdim,iclik_lmin,iclik_lmax,has_cl,bs_factor) 
    use CAMSPEC_EXTRA
    use temp_like_camspec
    implicit none
    integer,intent(in)::l_like_file,l_beam_file,l_tszxcib_file,l_sz143_file,l_data_vector,l_ksz_file,xdim,ipre_marged,icamspec_beam_mcmc_num,iclik_lmin,iclik_lmax,l_camspec_fiducial_foregrounds,l_camspec_fiducial_cl
    character(len=l_like_file)::like_file
    character(len=l_sz143_file)::sz143_file
    character(len=l_tszxcib_file)::tszxcib_file
    character(len=l_ksz_file)::ksz_file
    character(len=l_beam_file)::beam_file
    character(len=l_data_vector)::data_vector
    character(len=l_camspec_fiducial_foregrounds)::icamspec_fiducial_foregrounds
    character(len=l_camspec_fiducial_cl)::icamspec_fiducial_cl
    integer,intent(in),dimension(6)::lmaxs,lmins,spec_flag,has_cl
    logical::pre_marged
    integer::i
    real*8::bs_factor
    
    n_nuisance = xdim
    do i=1,6
        want_spec(i) = spec_flag(i)==1
        camspec_lmins(i) = lmins(i)
        camspec_lmaxs(i) = lmaxs(i)
    enddo

    camspec_beam_mcmc_num = icamspec_beam_mcmc_num
    
    camspec_fiducial_foregrounds = icamspec_fiducial_foregrounds
    camspec_fiducial_cl = icamspec_fiducial_cl

    pre_marged = ipre_marged==1

    call like_init(pre_marged,like_file, sz143_file, tszxcib_file, ksz_file, beam_file,data_vector)

    cam_version=2

    clik_lmin = iclik_lmin
    clik_lmax = iclik_lmax
    lmin = minval(lminX(1:4))
    lmax = maxval(lmaxX(1:4))
    lminte = lminX(5)
    lmaxte = lmaxX(5)
    lminee = lminX(6)
    lmaxee = lmaxX(6)

    ALLOCATE(cltt(0:lmax+1))
    ALLOCATE(clte(0:lmaxte+1))
    ALLOCATE(clee(0:lmaxee+1))

    allocate(nuisance(xdim))
    npar = (clik_lmax+1-clik_lmin)*has_cl(1) + (clik_lmax+1-clik_lmin)*has_cl(2) + (clik_lmax+1-clik_lmin)*has_cl(4) + xdim
    
    beam_factor = bs_factor

    l_has_cl(1) = has_cl(1)==1
    l_has_cl(2) = has_cl(2)==1
    l_has_cl(3) = has_cl(4)==1

END SUBROUTINE

SUBROUTINE   camspec_extra_init_v3(ipre_marged,like_file,l_like_file,sz143_file,l_sz143_file,tszxcib_file,l_tszxcib_file,ksz_file,l_ksz_file,beam_file,l_beam_file,data_vector,l_data_vector,&
                                   l_cib217_file,cib217_file,l_dust100_file,dust100_file,l_dust143_file,dust143_file,l_dust217_file,dust217_file,l_dust143x217_file,dust143x217_file,&
                                   icamspec_fiducial_foregrounds,l_camspec_fiducial_foregrounds,icamspec_fiducial_cl,l_camspec_fiducial_cl,lmins,lmaxs,spec_flag,icamspec_beam_mcmc_num,xdim,iclik_lmin,iclik_lmax,has_cl,bs_factor,sz_prior) 
    use CAMSPEC_EXTRA
    use temp_like_camspec3
    implicit none
    integer,intent(in)::l_like_file,l_beam_file,l_tszxcib_file,l_sz143_file,l_data_vector,l_ksz_file,xdim,ipre_marged,icamspec_beam_mcmc_num,iclik_lmin,iclik_lmax,l_camspec_fiducial_foregrounds,l_camspec_fiducial_cl
    integer,intent(in)::l_cib217_file,l_dust100_file,l_dust143_file,l_dust217_file,l_dust143x217_file,sz_prior
    
    character(len=l_like_file)::like_file
    character(len=l_sz143_file)::sz143_file
    character(len=l_tszxcib_file)::tszxcib_file
    character(len=l_ksz_file)::ksz_file
    character(len=l_beam_file)::beam_file
    character(len=l_data_vector)::data_vector
    character(len=l_camspec_fiducial_foregrounds)::icamspec_fiducial_foregrounds
    character(len=l_camspec_fiducial_cl)::icamspec_fiducial_cl
    character(len=l_cib217_file)::cib217_file
    character(len=l_dust217_file)::dust217_file
    character(len=l_dust100_file)::dust100_file
    character(len=l_dust143_file)::dust143_file
    character(len=l_dust143x217_file)::dust143x217_file

    integer,intent(in),dimension(6)::lmaxs,lmins,spec_flag,has_cl
    logical::pre_marged
    integer::i
    real*8::bs_factor
    
    n_nuisance = xdim
    do i=1,6
        want_spec(i) = spec_flag(i)==1
        camspec_lmins(i) = lmins(i)
        camspec_lmaxs(i) = lmaxs(i)
    enddo

    camspec_beam_mcmc_num = icamspec_beam_mcmc_num
    
    camspec_fiducial_foregrounds = icamspec_fiducial_foregrounds
    camspec_fiducial_cl = icamspec_fiducial_cl

    pre_marged = ipre_marged==1

    apply_tight_sz_prior = (sz_prior == 1)
    call like_init(pre_marged,like_file, sz143_file, tszxcib_file, ksz_file, beam_file,data_vector,cib217_file,dust100_file,dust143_file,dust217_file,dust143x217_file)

    cam_version=3

    clik_lmin = iclik_lmin
    clik_lmax = iclik_lmax
    lmin = minval(lminX(1:4))
    lmax = maxval(lmaxX(1:4))
    lminte = lminX(5)
    lmaxte = lmaxX(5)
    lminee = lminX(6)
    lmaxee = lmaxX(6)

    ALLOCATE(cltt(0:lmax+1))
    ALLOCATE(clte(0:lmaxte+1))
    ALLOCATE(clee(0:lmaxee+1))

    allocate(nuisance(xdim))
    
    npar = (clik_lmax+1-clik_lmin)*has_cl(1) + (clik_lmax+1-clik_lmin)*has_cl(2) + (clik_lmax+1-clik_lmin)*has_cl(4) + xdim
    
    beam_factor = bs_factor
    
    l_has_cl(1) = has_cl(1)==1
    l_has_cl(2) = has_cl(2)==1
    l_has_cl(3) = has_cl(4)==1
    
END SUBROUTINE

SUBROUTINE CAMSPEC_EXTRA_LKL_V2(LKL,CL)
    USE CAMSPEC_EXTRA
    use temp_like_camspec
    implicit none

    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:npar-1)::CL
    real(8)::zlike
    INTEGER::l,i,j,cnt,offset
    real(8)::tlkl

    offset = 0

    if (l_has_cl(1)) then
        cltt = 0
        DO l=lmin,lmax
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            cltt(l)=CL(l-clik_lmin)/2./3.141592653589793
        ENDDO
        offset = (clik_lmax+1-clik_lmin)
    endif

    if (l_has_cl(2)) then
        clee = 0
        DO l=lminee,lmaxee
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            clee(l)=CL(offset+l-clik_lmin)/2./3.141592653589793
        ENDDO
        offset = offset + (clik_lmax+1-clik_lmin)
    endif

    if (l_has_cl(3)) then
        clte = 0
        DO l=lminte,lmaxte
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            clte(l)=CL(offset+l-clik_lmin)/2./3.141592653589793
        ENDDO
        offset = offset + (clik_lmax+1-clik_lmin)
    endif
    
    do i=1,n_nuisance
        nuisance(i) = CL(offset + i-1)
    
    enddo
    
    call calc_like(tlkl,  cltt,clte,clee,nuisance)
    ! lkl is -2loglike clik returns loglik
    !print *,tlkl
    lkl = -tlkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL_V2

SUBROUTINE CAMSPEC_EXTRA_LKL_V3(LKL,CL)
    USE CAMSPEC_EXTRA
    use temp_like_camspec3
    implicit none

    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:npar-1)::CL
    real(8)::zlike
    INTEGER::l,i,j,cnt,offset
    real(8)::tlkl

    offset = 0

    if (l_has_cl(1)) then
        cltt = 0
        DO l=lmin,lmax
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            cltt(l)=CL(l-clik_lmin)/2./3.141592653589793
            !print *,l,cltt(l),CL(l-clik_lmin)
        ENDDo
        offset = (clik_lmax+1-clik_lmin)
    endif

    if (l_has_cl(2)) then
        clee = 0
        DO l=lminee,lmaxee
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            clee(l)=CL(offset+l-clik_lmin)/2./3.141592653589793
        ENDDO
        offset = offset + (clik_lmax+1-clik_lmin)
    endif

    if (l_has_cl(3)) then
        clte = 0
        DO l=lminte,lmaxte
            ! camspec expects cl/2pi !!! argl !
            !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
            clte(l)=CL(offset+l-clik_lmin)/2./3.141592653589793
        ENDDO
        offset = offset + (clik_lmax+1-clik_lmin)
    endif
    
    do i=1,n_nuisance
        nuisance(i) = CL(offset + i-1)
    
    enddo
    
    call calc_like(tlkl,  cltt,clte,clee,nuisance)
    ! lkl is -2loglike clik returns loglik
    !print *,tlkl
    lkl = -tlkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL_V3

SUBROUTINE CAMSPEC_EXTRA_FREE_V3()
    USE CAMSPEC_EXTRA
    use temp_like_camspec3

    DEALLOCATE(cltt)
    call camspec_clean()
    BOK =1
    ! forbid double instantiation of camspec. Need more work to cleanup correctly it seems...
END SUBROUTINE  CAMSPEC_EXTRA_FREE_v3

#endif

SUBROUTINE CAMSPEC_EXTRA_GETCASE(xc)
    USE CAMSPEC_EXTRA
    INTEGER::xc

    xc = xcase
END SUBROUTINE CAMSPEC_EXTRA_GETCASE

SUBROUTINE CAMSPEC_EXTRA_LKL(LKL,CL)
    USE CAMSPEC_EXTRA
    implicit none

    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:npar-1)::CL

    if (cam_version==1) then
        call CAMSPEC_EXTRA_LKL_V1(LKL,CL)
#ifndef CAMSPEC_V1
    else if (cam_version==2) then
        call CAMSPEC_EXTRA_LKL_V2(LKL,CL)
    else
        call CAMSPEC_EXTRA_LKL_V3(LKL,CL)
#endif
    endif

END SUBROUTINE CAMSPEC_EXTRA_LKL

SUBROUTINE CAMSPEC_EXTRA_LKL_V1(LKL,CL)
    USE CAMSPEC_EXTRA
    use temp_like
    implicit none

    REAL(8),INTENT(OUT)::LKL
    REAL(8),DIMENSION(0:npar-1)::CL
    real(8)::zlike, A_ps_100, A_ps_143, A_ps_217, A_cib_143, A_cib_217, A_sz, r_ps, r_cib, &
         cal0, cal1, cal2, xi, A_ksz
    real*8, dimension(1:beam_Nspec,1:num_modes_per_beam) :: beam_coeffs
    INTEGER::l,i,j,cnt,offset
    real(8)::tlkl

    cltt(:lmin-1) = 0
    DO l=lmin,lmax
        ! camspec expects cl/2pi !!! argl !
        !cltt(l)=CL(l-lmin)/2./3.14159265358979323846264338328
        cltt(l)=CL(l-lmin)/2./3.141592653589793
    ENDDO
    offset = (lmax+1-lmin)

    do i=1,num_non_beam
        nuisance(i) = CL(offset + i-1)
    enddo
    !A_ps_100  = CL(lmax+1-lmin + 0)
    !A_ps_143  = CL(lmax+1-lmin + 1)
    !A_ps_217  = CL(lmax+1-lmin + 2)
    !A_cib_143 = CL(lmax+1-lmin + 3)
    !A_cib_217 = CL(lmax+1-lmin + 4)
    !A_sz      = CL(lmax+1-lmin + 5)
    !r_ps      = CL(lmax+1-lmin + 6)
    !r_cib     = CL(lmax+1-lmin + 7) 
    !cal0      = CL(lmax+1-lmin + 8) 
    !cal1      = CL(lmax+1-lmin + 9) 
    !cal2      = CL(lmax+1-lmin + 10)     
    !xi        = CL(lmax+1-lmin + 11)     
    !A_ksz     = CL(lmax+1-lmin + 12)

    !print *,CL(lmax+1-lmin+14)

    cnt = 1
    DO i=1,keep_num
            nuisance(cnt+num_non_beam) = CL(offset+num_non_beam+cnt-1)
            cnt = cnt + 1
    ENDDO


    call calc_like(tlkl,  cltt,nuisance)
    ! lkl is -2loglike clik returns loglik
    !print *,tlkl
    lkl = -tlkl/2.

END SUBROUTINE CAMSPEC_EXTRA_LKL_V1


SUBROUTINE CAMSPEC_EXTRA_FG(rCL_FG,NUIS,lm)
    USE CAMSPEC_EXTRA
    use temp_like
    implicit none
    REAL(8),DIMENSION(1:num_non_beam+beam_Nspec*num_modes_per_beam)::NUIS
    REAL(8),dimension(4,0:lm)::rCL_FG
    integer::lm



    call COMPUTE_FG(rCL_FG,NUIS)
    
END SUBROUTINE


