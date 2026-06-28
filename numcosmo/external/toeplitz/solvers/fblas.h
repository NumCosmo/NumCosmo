/* collection of external declarations of Fortran 77 blas routines */

/* including the needed types */
#include "ftypes.h"

#ifdef  __cplusplus
extern "C" {
#endif

extern void F77CALL (caxpy) (const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (ccopy) (const finteger*,const fsinglecomplex*,
      const finteger*,fsinglecomplex*,const finteger*);

#if (defined(USE_F2C) || defined(RETURN_COMPLEX_WITH_POINTER))
extern void F77CALL (cdotc) (fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*);

extern void F77CALL (cdotu) (fsinglecomplex *,const finteger*,
      const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*);
#else
extern fsinglecomplex F77CALL (cdotc) (const finteger*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,const finteger*);

extern fsinglecomplex F77CALL (cdotu) (const finteger*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,const finteger*);
#endif

extern void F77CALL (cgbmv) (const char*,const finteger*,const finteger*,
      const finteger*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (cgemm) (const char*,const char*,const finteger*,
      const finteger*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (cgemv) (const char*,const finteger*,const finteger*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (cgerc) (const finteger*,const finteger*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (cgeru) (const finteger*,const finteger*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (chbmv) (const char*,const finteger*,const finteger*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (chemm) (const char*,const char*,const finteger*,
      const finteger*,const fsinglecomplex*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (chemv) (const char*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (cher) (const char*,const finteger*,const fsingle*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (cher2) (const char*,const finteger*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (cher2k) (const char*,const char*,const finteger*,
      const finteger*,const fsinglecomplex*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,const finteger*,const fsingle*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (cherk) (const char*,const char*,const finteger*,
      const finteger*,const fsingle*,const fsinglecomplex*,const finteger*,
      const fsingle*,fsinglecomplex*,const finteger*);

extern void F77CALL (chpmv) (const char*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (chpr) (const char*,const finteger*,const fsingle*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*);

extern void F77CALL (chpr2) (const char*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      const finteger*,fsinglecomplex*);

extern void F77CALL (crotg) (fsinglecomplex*,fsinglecomplex*,fsingle*,
      fsinglecomplex*);

extern void F77CALL (cscal) (const finteger*,const fsinglecomplex*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (csscal) (const finteger*,const fsingle*,fsinglecomplex*,
      const finteger*);

extern void F77CALL (cswap) (const finteger*,fsinglecomplex*,const finteger*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (csymm) (const char*,const char*,const finteger*,
      const finteger*,fsinglecomplex*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,const finteger*,const fsinglecomplex*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (csyr2k) (const char*,const char*,const finteger*,
      const finteger*,const fsinglecomplex*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,const finteger*,
      const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (csyrk) (const char*,const char*,const finteger*,
      const finteger*,const fsinglecomplex*,const fsinglecomplex*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (ctbmv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fsinglecomplex*,const finteger*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (ctbsv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fsinglecomplex*,const finteger*,
      fsinglecomplex*,const finteger*);

extern void F77CALL (ctpmv) (const char*,const char*,const char*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (ctpsv) (const char*,const char*,const char*,
      const finteger*,const fsinglecomplex*,fsinglecomplex*,const finteger*);

extern void F77CALL (ctrmm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (ctrmv) (const char*,const char*,const char*,
      const finteger*,const fsinglecomplex*,const finteger*,fsinglecomplex*,
      const finteger*);

extern void F77CALL (ctrsm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fsinglecomplex*,
      const fsinglecomplex*,const finteger*,fsinglecomplex*,const finteger*);

extern void F77CALL (ctrsv) (const char*,const char*,const char*,
      const finteger*,const fsinglecomplex*,const finteger*,fsinglecomplex*,
      const finteger*);

extern fdouble F77CALL (dasum) (const finteger*,const fdouble*,
      const finteger*);

extern void F77CALL (daxpy) (const finteger*,const fdouble*,const fdouble*,
      const finteger*,fdouble*,const finteger*);

extern fdouble F77CALL (dcabs1) (const fdoublecomplex*);

extern void F77CALL (dcopy) (const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern fdouble F77CALL (ddot) (const finteger*,const fdouble*,const finteger*,
      const fdouble*,const finteger*);

extern void F77CALL (dgbmv) (const char*,const finteger*,const finteger*,
      const finteger*,const finteger*,const fdouble*,const fdouble*,
      const finteger*,const fdouble*,const finteger*,const fdouble*,fdouble*,
      const finteger*);

extern void F77CALL (dgemm) (const char*,const char*,const finteger*,
      const finteger*,const finteger*,const fdouble*,const fdouble*,
      const finteger*,const fdouble*,const finteger*,const fdouble*,fdouble*,
      const finteger*);

extern void F77CALL (dgemv) (const char*,const finteger*,const finteger*,
      const fdouble*,const fdouble*,const finteger*,const fdouble*,
      const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dger) (const finteger*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,const fdouble*,const finteger*,
      const fdouble*,const finteger*);

extern fdouble F77CALL (dnrm2) (const finteger*,const fdouble*,
      const finteger*);

extern void F77CALL (drot) (const finteger*,fdouble*,const finteger*,fdouble*,
      const finteger*,const fdouble*,const fdouble*);

extern void F77CALL (drotg) (fdouble*,fdouble*,fdouble*,fdouble*);

extern void F77CALL (dsbmv) (const char*,const finteger*,const finteger*,
      const fdouble*,const fdouble*,const finteger*,const fdouble*,
      const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dscal) (const finteger*,const fdouble*,fdouble*,
      const finteger*);

extern void F77CALL (dspmv) (const char*,const finteger*,const fdouble*,
      const fdouble*,const fdouble*,const finteger*,const fdouble*,
      fdouble*,const finteger*);

extern void F77CALL (dspr) (const char*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,fdouble*);

extern void F77CALL (dspr2) (const char*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,const fdouble*,const finteger*,fdouble*);

extern void F77CALL (dswap) (const finteger*,fdouble*,const finteger*,
      fdouble*,const finteger*);

extern void F77CALL (dsymm) (const char*,const char*,const finteger*,
      const finteger*,const fdouble*,const fdouble*,const finteger*,
      const fdouble*,const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dsymv) (const char*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,const fdouble*,const finteger*,
      const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dsyr) (const char*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,fdouble*,const finteger*);

extern void F77CALL (dsyr2) (const char*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern void F77CALL (dsyr2k) (const char*,const char*,const finteger*,
      const finteger*,const fdouble*,const fdouble*,const finteger*,
      const fdouble*,const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dsyrk) (const char*,const char*,const finteger*,
      const finteger*,const fdouble*,const fdouble*,const finteger*,
      const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dtbmv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern void F77CALL (dtbsv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern void F77CALL (dtpmv) (const char*,const char*,const char*,
      const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dtpsv) (const char*,const char*,const char*,
      const finteger*,const fdouble*,fdouble*,const finteger*);

extern void F77CALL (dtrmm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fdouble*,const fdouble*,
      const finteger*,fdouble*,const finteger*);

extern void F77CALL (dtrmv) (const char*,const char*,const char*,
      const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern void F77CALL (dtrsm) (const char*,const char*,const char*,
      const char*,const finteger*,const finteger*,const fdouble*,
      const fdouble*,const finteger*,fdouble*,const finteger*);

extern void F77CALL (dtrsv) (const char*,const char*,const char*,
      const finteger*,const fdouble*,const finteger*,
      fdouble*,const finteger*);

extern fdouble F77CALL (dzasum) (const finteger*,const fdoublecomplex*,
      const finteger*);

extern fdouble F77CALL (dznrm2) (const finteger*,const fdoublecomplex*,
      const finteger*);

extern finteger F77CALL (icamax) (const finteger*,const fsinglecomplex*,
      const finteger*);

extern finteger F77CALL (idamax) (const finteger*,const fdouble*,
      const finteger*);

extern finteger F77CALL (isamax) (const finteger*,const fsingle*,
      const finteger*);

extern finteger F77CALL (izamax) (const finteger*,const fdoublecomplex*,
      const finteger*);

extern flogical F77CALL (lsame) (const char*,const char*);

#ifdef USE_F2C
extern fdouble F77CALL (sasum) (const finteger*,const fsingle*,
      const finteger*);
#else
extern fsingle F77CALL (sasum) (const finteger*,const fsingle*,
      const finteger*);
#endif

extern void F77CALL (saxpy) (const finteger*,const fsingle*,const fsingle*,
      const finteger*,fsingle*,const finteger*);

#ifdef USE_F2C
extern fdouble F77CALL (scasum) (const finteger*,const fsinglecomplex*,
      const finteger*);

extern fdouble F77CALL (scnrm2) (const finteger*,const fsinglecomplex*,
      const finteger*);
#else
extern fsingle F77CALL (scasum) (const finteger*,const fsinglecomplex*,
      const finteger*);

extern fsingle F77CALL (scnrm2) (const finteger*,const fsinglecomplex*,
      const finteger*);
#endif

extern void F77CALL (scopy) (const finteger*,const fsingle*,const finteger*,
      fsingle*,const finteger*);

#ifdef USE_F2C
extern fdouble F77CALL (sdot) (const finteger*,const fsingle*,const finteger*,
      const fsingle*,const finteger*);
#else
extern fsingle F77CALL (sdot) (const finteger*,const fsingle*,const finteger*,
      const fsingle*,const finteger*);
#endif

extern void F77CALL (sgbmv) (const char*,const finteger*,const finteger*,
      const finteger*,const finteger*,const fsingle*,const fsingle*,
      const finteger*,const fsingle*,const finteger*,const fsingle*,fsingle*,
      const finteger*);

extern void F77CALL (sgemm) (const char*,const char*,const finteger*,
      const finteger*,const finteger*,const fsingle*,const fsingle*,
      const finteger*,const fsingle*,const finteger*,const fsingle*,fsingle*,
      const finteger*);

extern void F77CALL (sgemv) (const char*,const finteger*,const finteger*,
      const fsingle*,const fsingle*,const finteger*,const fsingle*,
      const finteger*,const fsingle*,fsingle*,const finteger*);

extern void F77CALL (sger) (const finteger*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,const fsingle*,const finteger*,
      fsingle*,const finteger*);

#ifdef USE_F2C
extern fdouble F77CALL (snrm2) (const finteger*,const fsingle*,
      const finteger*);
#else
extern fsingle F77CALL (snrm2) (const finteger*,const fsingle*,
      const finteger*);
#endif

extern void F77CALL (srot) (const finteger*,fsingle*,const finteger*,
      fsingle*,const finteger*,const fsingle*,const fsingle*);

extern void F77CALL (srotg) (fsingle*,fsingle*,fsingle*,fsingle*);

extern void F77CALL (ssbmv) (const char*,const finteger*,const finteger*,
      const fsingle*,const fsingle*,const finteger*,const fsingle*,
      const finteger*,const fsingle*,fsingle*,const finteger*);

extern void F77CALL (sscal) (const finteger*,const fsingle*,fsingle*,
      const finteger*);

extern void F77CALL (sspmv) (const char*,const finteger*,const fsingle*,
      const fsingle*,const fsingle*,const finteger*,const fsingle*,
      fsingle*,const finteger*);

extern void F77CALL (sspr) (const char*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,fsingle*);

extern void F77CALL (sspr2) (const char*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,const fsingle*,const finteger*,fsingle*);

extern void F77CALL (sswap) (const finteger*,fsingle*,const finteger*,
      fsingle*,const finteger*);

extern void F77CALL (ssymm) (const char*,const char*,const finteger*,
      const finteger*,const fsingle*,const fsingle*,const finteger*,
      const fsingle*,const finteger*,const fsingle*,fsingle*,const finteger*);

extern void F77CALL (ssymv) (const char*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,const fsingle*,const finteger*,
      const fsingle*,fsingle*,const finteger*);

extern void F77CALL (ssyr) (const char*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,fsingle*,const finteger*);

extern void F77CALL (ssyr2) (const char*,const finteger*,const fsingle*,
      const fsingle*,const finteger*,const fsingle*,const finteger*,
      fsingle*,const finteger*);

extern void F77CALL (ssyr2k) (const char*,const char*,const finteger*,const finteger*,fsingle*,
      fsingle*,const finteger*,fsingle*,const finteger*,fsingle*,fsingle*,const finteger*);

extern void F77CALL (ssyrk) (const char*,const char*,const finteger*,
      const finteger*,const fsingle*,const fsingle*,const finteger*,
      const fsingle*,fsingle*,const finteger*);

extern void F77CALL (stbmv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fsingle*,
      const finteger*,fsingle*,const finteger*);

extern void F77CALL (stbsv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fsingle*,const finteger*,
      fsingle*,const finteger*);

extern void F77CALL (stpmv) (const char*,const char*,const char*,
      const finteger*,const fsingle*,fsingle*,const finteger*);

extern void F77CALL (stpsv) (const char*,const char*,const char*,
      const finteger*,const fsingle*,fsingle*,const finteger*);

extern void F77CALL (strmm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fsingle*,const fsingle*,
      const finteger*,fsingle*,const finteger*);

extern void F77CALL (strmv) (const char*,const char*,const char*,
      const finteger*,const fsingle*,const finteger*,fsingle*,const finteger*);

extern void F77CALL (strsm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fsingle*,const fsingle*,
      const finteger*,fsingle*,const finteger*);

extern void F77CALL (strsv) (const char*,const char*,const char*,
      const finteger*,const fsingle*,const finteger*,fsingle*,const finteger*);

extern void F77CALL (xerbla) (const char*,const finteger*);

extern void F77CALL (zaxpy) (const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,const finteger*);

extern void F77CALL (zcopy) (const finteger*,const fdoublecomplex*,
      const finteger*,fdoublecomplex*,const finteger*);

#if (defined(USE_F2C) || defined(RETURN_COMPLEX_WITH_POINTER))
extern void F77CALL (zdotc) (fdoublecomplex *,const finteger*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*);

extern void F77CALL (zdotu) (fdoublecomplex *,const finteger*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*);
#else
extern fdoublecomplex F77CALL (zdotc) (const finteger*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*);

extern fdoublecomplex F77CALL (zdotu) (const finteger*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*);
#endif

extern void F77CALL (zdscal) (const finteger*,const fdouble*,fdoublecomplex*,
      const finteger*);

extern void F77CALL (zgbmv) (const char*,const finteger*,const finteger*,
      const finteger*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zgemm) (const char*,const char*,const finteger*,
      const finteger*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zgemv) (const char*,const finteger*,const finteger*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (zgerc) (const finteger*,const finteger*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,
      const finteger*);

extern void F77CALL (zgeru) (const finteger*,const finteger*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,
      const finteger*);

extern void F77CALL (zhbmv) (const char*,const finteger*,const finteger*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (zhemm) (const char*,const char*,const finteger*,
      const finteger*,const fdoublecomplex*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zhemv) (const char*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,
      const finteger*);

extern void F77CALL (zher) (const char*,const finteger*,const fdouble*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,const finteger*);

extern void F77CALL (zher2) (const char*,const finteger*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,const finteger*);

extern void F77CALL (zher2k) (const char*,const char*,const finteger*,
      const finteger*,const fdoublecomplex*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*,const fdouble*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (zherk) (const char*,const char*,const finteger*,
      const finteger*,const fdouble*,const fdoublecomplex*,const finteger*,
      const fdouble*,fdoublecomplex*,const finteger*);

extern void F77CALL (zhpmv) (const char*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zhpr) (const char*,const finteger*,const fdouble*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*);

extern void F77CALL (zhpr2) (const char*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,const fdoublecomplex*,
      const finteger*,fdoublecomplex*);

extern void F77CALL (zrotg) (fdoublecomplex*,fdoublecomplex*,fdouble*,
      fdoublecomplex*);

extern void F77CALL (zscal) (const finteger*,const fdoublecomplex*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (zswap) (const finteger*,fdoublecomplex*,const finteger*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (zsymm) (const char*,const char*,const finteger*,
      const finteger*,const fdoublecomplex*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zsyr2k) (const char*,const char*,const finteger*,
      const finteger*,const fdoublecomplex*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,const finteger*,
      const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (zsyrk) (const char*,const char*,const finteger*,
      const finteger*,const fdoublecomplex*,const fdoublecomplex*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (ztbmv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fdoublecomplex*,const finteger*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (ztbsv) (const char*,const char*,const char*,
      const finteger*,const finteger*,const fdoublecomplex*,const finteger*,
      fdoublecomplex*,const finteger*);

extern void F77CALL (ztpmv) (const char*,const char*,const char*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (ztpsv) (const char*,const char*,const char*,
      const finteger*,const fdoublecomplex*,fdoublecomplex*,const finteger*);

extern void F77CALL (ztrmm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,const finteger*);

extern void F77CALL (ztrmv) (const char*,const char*,const char*,
      const finteger*,const fdoublecomplex*,const finteger*,fdoublecomplex*,
      const finteger*);

extern void F77CALL (ztrsm) (const char*,const char*,const char*,const char*,
      const finteger*,const finteger*,const fdoublecomplex*,
      const fdoublecomplex*,const finteger*,fdoublecomplex*,const finteger*);

extern void F77CALL (ztrsv) (const char*,const char*,const char*,
      const finteger*,const fdoublecomplex*,const finteger*,fdoublecomplex*,
      const finteger*);

#ifdef  __cplusplus
}
#endif
