/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to one of `_getb67', `GETB67', `getb67' for Cray-2 and Cray-YMP
   systems. This function is required for `alloca.c' support on those systems.
   */
/* #undef CRAY_STACKSEG_END */

/* Define to 1 if using `alloca.c'. */
/* #undef C_ALLOCA */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define F77_FUNC(name,NAME) name ## _

/* As F77_FUNC, but for C identifiers containing underscores. */
#define F77_FUNC_(name,NAME) name ## _

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* optimize gsl access functions */
#define GSL_RANGE_CHECK_OFF /**/

/* Disable GObject cast checks */
#define G_DISABLE_CAST_CHECKS 1

/* force glib to disable inline */
/* #undef G_IMPLEMENT_INLINES */

/* Define to 1 to support Advanced Bit Manipulation */
/* #undef HAVE_ABM */

/* Define to 1 if you have the `acb_exp' function. */
/* #undef HAVE_ACB_EXP */

/* Define to 1 if you have the <acb.h> header file. */
/* #undef HAVE_ACB_H */

/* Define to 1 if you have the `acb_lgamma' function. */
/* #undef HAVE_ACB_LGAMMA */

/* Define to 1 if you have the `acb_set' function. */
/* #undef HAVE_ACB_SET */

/* Define to 1 if you have the `acb_sin' function. */
/* #undef HAVE_ACB_SIN */

/* Define to 1 to support Multi-Precision Add-Carry Instruction Extensions */
/* #undef HAVE_ADX */

/* Define to 1 to support Advanced Encryption Standard New Instruction Set
   (AES-NI) */
/* #undef HAVE_AES */

/* Define to 1 if you have `alloca', as a function or macro. */
#define HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#define HAVE_ALLOCA_H 1

/* Support Altivec instructions */
/* #undef HAVE_ALTIVEC */

/* Define to 1 if you have the <arb.h> header file. */
/* #undef HAVE_ARB_H */

/* Define to 1 if you have the `arb_set' function. */
/* #undef HAVE_ARB_SET */

/* Define to 1 to support Advanced Vector Extensions */
/* #undef HAVE_AVX */

/* Define to 1 to support Advanced Vector Extensions 2 */
/* #undef HAVE_AVX2 */

/* Define to 1 to support AVX-512 Byte and Word Instructions */
/* #undef HAVE_AVX512_BW */

/* Define to 1 to support AVX-512 Conflict Detection Instructions */
/* #undef HAVE_AVX512_CD */

/* Define to 1 to support AVX-512 Doubleword and Quadword Instructions */
/* #undef HAVE_AVX512_DQ */

/* Define to 1 to support AVX-512 Exponential & Reciprocal Instructions */
/* #undef HAVE_AVX512_ER */

/* Define to 1 to support AVX-512 Foundation Extensions */
/* #undef HAVE_AVX512_F */

/* Define to 1 to support AVX-512 Integer Fused Multiply Add Instructions */
/* #undef HAVE_AVX512_IFMA */

/* Define to 1 to support AVX-512 Conflict Prefetch Instructions */
/* #undef HAVE_AVX512_PF */

/* Define to 1 to support AVX-512 Vector Byte Manipulation Instructions */
/* #undef HAVE_AVX512_VBMI */

/* Define to 1 to support AVX-512 Vector Length Extensions */
/* #undef HAVE_AVX512_VL */

/* Define to 1 if you have the `backtrace' function. */
#define HAVE_BACKTRACE 1

/* Define to 1 if you have the `backtrace_symbols' function. */
#define HAVE_BACKTRACE_SYMBOLS 1

/* use GSL's BLAS */
#define HAVE_BLAS 1

/* Define to 1 to support Bit Manipulation Instruction Set 1 */
/* #undef HAVE_BMI1 */

/* Define to 1 to support Bit Manipulation Instruction Set 2 */
/* #undef HAVE_BMI2 */

/* Define to 1 if you have the <cblas.h> header file. */
#define HAVE_CBLAS_H 1

/* Define to 1 if the system has the type `CBLAS_ORDER'. */
#define HAVE_CBLAS_ORDER 1

/* Have the cfitsio package */
#define HAVE_CFITSIO 1

/* have clapack support */
/* #undef HAVE_CLAPACK */

/* have clapack.h header */
/* #undef HAVE_CLAPACK_H */

/* Define to 1 if you have the `cos' function. */
#define HAVE_COS 1

/* Define to 1 if you have the declaration of `isfinite', and to 0 if you
   don't. */
#define HAVE_DECL_ISFINITE 1

/* Define to 1 if you have the declaration of `lgamma_r', and to 0 if you
   don't. */
#define HAVE_DECL_LGAMMA_R 1

/* Define to 1 if you have the `dgeevx_' function. */
#define HAVE_DGEEVX_ 1

/* Define to 1 if you have the `dgeev_' function. */
#define HAVE_DGEEV_ 1

/* Define to 1 if you have the `dgelqf_' function. */
#define HAVE_DGELQF_ 1

/* Define to 1 if you have the `dgels_' function. */
#define HAVE_DGELS_ 1

/* Define to 1 if you have the `dgeqlf_' function. */
#define HAVE_DGEQLF_ 1

/* Define to 1 if you have the `dgeqrf_' function. */
#define HAVE_DGEQRF_ 1

/* Define to 1 if you have the `dgerqf_' function. */
#define HAVE_DGERQF_ 1

/* Define to 1 if you have the `dggglm_' function. */
#define HAVE_DGGGLM_ 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `dposv_' function. */
#define HAVE_DPOSV_ 1

/* Define to 1 if you have the `dpotrf_' function. */
#define HAVE_DPOTRF_ 1

/* Define to 1 if you have the `dpotri_' function. */
#define HAVE_DPOTRI_ 1

/* Define to 1 if you have the `dpotrs_' function. */
#define HAVE_DPOTRS_ 1

/* Define to 1 if you have the `dptsv_' function. */
#define HAVE_DPTSV_ 1

/* Define to 1 if you have the `dsyevd_' function. */
#define HAVE_DSYEVD_ 1

/* Define to 1 if you have the `dsyevr_' function. */
#define HAVE_DSYEVR_ 1

/* Define to 1 if you have the `dsysvxx_' function. */
/* #undef HAVE_DSYSVXX_ */

/* Define to 1 if you have the `dsysvx_' function. */
#define HAVE_DSYSVX_ 1

/* Define to 1 if you have the `dsysv_' function. */
#define HAVE_DSYSV_ 1

/* Define to 1 if you have the `dsytrf_' function. */
#define HAVE_DSYTRF_ 1

/* Define to 1 if you have the `dsytri_' function. */
#define HAVE_DSYTRI_ 1

/* Define to 1 if you have the `dsytrs_' function. */
#define HAVE_DSYTRS_ 1

/* Define to 1 if you have the `dtrsv_' function. */
#define HAVE_DTRSV_ 1

/* Define to 1 if you have the `erf' function. */
#define HAVE_ERF 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the `exp10' function. */
#define HAVE_EXP10 1

/* Have the fast fourier package */
#define HAVE_FFTW3 1

/* Have the fast fourier package (float) */
#define HAVE_FFTW3F 1

/* fftw has alloc functions */
#define HAVE_FFTW3_ALLOC 1

/* Define to 1 if you have the `finite' function. */
#define HAVE_FINITE 1

/* Define to 1 if you have the `fma' function. */
#define HAVE_FMA 1

/* Define to 1 to support Fused Multiply-Add Extensions 3 */
/* #undef HAVE_FMA3 */

/* Define to 1 to support Fused Multiply-Add Extensions 4 */
/* #undef HAVE_FMA4 */

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Whether the compiler supports GCC vector extensions */
#define HAVE_GCC_VECTOR_EXTENSIONS /**/

/* gsl 2.2 */
#define HAVE_GSL_2_2 1

/* gsl 2.4 */
#define HAVE_GSL_2_4 1

/* use GSL's BLAS header */
/* #undef HAVE_GSL_CBLAS_H */

/* Defined if you have HDF5 support */
/* #undef HAVE_HDF5 */

/* use inline functions in GSL */
#define HAVE_INLINE 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* have lapacke support */
/* #undef HAVE_LAPACKE */

/* have lapacke.h header */
/* #undef HAVE_LAPACKE_H */

/* Define to 1 if you have the `lgamma_r' function. */
#define HAVE_LGAMMA_R 1

/* Define to 1 if you have the `arb' library (-larb). */
/* #undef HAVE_LIBARB */

/* Have any version of libcuba */
#define HAVE_LIBCUBA 1

/* Have the libcuba >= 3.1 */
/* #undef HAVE_LIBCUBA_3_1 */

/* Have the libcuba >= 3.3 */
/* #undef HAVE_LIBCUBA_3_3 */

/* Have the libcuba >= 4.0 */
#define HAVE_LIBCUBA_4_0 1

/* Define to 1 if you have the `gmp' library (-lgmp). */
#define HAVE_LIBGMP 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `mpfr' library (-lmpfr). */
#define HAVE_LIBMPFR 1

/* Define to 1 if the type `long double' works and has more range or precision
   than `double'. */
#define HAVE_LONG_DOUBLE 1

/* Define to 1 if the type `long double' works and has more range or precision
   than `double'. */
#define HAVE_LONG_DOUBLE_WIDER 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <mkl_cblas.h> header file. */
/* #undef HAVE_MKL_CBLAS_H */

/* have mkl_lapacke.h header */
/* #undef HAVE_MKL_LAPACKE_H */

/* Define to 1 if you have the <mkl_lapack.h> header file. */
/* #undef HAVE_MKL_LAPACK_H */

/* Define to 1 if you have the `MKL_Set_Num_Threads' function. */
/* #undef HAVE_MKL_SET_NUM_THREADS */

/* Define to 1 to support Multimedia Extensions */
/* #undef HAVE_MMX */

/* Define if you have the MPI library. */
/* #undef HAVE_MPI */

/* Define to 1 to support Memory Protection Extensions */
/* #undef HAVE_MPX */

/* Have the NLopt package */
#define HAVE_NLOPT 1

/* If nlopt version is bigger than 2.2 */
#define HAVE_NLOPT_2_2 1

/* Define to 1 if you have the `openblas_set_num_threads' function. */
#define HAVE_OPENBLAS_SET_NUM_THREADS 1

/* Define to 1 if you have the `powl' function. */
#define HAVE_POWL 1

/* Define to 1 to support Prefetch Vector Data Into Caches WT1 */
/* #undef HAVE_PREFETCHWT1 */

/* Define to 1 to support Digital Random Number Generator */
/* #undef HAVE_RDRND */

/* Define to 1 to support Secure Hash Algorithm Extension */
/* #undef HAVE_SHA */

/* Define to 1 if you have the `sin' function. */
#define HAVE_SIN 1

/* Define to 1 if you have the `sincos' function. */
#define HAVE_SINCOS 1

/* Define to 1 to support Streaming SIMD Extensions */
/* #undef HAVE_SSE */

/* Define to 1 to support Streaming SIMD Extensions */
/* #undef HAVE_SSE2 */

/* Define to 1 to support Streaming SIMD Extensions 3 */
/* #undef HAVE_SSE3 */

/* Define to 1 to support Streaming SIMD Extensions 4.1 */
/* #undef HAVE_SSE4_1 */

/* Define to 1 to support Streaming SIMD Extensions 4.2 */
/* #undef HAVE_SSE4_2 */

/* Define to 1 to support AMD Streaming SIMD Extensions 4a */
/* #undef HAVE_SSE4a */

/* Define to 1 to support Supplemental Streaming SIMD Extensions 3 */
/* #undef HAVE_SSSE3 */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `vfork' function. */
#define HAVE_VFORK 1

/* Define to 1 if you have the <vfork.h> header file. */
/* #undef HAVE_VFORK_H */

/* Support VSX instructions */
/* #undef HAVE_VSX */

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

/* Define to 1 to support eXtended Operations Extensions */
/* #undef HAVE_XOP */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Number of cores of the current machine */
#define NCM_NCORES 4

/* Number of threads in the pool */
#define NCM_THREAD_POOL_MAX 4

/* Maximum number of dimensions */
/* #undef NCOMP */

/* Maximum number of components */
/* #undef NDIM */

/* libcuba C interface */
#define NOUNDERSCORE 1

/* Define to the NUMCOSMO binary age */
#define NUMCOSMO_BINARY_AGE 1504

/* Define to the NUMCOSMO interface age */
#define NUMCOSMO_INTERFACE_AGE 0

/* Define to the NUMCOSMO major version */
#define NUMCOSMO_MAJOR_VERSION 0

/* Define to the NUMCOSMO micro version */
#define NUMCOSMO_MICRO_VERSION 4

/* Define to the NUMCOSMO minor version */
#define NUMCOSMO_MINOR_VERSION 15

/* Name of package */
#define PACKAGE "numcosmo"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://savannah.nongnu.org/bugs/?func=additem&group=numcosmo"

/* PACKAGE_DATA_DIR */
#define PACKAGE_DATA_DIR "/usr/local/share/numcosmo-0.15.4"

/* Define to the full name of this package. */
#define PACKAGE_NAME "numcosmo"

/* PACKAGE_SOURCE_DIR */
#define PACKAGE_SOURCE_DIR "/home/cinthia/NumCosmo"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "numcosmo 0.15.4"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "numcosmo"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.nongnu.org/numcosmo/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.15.4"

/* If using the C implementation of alloca, define if you know the
   direction of stack growth for your system; otherwise it will be
   automatically deduced at runtime.
	STACK_DIRECTION > 0 => grows toward higher addresses
	STACK_DIRECTION < 0 => grows toward lower addresses
	STACK_DIRECTION = 0 => direction of growth unknown */
/* #undef STACK_DIRECTION */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Enable extensions on AIX 3, Interix.  */
#ifndef _ALL_SOURCE
# define _ALL_SOURCE 1
#endif
/* Enable GNU extensions on systems that have them.  */
#ifndef _GNU_SOURCE
# define _GNU_SOURCE 1
#endif
/* Enable threading extensions on Solaris.  */
#ifndef _POSIX_PTHREAD_SEMANTICS
# define _POSIX_PTHREAD_SEMANTICS 1
#endif
/* Enable extensions on HP NonStop.  */
#ifndef _TANDEM_SOURCE
# define _TANDEM_SOURCE 1
#endif
/* Enable general extensions on Solaris.  */
#ifndef __EXTENSIONS__
# define __EXTENSIONS__ 1
#endif


/* Version number of package */
#define VERSION "0.15.4"

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* fffree */
/* #undef fffree */

/* fits_free_memory */
/* #undef fits_free_memory */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to `int' if <sys/types.h> does not define. */
/* #undef pid_t */

/* Define to the equivalent of the C99 'restrict' keyword, or to
   nothing if this is not supported.  Do not define if restrict is
   supported directly.  */
#define restrict __restrict
/* Work around a bug in Sun C++: it does not support _Restrict or
   __restrict__, even though the corresponding Sun C compiler ends up with
   "#define restrict _Restrict" or "#define restrict __restrict__" in the
   previous line.  Perhaps some future version of Sun C++ will work with
   restrict; if so, hopefully it defines __RESTRICT like Sun C does.  */
#if defined __SUNPRO_CC && !defined __RESTRICT
# define _Restrict
# define __restrict__
#endif

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */

/* Define to empty if the keyword `volatile' does not work. Warning: valid
   code using `volatile' can become incorrect without. Disable with care. */
/* #undef volatile */

#include <numcosmo/config_extra.h>
