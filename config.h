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

/* Define to 1 if you have the `acb_exp' function. */
#define HAVE_ACB_EXP 1

/* Define to 1 if you have the <acb.h> header file. */
#define HAVE_ACB_H 1

/* Define to 1 if you have the `acb_lgamma' function. */
#define HAVE_ACB_LGAMMA 1

/* Define to 1 if you have the `acb_set' function. */
#define HAVE_ACB_SET 1

/* Define to 1 if you have the `acb_sin' function. */
#define HAVE_ACB_SIN 1

/* Define to 1 if you have `alloca', as a function or macro. */
#define HAVE_ALLOCA 1

/* Define to 1 if you have <alloca.h> and it should be used (not on Ultrix).
   */
#define HAVE_ALLOCA_H 1

/* Define to 1 if you have the <arb.h> header file. */
#define HAVE_ARB_H 1

/* Define to 1 if you have the `arb_set' function. */
#define HAVE_ARB_SET 1

/* Define to 1 if you have the `backtrace' function. */
#define HAVE_BACKTRACE 1

/* Define to 1 if you have the `backtrace_symbols' function. */
#define HAVE_BACKTRACE_SYMBOLS 1

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define to 1 if you have the <cblas.h> header file. */
/* #undef HAVE_CBLAS_H */

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

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the `erf' function. */
#define HAVE_ERF 1

/* Define to 1 if you have the <execinfo.h> header file. */
#define HAVE_EXECINFO_H 1

/* Define to 1 if you have the `exp10' function. */
/* #undef HAVE_EXP10 */

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

/* Define to 1 if you have the `fork' function. */
#define HAVE_FORK 1

/* Whether the compiler supports GCC vector extensions */
#define HAVE_GCC_VECTOR_EXTENSIONS /**/

/* gsl 2.0 */
#define HAVE_GSL_2_0 1

/* gsl 2.2 */
#define HAVE_GSL_2_2 1

/* use GSL's BLAS */
/* #undef HAVE_GSL_CBLAS_H */

/* gsl support fixed gauss legendre rules */
#define HAVE_GSL_GLF 1

/* gsl support odeiv2 ode algorightms */
#define HAVE_GSL_ODEIV2 1

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

/* Define to 1 if you have the `arb' library (-larb). */
#define HAVE_LIBARB 1

/* Define to 1 if you have the `cuba' library (-lcuba). */
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

/* Define to 1 if you have the <mkl_blas.h> header file. */
/* #undef HAVE_MKL_BLAS_H */

/* Define to 1 if you have the <mkl_cblas.h> header file. */
/* #undef HAVE_MKL_CBLAS_H */

/* have mkl_lapacke.h header */
/* #undef HAVE_MKL_LAPACKE_H */

/* Define to 1 if you have the <mkl_lapack.h> header file. */
/* #undef HAVE_MKL_LAPACK_H */

/* Define to 1 if you have the `MKL_Set_Num_Threads' function. */
/* #undef HAVE_MKL_SET_NUM_THREADS */

/* Have the NLopt package */
#define HAVE_NLOPT 1

/* If nlopt version is bigger than 2.2 */
#define HAVE_NLOPT_2_2 1

/* Define to 1 if you have the `openblas_set_num_threads' function. */
/* #undef HAVE_OPENBLAS_SET_NUM_THREADS */

/* Define to 1 if you have the `powl' function. */
#define HAVE_POWL 1

/* Define to 1 if you have the `sin' function. */
#define HAVE_SIN 1

/* Define to 1 if you have the `sincos' function. */
/* #undef HAVE_SINCOS */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* SUNDIALS presence */
#define HAVE_SUNDIALS 1

/* SUNDIALS version */
/* #undef HAVE_SUNDIALS_2_5_0 */

/* SUNDIALS version */
/* #undef HAVE_SUNDIALS_2_6_0 */

/* SUNDIALS version */
#define HAVE_SUNDIALS_2_7_0 1

/* SUNDIALS has ARKode */
#define HAVE_SUNDIALS_ARKODE 1

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

/* Define to 1 if `fork' works. */
#define HAVE_WORKING_FORK 1

/* Define to 1 if `vfork' works. */
#define HAVE_WORKING_VFORK 1

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
#define NUMCOSMO_BINARY_AGE 1303

/* Define to the NUMCOSMO interface age */
#define NUMCOSMO_INTERFACE_AGE 0

/* Define to the NUMCOSMO major version */
#define NUMCOSMO_MAJOR_VERSION 0

/* Define to the NUMCOSMO micro version */
#define NUMCOSMO_MICRO_VERSION 3

/* Define to the NUMCOSMO minor version */
#define NUMCOSMO_MINOR_VERSION 13

/* Name of package */
#define PACKAGE "numcosmo"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://savannah.nongnu.org/bugs/?func=additem&group=numcosmo"

/* PACKAGE_DATA_DIR */
#define PACKAGE_DATA_DIR "/usr/local/share/numcosmo-0.13.3"

/* Define to the full name of this package. */
#define PACKAGE_NAME "numcosmo"

/* PACKAGE_SOURCE_DIR */
#define PACKAGE_SOURCE_DIR "/Users/Cosmology/Physics/Library/NumCosmo"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "numcosmo 0.13.3"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "numcosmo"

/* Define to the home page for this package. */
#define PACKAGE_URL "http://www.nongnu.org/numcosmo/"

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.13.3"

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
#define VERSION "0.13.3"

/* Define to 1 if on MINIX. */
/* #undef _MINIX */

/* gint */
#define _NCM_SUNDIALS_INT_TYPE glong

/* Define to 2 if the system does not provide POSIX.1 features except with
   this defined. */
/* #undef _POSIX_1_SOURCE */

/* Define to 1 if you need to in order for `stat' and other things to work. */
/* #undef _POSIX_SOURCE */

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

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* Define to `int' if <sys/types.h> does not define. */
/* #undef ssize_t */

/* Define as `fork' if `vfork' does not work. */
/* #undef vfork */
