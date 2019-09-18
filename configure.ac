dnl ***************************************************************************
dnl Process this file with autoconf to produce a configure script.
dnl ***************************************************************************

#
# The following version number definitions apply to GLib, GModule, GObject,
# GThread and GIO as a whole, so if changes occurred in any of them, they are
# all treated with the same interface and binary age.
#
# Making releases:
#   glib_micro_version += 1;
#   glib_interface_age += 1;
#   glib_binary_age += 1;
# if any functions have been added, set glib_interface_age to 0.
# if backwards compatibility has been broken,
# set glib_binary_age _and_ glib_interface_age to 0.
#
# remember to add a GLIB_VERSION_2_xx macro every time the minor version is
# bumped, as well as the GLIB_DEPRECATED_IN and GLIB_AVAILABLE_IN macros
# for that version - see gversion.h for further information.
#
# in easier to understand terms:
#
# <mclasen> on the stable branch, interface age == micro
# <mclasen> on the unstable (ie master), interface age = 0

m4_define([numcosmo_major_version], [0])
m4_define([numcosmo_minor_version], [14])
m4_define([numcosmo_micro_version], [2])
m4_define([numcosmo_interface_age], [0])
m4_define([numcosmo_binary_age],
          [m4_eval(100 * numcosmo_minor_version + numcosmo_micro_version)])
m4_define([numcosmo_version],
          [numcosmo_major_version.numcosmo_minor_version.numcosmo_micro_version])

# libtool version related macros
m4_define([numcosmo_lt_release], [numcosmo_major_version.numcosmo_minor_version])
m4_define([numcosmo_lt_current],
          [m4_eval(100 * numcosmo_minor_version + numcosmo_micro_version - numcosmo_interface_age)])
m4_define([numcosmo_lt_revision], [numcosmo_interface_age])
m4_define([numcosmo_lt_age], [m4_eval(numcosmo_binary_age - numcosmo_interface_age)])
m4_define([numcosmo_lt_current_minus_age],
          [m4_eval(numcosmo_lt_current - numcosmo_lt_age)])                    
          
AC_INIT([numcosmo],[numcosmo_version],[https://savannah.nongnu.org/bugs/?func=additem&group=numcosmo],[numcosmo],[http://www.nongnu.org/numcosmo/])
AC_CONFIG_SRCDIR([numcosmo/numcosmo.h])

dnl ***************************************************************************
dnl add the correct m4 directory location.
dnl ***************************************************************************
AC_CONFIG_MACRO_DIR([m4])

AC_PREREQ([2.69])
AM_INIT_AUTOMAKE([1.11 gnu -Wno-portability subdir-objects])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

PACKAGE=numcosmo
NUMCOSMO_MAJOR_VERSION=numcosmo_major_version
NUMCOSMO_MINOR_VERSION=numcosmo_minor_version
NUMCOSMO_MICRO_VERSION=numcosmo_micro_version
NUMCOSMO_INTERFACE_AGE=numcosmo_interface_age
NUMCOSMO_BINARY_AGE=numcosmo_binary_age
NUMCOSMO_VERSION=numcosmo_version

AC_SUBST(NUMCOSMO_MAJOR_VERSION)
AC_SUBST(NUMCOSMO_MINOR_VERSION)
AC_SUBST(NUMCOSMO_MICRO_VERSION)
AC_SUBST(NUMCOSMO_VERSION)
AC_SUBST(NUMCOSMO_INTERFACE_AGE)
AC_SUBST(NUMCOSMO_BINARY_AGE)

AC_DEFINE(NUMCOSMO_MAJOR_VERSION, [numcosmo_major_version],
          [Define to the NUMCOSMO major version])
AC_DEFINE(NUMCOSMO_MINOR_VERSION, [numcosmo_minor_version],
          [Define to the NUMCOSMO minor version])
AC_DEFINE(NUMCOSMO_MICRO_VERSION, [numcosmo_micro_version],
          [Define to the NUMCOSMO micro version])
AC_DEFINE(NUMCOSMO_INTERFACE_AGE, [numcosmo_interface_age],
          [Define to the NUMCOSMO interface age])
AC_DEFINE(NUMCOSMO_BINARY_AGE, [numcosmo_binary_age],
          [Define to the NUMCOSMO binary age])

# libtool versioning
LT_RELEASE=numcosmo_lt_release
LT_CURRENT=numcosmo_lt_current
LT_REVISION=numcosmo_lt_revision
LT_AGE=numcosmo_lt_age
LT_CURRENT_MINUS_AGE=numcosmo_lt_current_minus_age
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)
AC_SUBST(LT_CURRENT_MINUS_AGE)

AC_CONFIG_HEADERS([config.h])
AH_BOTTOM([#include <numcosmo/config_extra.h>])

AM_MAINTAINER_MODE

AX_CODE_COVERAGE
  
dnl ***************************************************************************
dnl things required by automake
dnl AC_ARG_PROGRAM
dnl ***************************************************************************
AC_PROG_MAKE_SET

dnl ***************************************************************************
dnl Check for which system.
dnl ***************************************************************************
AC_CANONICAL_HOST

dnl ***************************************************************************
dnl Checking from mpi support 
dnl ***************************************************************************
dnl If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
dnl If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
dnl If --with-mpi=no is used, use a standard C compiler instead.

AC_ARG_WITH([mpi], [AS_HELP_STRING([--with-mpi],
  [compile with MPI (parallelization) support. If none is found,
  MPI is not used. Default: no])
  ],,[with_mpi=no])

if test x$with_mpi != "xno"; then
  AC_CHECK_TOOLS([CC], [mpicc mpixlc_r mpixlc hcc mpxlc_r mpxlc sxmpicc mpifcc mpgcc mpcc cmpicc cc gcc])
  AC_CHECK_TOOLS([FC], [mpif95 mpxlf95_r mpxlf95 ftn mpif90 mpxlf90_r mpxlf90 mpf90 cmpif90c sxmpif90 mpif77 hf77 mpxlf_r mpxlf mpifrt mpf77 cmpifc xlf95 pgf95 pathf95 ifort g95 f95 fort ifc efc openf95 sunf95 crayftn gfortran lf95 ftn xlf90 f90 pgf90 pghpf pathf90 epcf90 sxf90 openf90 sunf90 xlf f77 frt pgf77 pathf77 g77 cf77 fort77 fl32 af77])
fi

dnl ***************************************************************************
dnl Enable gnu extensions when available
dnl ***************************************************************************
AC_USE_SYSTEM_EXTENSIONS

dnl ***************************************************************************
dnl Compilers
dnl ***************************************************************************
AC_LANG(C)

AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])

AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
  use_mpi=no
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])

if test x$use_mpi != "xno"; then
  AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])
  have_mpi_support="#define NUMCOSMO_HAVE_MPI 1"
  AC_SUBST(have_mpi_support)
fi

AC_PROG_CC_C99
AC_PROG_F77

dnl ***************************************************************************
dnl Checks for programs.
dnl ***************************************************************************
AC_PROG_INSTALL   
AC_PROG_LN_S
AC_PROG_CC_STDC
AC_PROG_GREP
AC_PROG_SED
AM_PROG_VALAC
dnl LT_INIT([disable-shared])
dnl AC_DISABLE_SHARED
LT_INIT

AX_COMPILER_VENDOR()
AX_COMPILER_VERSION()

dnl ***************************************************************************
dnl C keywords
dnl ***************************************************************************
AC_C_CONST()
AC_C_RESTRICT()
AC_C_VOLATILE()

dnl ***************************************************************************
dnl Makefile.decl variable dependency (here glib_native_win32 is null)
dnl ***************************************************************************
AM_CONDITIONAL(OS_UNIX, [test "$glib_native_win32" != "yes"])

dnl ***************************************************************************
dnl Gobject Introspection
dnl ***************************************************************************
GOBJECT_INTROSPECTION_CHECK([0.6.7])

dnl ***************************************************************************
dnl Compiler flags for C
dnl ***************************************************************************
NC_CFLAGS="$CFLAGS"
ac_[]_AC_LANG_ABBREV[]_werror_flag=yes; dnl Push AC_LANG_WERROR    
AC_ARG_ENABLE([opt-cflags], 
  [AS_HELP_STRING([--enable-opt-cflags  Enable all safe and useful speedups compiler flags (change any -O? to -O3 from CFLAGS)])], 
  [enable_comp_flags="yes"],[])
if test x$enable_comp_flags = "xyes"; then
  NC_CFLAGS=$CFLAGS
  NC_CFLAGS=`echo $NC_CFLAGS | $SED 's/\-O[[0123456789]][[ ]]*//g'`
dnl  NC_CFLAGS=`echo $NC_CFLAGS | $SED 's/\-g[[ ]]*//g'`
  CFLAGS=""
  ac_test_CFLAGS_old=$ac_test_CFLAGS
  ac_test_CFLAGS=""
  AX_CC_MAXOPT()
  CFLAGS="$NC_CFLAGS $CFLAGS"
  CFLAGS=`echo $CFLAGS | $SED 's/\-ffast-math[[ ]]*//g'`
  ac_test_CFLAGS=$ac_test_CFLAGS_old
fi

AC_ARG_ENABLE([warn-cflags], 
  [AS_HELP_STRING([--enable-warn-cflags  Enable all useful warning compiler flags])], 
  [enable_warn_flags="yes"],[])
if test x$enable_warn_flags = "xyes"; then
  NC_CFLAGS=`echo $NC_CFLAGS | $SED 's/\-Wall//g'`
  AX_CHECK_COMPILE_FLAG([-Wall],[NC_CFLAGS="$NC_CFLAGS -Wall"])
dnl  AX_CHECK_COMPILE_FLAG([-Wundef],[NC_CFLAGS="$NC_CFLAGS -Wundef"])
dnl  AX_CHECK_COMPILE_FLAG([-Wshadow],[NC_CFLAGS="$NC_CFLAGS -Wshadow"])
dnl  AX_CHECK_COMPILE_FLAG([-Wpointer-arith],[NC_CFLAGS="$NC_CFLAGS -Wpointer-arith"])
dnl  AX_CHECK_COMPILE_FLAG([-Wcast-align],[NC_CFLAGS="$NC_CFLAGS -Wcast-align"])
dnl  AX_CHECK_COMPILE_FLAG([-Wstrict-prototypes],[NC_CFLAGS="$NC_CFLAGS -Wstrict-prototypes"])
dnl  AX_CHECK_COMPILE_FLAG([-Wwrite-strings],[NC_CFLAGS="$NC_CFLAGS -Wwrite-strings"])
dnl  AX_CHECK_COMPILE_FLAG([-Wunreachable-code],[NC_CFLAGS="$NC_CFLAGS -Wunreachable-code"])
fi

PLC_NOWARN_FLAGS=""
AX_CHECK_COMPILE_FLAG([-Wno-unused-variable], [PLC_NOWARN_FLAGS="$PLC_NOWARN_FLAGS -Wno-unused-variable"])
AX_CHECK_COMPILE_FLAG([-Wno-unused-result],   [PLC_NOWARN_FLAGS="$PLC_NOWARN_FLAGS -Wno-unused-result"])

case $ax_cv_c_compiler_vendor in #(
  gnu) 
    AX_CHECK_COMPILE_FLAG([-Wno-maybe-uninitialized],    [PLC_NOWARN_FLAGS="$PLC_NOWARN_FLAGS -Wno-maybe-uninitialized"])
    AX_CHECK_COMPILE_FLAG([-Wno-unused-but-set-variable],[PLC_NOWARN_FLAGS="$PLC_NOWARN_FLAGS -Wno-unused-but-set-variable"])
  ;;
esac

AC_SUBST(PLC_NOWARN_FLAGS)

AC_ARG_ENABLE([simd-ext-flags], 
  [AS_HELP_STRING([--enable-simd-ext-flags  Enable all SIMD and extension flags])], 
  [enable_simd_ext_flags="yes"],[])
  
if test x$enable_simd_ext_flags = "xyes"; then
  AX_EXT()
  AC_SUBST(SIMD_FLAGS)
  AC_SUBST(CPUEXT_FLAGS)
  CFLAGS="$CFLAGS $SIMD_FLAGS $CPUEXT_FLAGS"
fi

dnl checking if -fvisibility=hidden is available
VISIBILITY_HIDDEN_FLAG=""
AX_CHECK_COMPILE_FLAG([-fvisibility=hidden],[VISIBILITY_HIDDEN_FLAG="-fvisibility=hidden"])
AC_SUBST(VISIBILITY_HIDDEN_FLAG)

dnl checking if -no-undefined is available
NOUNDEFINED_FLAG=""
AX_CHECK_LINK_FLAG([-Wl,-no-undefined],[NOUNDEFINED_FLAG="-no-undefined"])
AC_SUBST(NOUNDEFINED_FLAG)

ac_[]_AC_LANG_ABBREV[]_werror_flag=; dnl Pop AC_LANG_WERROR

dnl ***************************************************************************
dnl Fortran flags
dnl ***************************************************************************

AX_GFORTRAN_VERSION()
AC_FC_FREEFORM()
AC_FC_LINE_LENGTH([unlimited])

AC_F77_LIBRARY_LDFLAGS()
AC_FC_LIBRARY_LDFLAGS()

if test "x$FLIBS" != x; then
  FLIBS_WL="-Wl"
  for fc_opt in $FLIBS; do
    FLIBS_WL="$FLIBS_WL,${fc_opt}"
  done
  FLIBS="$FLIBS_WL"
fi

if test "x$FCLIBS" != x; then
  FCLIBS_WL="-Wl"
  for fc_opt in $FCLIBS; do
    FCLIBS_WL="$FCLIBS_WL,${fc_opt}"
  done
  FCLIBS="$FCLIBS_WL"
  LIBS="$FCLIBS $LIBS"
fi

AC_F77_WRAPPERS()

if test x$GFORTRAN_VERSION != x; then
  GFORTRAN_VERSION_MAJOR=$(echo $GFORTRAN_VERSION | cut -d. -f1)
fi

AC_LANG_PUSH([Fortran])
dnl  AX_APPEND_COMPILE_FLAGS([-finit-local-zero -finit-integer=0 -finit-real=zero -finit-logical=false -fno-automatic], [FCFLAGS])
dnl ***************************************************************************
dnl For intel family compilers we must -diag-error=all to 
dnl make it fails on bad flags (not enough, see below)
dnl ***************************************************************************
  FCFLAGS_CHECK=""
  
  case $FC in #(
    *ifort*) 
      AX_APPEND_COMPILE_FLAGS([-diag-error=all], [FCFLAGS_CHECK])
    ;;
  esac

dnl ***************************************************************************
dnl Ensure any optimization flags we maintain produce no warnings
dnl Required because Intel ifort barfs on bad flags
dnl ***************************************************************************
  ac_[]_AC_LANG_ABBREV[]_werror_flag=yes; dnl Push AC_LANG_WERROR    
  AX_APPEND_COMPILE_FLAGS([-reentrancy -threads -recursive -automatic -fautomatic -frecursive -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-function -Wno-conversion], [FCFLAGS], [$FCFLAGS_CHECK])

  if test $GFORTRAN_VERSION_MAJOR -gt 4; then
    AX_APPEND_COMPILE_FLAGS([-Wno-tabs], [FCFLAGS], [$FCFLAGS_CHECK])
  else
    AX_APPEND_COMPILE_FLAGS([-Wtabs], [FCFLAGS], [$FCFLAGS_CHECK])
  fi

  AC_OPENMP()
  AC_SUBST(OPENMP_FCFLAGS)
  
  ac_[]_AC_LANG_ABBREV[]_werror_flag=; dnl Pop AC_LANG_WERROR

AC_LANG_POP([Fortran])

dnl ***************************************************************************
dnl Check for vapigen
dnl ***************************************************************************
AC_ARG_ENABLE(vapi,
              [AS_HELP_STRING([--enable-vapi],[Generate vala bindings, experimental [default=no]])],enable_vapi=yes,
              enable_vapi=no)

if test "${enable_vapi}" != no; then
  AC_PATH_PROG(VAPIGEN, vapigen, [no])
  AC_SUBST(VAPIGEN)
  VAPIDIR="${datadir}/vala/vapi"
  AC_SUBST(VAPIDIR)
fi

AM_CONDITIONAL(HAVE_VAPIGEN, [test "x$VAPIGEN" != xno -a x$enable_vapi != xno])

dnl ***************************************************************************
dnl check for gtk-doc
dnl ***************************************************************************
m4_ifdef([GTK_DOC_CHECK], [
GTK_DOC_CHECK([1.14],[--flavour no-tmpl])
],[
AM_CONDITIONAL([ENABLE_GTK_DOC], false)
])

AC_ARG_ENABLE(man,
              [AS_HELP_STRING([--enable-man],[regenerate man pages from Docbook [default=no]])],enable_man=yes,
              enable_man=no)

if test "${enable_man}" != no; then

dnl ***************************************************************************
dnl Check for xsltproc
dnl ***************************************************************************
AC_PATH_PROG([XSLTPROC], [xsltproc])
  if test -z "$XSLTPROC"; then
    enable_man=no
  fi
fi

AM_CONDITIONAL(ENABLE_MAN, test x$enable_man != xno)

dnl ***************************************************************************
dnl Set PACKAGE_DATA_DIR in config.h.
dnl ***************************************************************************
if test "x${prefix}" = "xNONE"; then
  PACKAGE_DATA_DIR="${ac_default_prefix}/share/${PACKAGE}-${VERSION}"
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${ac_default_prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
else
  PACKAGE_DATA_DIR="${prefix}/share/${PACKAGE}-${VERSION}"
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
fi
AC_SUBST(PACKAGE_DATA_DIR)

dnl ***************************************************************************
dnl Set PACKAGE_SOURCE_DIR in config.h.
dnl ***************************************************************************
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}",[PACKAGE_SOURCE_DIR])

dnl ***************************************************************************
dnl Dependencies
dnl ***************************************************************************
GLIB_REQUIRED_VERSION=2.38.0
GSL_REQUIRED_VERSION=2.0
FFTW3_REQUIRED_VERSION=3.1.2
FFTW3_ALLOC_VERSION=3.3.2
CFITSIO_REQUIRED_VERSION=3.25

AC_SUBST(GLIB_REQUIRED_VERSION)
AC_SUBST(GSL_REQUIRED_VERSION)

dnl ***************************************************************************
dnl Internatinalization
dnl ***************************************************************************
dnl GETTEXT_PACKAGE=numcosmo
dnl AC_SUBST(GETTEXT_PACKAGE)
dnl AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,"$GETTEXT_PACKAGE", [GETTEXT package name])
dnl AM_GLIB_GNU_GETTEXT

dnl ***************************************************************************
dnl Check for inline function usage
dnl ***************************************************************************

AC_C_INLINE

AC_MSG_CHECKING(--enable-inline)
AC_ARG_ENABLE(inline,
  [AS_HELP_STRING([--enable-inline],[Enable inline functions, disable to debug [[default=yes]]])],
  [enable_inline=$enableval],[enable_inline="yes"])
AC_MSG_RESULT($enable_inline)

if test "x$enable_inline" = "xyes"; then
  if test "$ac_cv_c_inline" != no ; then
    AC_DEFINE([HAVE_INLINE], [1], [use inline functions in GSL])
    AC_SUBST(HAVE_INLINE)
    have_inline_support="#define NUMCOSMO_HAVE_INLINE 1"
    have_inline_macro="#define NCM_INLINE static inline"
    AC_SUBST(have_inline_support)
    AC_SUBST(have_inline_macro)
  fi
else
  AC_DEFINE([G_IMPLEMENT_INLINES], [1], [force glib to disable inline])
  AC_SUBST(G_IMPLEMENT_INLINES)
fi

dnl ***************************************************************************
dnl GCC vector extensions
dnl ***************************************************************************

support_for_gcc_vector_extensions=no

AC_MSG_CHECKING([for GCC vector extensions])
AC_LINK_IFELSE([AC_LANG_SOURCE([[
unsigned int __attribute__ ((vector_size(16))) e, a, b;
int main (void) { e = a - ((b << 27) + (b >> (32 - 27))) + 1; return e[0]; }
]])], support_for_gcc_vector_extensions=yes)

if test x$support_for_gcc_vector_extensions = xyes; then
  AC_DEFINE([HAVE_GCC_VECTOR_EXTENSIONS], [],
            [Whether the compiler supports GCC vector extensions])
fi

AC_MSG_RESULT($support_for_gcc_vector_extensions)

dnl ***************************************************************************
dnl Disable GObject cast checks for speed improvements
dnl ***************************************************************************

AC_MSG_CHECKING(--enable-gobject-cast-check)
AC_ARG_ENABLE(gobject-cast-check,
  [AS_HELP_STRING([--enable-gobject-cast-check],[Enable GObject cast checks, good to debug but slow])],
  [enable_gobject_cast_check=$enableval],[enable_gobject_cast_check="no"])
AC_MSG_RESULT($enable_gobject_cast_check)
if test "x$enable_gobject_cast_check" = "xno"; then
  AC_DEFINE([G_DISABLE_CAST_CHECKS], [1], [Disable GObject cast checks])
  AC_SUBST(G_DISABLE_CAST_CHECKS)    
fi

dnl ***************************************************************************
dnl Enable range checking for debug purposes
dnl ***************************************************************************

AC_MSG_CHECKING(--enable-range-check)
AC_ARG_ENABLE(range-check,
  [AS_HELP_STRING([--enable-range-check],[Enable range checking for debug purposes])],
  [enable_range_check=$enableval],[enable_range_check="no"])
AC_MSG_RESULT($enable_range_check)
if test "x$enable_range_check" = "xno"; then
  AC_DEFINE([GSL_RANGE_CHECK_OFF], [], [optimize gsl access functions])
  AC_SUBST(GSL_RANGE_CHECK_OFF)
fi

dnl ***************************************************************************
dnl Check for multi-thread.
dnl ***************************************************************************

AX_COUNT_CPUS()
AC_MSG_CHECKING([Number of cores])
if test "x$CPU_COUNT" = "x"; then
  AC_MSG_RESULT([cannot detect, using default `4'])
  CPU_COUNT=4
else
  AC_MSG_RESULT([$CPU_COUNT])
fi

AC_DEFINE_UNQUOTED([NCM_NCORES],[$CPU_COUNT], [Number of cores of the current machine])

AC_ARG_WITH([thread-pool-max],
            AS_HELP_STRING([--with-thread-pool-max=N],[Number of threads in the pool]), 
            [thread_pool_max='yes'], [thread_pool_max='no'])

AC_MSG_CHECKING([Cheking for the number of threads in the global thread pool])
if test "x$thread_pool_max" = "xyes"; then
  AC_MSG_RESULT([$with_thread_pool_max])
  AC_DEFINE_UNQUOTED([NCM_THREAD_POOL_MAX],[$with_thread_pool_max], [Number of threads in the pool])
else
  AC_MSG_RESULT([not defined, using number of cores: $CPU_COUNT])
  AC_DEFINE_UNQUOTED([NCM_THREAD_POOL_MAX],[$CPU_COUNT], [Number of threads in the pool])  
fi

AC_OPENMP()
AC_SUBST(OPENMP_CFLAGS)

AM_CONDITIONAL(COMPILE_WITH_MP, [test "x$ac_cv_prog_c_openmp" != "xunsupported"])

if test "x$ac_cv_prog_c_openmp" = "xunsupported"; then
  OPENMP_FCFLAGS=""
fi

dnl ***************************************************************************
dnl Check for backtrace support.
dnl ***************************************************************************
AC_CHECK_HEADERS([execinfo.h],[
		 AC_CHECK_FUNCS([backtrace backtrace_symbols])
],[])

dnl ***************************************************************************
dnl Check for math.h
dnl ***************************************************************************

AC_CHECK_HEADERS([math.h],, AC_MSG_ERROR([math.h header files missing]))
AC_CHECK_LIB([m],[cos])
AC_CHECK_FUNCS([cos sin sincos erf powl exp10 fma finite lgamma_r])
AC_CHECK_DECLS([isfinite, lgamma_r],[],[],[[#include <math.h>]])

dnl ***************************************************************************
dnl Check for dlfcn.h
dnl ***************************************************************************

AC_CHECK_HEADERS([dlfcn.h],, AC_MSG_ERROR([dlfcn.h header files missing]))
AC_CHECK_DECL([RTLD_DEFAULT],
	[has_RTLD_DEFAULT="yes"],
	[has_RTLD_DEFAULT="no"],
	[#include <dlfcn.h>])

dnl ***************************************************************************
dnl Checks for internal libcuba
dnl ***************************************************************************

AC_FUNC_ALLOCA
AC_FUNC_FORK
AC_TYPE_LONG_DOUBLE_WIDER
ac_cv_c_long_double=$ac_cv_type_long_double_wider
if test $ac_cv_c_long_double = yes; then
  AC_DEFINE([HAVE_LONG_DOUBLE],[1],[Define to 1 if the type `long double' works and has more range or precision than `double'.])
fi
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_PID_T

AC_DEFINE([NOUNDERSCORE],[1],[libcuba C interface])

MAXDIM=${MAXDIM:-16}
AC_ARG_WITH(maxdim, 
  [AS_HELP_STRING([--with-libcuba-maxdim=N],
    [the maximum dimension for integration,
     if variable-size array are not supported])],
  [MAXDIM=$withval])

MAXCOMP=${MAXCOMP:-4}
AC_ARG_WITH(maxcomp,
  [AS_HELP_STRING([--with-libcuba-maxcomp=N],
    [the maximum number of components of the integrand,
     if variable-size array are not supported])],
  [MAXCOMP=$withval])

AC_MSG_CHECKING([for variable-size arrays])
AC_COMPILE_IFELSE([AC_LANG_SOURCE(,[[
        void test(int n) {
          char s[n];
        }
  ]])],  
  [AC_MSG_RESULT([yes])],
  [AC_MSG_RESULT([no, using MAXDIM=$MAXDIM and MAXCOMP=$MAXCOMP])
   AC_DEFINE_UNQUOTED([NDIM], [$MAXDIM], [Maximum number of components])
   AC_DEFINE_UNQUOTED([NCOMP], [$MAXCOMP], [Maximum number of dimensions])]
)

dnl ***************************************************************************
dnl Check for glib.
dnl ***************************************************************************

PKG_CHECK_MODULES(GLIB, [glib-2.0 >= $GLIB_REQUIRED_VERSION gthread-2.0 >= $GLIB_REQUIRED_VERSION gobject-2.0 >= $GLIB_REQUIRED_VERSION gio-2.0 >= $GLIB_REQUIRED_VERSION])
AC_SUBST(GLIB_CFLAGS)
AC_SUBST(GLIB_LIBS)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $GLIB_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $GLIB_CFLAGS"

# Setup GLIB_MKENUMS to use glib-mkenums even if GLib is uninstalled.
GLIB_MKENUMS=`$PKG_CONFIG --variable=glib_mkenums glib-2.0`
AC_SUBST(GLIB_MKENUMS)

dnl ***************************************************************************
dnl Check for gtk+-3.0.
dnl ***************************************************************************

PKG_CHECK_MODULES(GTK, [gtk+-3.0 >= 3.6],,[GTK='no'])
AC_SUBST(GTK_CFLAGS)
AC_SUBST(GTK_LIBS)
AM_CONDITIONAL(HAVE_GTK, [test "x$GTK" != xno])
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $GTK_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $GTK_CFLAGS"

dnl ***************************************************************************
dnl Check for gsl
dnl ***************************************************************************

PKG_CHECK_MODULES(GSL, [gsl >= $GSL_REQUIRED_VERSION])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $GSL_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $GSL_CFLAGS"

PKG_CHECK_EXISTS([gsl >= 2.2], [AC_DEFINE([HAVE_GSL_2_2],[1], [gsl 2.2])], [])
PKG_CHECK_EXISTS([gsl >= 2.4], [AC_DEFINE([HAVE_GSL_2_4],[1], [gsl 2.4])], [])

dnl ***************************************************************************
dnl Check for gmp and mpfr
dnl ***************************************************************************

AC_CHECK_FUNC([__gmpz_init],[],[AC_CHECK_LIB([gmp], [__gmpz_init],,[AC_MSG_ERROR([GMP not found. See README.])])])
AC_CHECK_FUNC([mpfr_init],[],[AC_CHECK_LIB([mpfr], [mpfr_init],,[AC_MSG_ERROR([MPFR not found. See README.])])])

dnl ***************************************************************************
dnl Check for arb library
dnl ***************************************************************************

AC_CHECK_HEADERS([arb.h acb.h])
AC_CHECK_LIB([arb],[arb_set])
if test "x$ac_cv_lib_arb_arb_set" = "x"; then
  AC_CHECK_LIB([flint-arb],[arb_set])
  if test "x$ac_cv_lib_flint_arb_arb_set" != "x"; then
    HAVE_ARB=yes
  fi
else
  HAVE_ARB=yes  
fi

if test "x$HAVE_ARB" = "xyes"; then
  AC_CHECK_FUNCS([arb_set acb_set acb_lgamma acb_exp acb_sin])
fi

dnl ***************************************************************************
dnl Check for hdf5 support
dnl ***************************************************************************
AX_LIB_HDF5()
AC_SUBST(HDF5_CPPFLAGS)
AC_SUBST(HDF5_CFLAGS)
AC_SUBST(HDF5_LDFLAGS)
AC_SUBST(HDF5_LIBS)

dnl ***************************************************************************
dnl Check for fftw3
dnl ***************************************************************************

PKG_CHECK_MODULES(FFTW3, [fftw3 >= $FFTW3_REQUIRED_VERSION],[
  AC_DEFINE([HAVE_FFTW3],[1], [Have the fast fourier package])
  have_fftw3_support="#define NUMCOSMO_HAVE_FFTW3 1"],[
  have_fftw3_support=""])
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
AC_SUBST(have_fftw3_support)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $FFTW3_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $FFTW3_CFLAGS"

have_fftw3_alloc="no"
AC_MSG_CHECKING([fftw3 alloc])
PKG_CHECK_EXISTS([fftw3 >= $FFTW3_ALLOC_VERSION], [
  AC_DEFINE([HAVE_FFTW3_ALLOC],[1], [fftw has alloc functions])
  have_fftw3_alloc="yes"], [])
AC_MSG_RESULT($have_fftw3_alloc)

AC_MSG_CHECKING([fftw3 float])
PKG_CHECK_EXISTS([fftw3f >= $FFTW3_REQUIRED_VERSION], [
  AC_MSG_RESULT(yes)
  PKG_CHECK_MODULES(FFTW3F, [fftw3f >= $FFTW3_REQUIRED_VERSION],[
    AC_DEFINE([HAVE_FFTW3F],[1], [Have the fast fourier package (float)])
    have_fftw3f_support="#define NUMCOSMO_HAVE_FFTW3F 1"],[
    have_fftw3f_support=""
    have_fftw3_float="yes"])
  AC_SUBST(FFTW3F_CFLAGS)
  AC_SUBST(FFTW3F_LIBS)
  AC_SUBST(have_fftw3f_support)
  TEST_SHARED_LIBS="$TEST_SHARED_LIBS $FFTW3F_LIBS"
  TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $FFTW3F_CFLAGS"
  ],[
  AC_MSG_RESULT(no)
  ])

dnl ***************************************************************************
dnl Check for BLAS
dnl ***************************************************************************

AX_BLAS([],[have_blas="no"])
AC_MSG_CHECKING([blas support])
if test "x$have_blas" != "xno"; then
  blas_mkl=`echo $BLAS_LIBS | grep -i mkl`
  blas_atlas=`echo $BLAS_LIBS | grep -i atlas`
  blas_openblas=`echo $BLAS_LIBS | grep -i openblas`
  
  if test "x$blas_mkl" != "x"; then
    blas_support="Intel's MKL"
  elif test "x$blas_atlas" != "x"; then
    blas_support="ATLAS"
  elif test "x$blas_openblas" != "x"; then 
    blas_support="OpenBLAS"
  else
    blas_support="Other"
  fi
  AC_MSG_RESULT([$blas_support])

  if test "x$blas_mkl" != "x"; then
    AC_CHECK_HEADERS([mkl_cblas.h],[],[AC_MSG_ERROR([unable to find mkl_cblas.h])])
    AC_CHECK_TYPES([CBLAS_ORDER],[blas_types=""],[blas_types="#define NCM_BLAS_NOT_TYPEDEFED 1"],[[#include <mkl_cblas.h>]])
    numcosmo_blas_include="\
#include <mkl_cblas.h>
#define CBLAS_H
#define __GSL_CBLAS_H__
$blas_types
#include <numcosmo/math/ncm_gsl_blas_types.h>
"
  else
    AC_CHECK_HEADERS([cblas.h],[have_cblas_h="yes"],[have_cblas_h="no"])
    if test "x$have_cblas_h" = "xyes"; then
      AC_CHECK_TYPES([CBLAS_ORDER],[blas_types=""],[blas_types="#define NCM_BLAS_NOT_TYPEDEFED 1"],[[#include <cblas.h>]])
      numcosmo_blas_include="\
#include <cblas.h>
#define __MKL_CBLAS_H__
#define __GSL_CBLAS_H__
$blas_types
#include <numcosmo/math/ncm_gsl_blas_types.h>
"
    else
      AC_DEFINE([HAVE_GSL_CBLAS_H], [1], [use GSL's BLAS header])
      numcosmo_blas_include="\
#include <gsl/gsl_cblas.h>
#define __MKL_CBLAS_H__
#define CBLAS_H
#define NCM_BLAS_NOT_TYPEDEFED 1
#include <numcosmo/math/ncm_gsl_blas_types.h>
"
    fi
  fi
  LIBS="$BLAS_LIBS $LIBS"

  AC_CHECK_FUNC([cblas_dnrm2], [have_cblas="yes"], [])
  AC_MSG_CHECKING([if cblas is already supported])

  if test "x$have_cblas" != "xyes"; then
    AC_MSG_RESULT([no])
    AC_SEARCH_LIBS([cblas_dgemm],[cblas],[have_cblas="yes"],[AC_MSG_ERROR([unable to link to blas])])  
  else
    AC_MSG_RESULT([yes])
  fi
  
  AC_MSG_CHECKING([which cblas to use])
  if test "x$have_cblas" = "xyes"; then
    AC_MSG_RESULT([cblas, removing gslcblas])
    GSL_LIBS=`echo $GSL_LIBS | $SED 's/\-lgslcblas//g'`
  else
    AC_MSG_RESULT([none found, using gslcblas])
fi  
else
  AC_MSG_RESULT([GSL's BLAS])
  AC_DEFINE([HAVE_GSL_CBLAS_H], [1], [use GSL's BLAS header])
  AC_DEFINE([HAVE_BLAS], [1], [use GSL's BLAS])
  numcosmo_blas_include="\
#include <gsl/gsl_cblas.h>
#define __MKL_CBLAS_H__
#define CBLAS_H
#define NCM_BLAS_NOT_TYPEDEFED 1
#include <numcosmo/math/ncm_gsl_blas_types.h>
"    
fi
AC_SUBST(numcosmo_blas_include)

AC_CHECK_FUNCS([openblas_set_num_threads MKL_Set_Num_Threads])

AC_ARG_WITH(lapacke,
        [AS_HELP_STRING([--with-lapacke=<lib>], [use LAPACKe library <lib>])])
case $with_lapacke in
        yes | "") ;;
        no) lapacke_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACKE_LIBS="$with_lapacke" ;;
        *) LAPACKE_LIBS="-l$with_lapacke" ;;
esac

AX_LAPACK([],[have_lapack="no"])
AC_MSG_CHECKING([if lapack is available $LAPACK_LIBS])
if test "x$have_lapack" != "xno"; then
  AC_MSG_RESULT([yes])
  AC_MSG_CHECKING([lapack support])
  lapack_mkl=`echo $BLAS_LIBS $LAPACK_LIBS | grep -i mkl`
  lapack_atlas=`echo $BLAS_LIBS $LAPACK_LIBS | grep -i atlas`
  lapack_openblas=`echo $BLAS_LIBS $LAPACK_LIBS | grep -i openblas`

  if test "x$lapack_mkl" != "x"; then
    lapack_support="Intel's MKL"
  elif test "x$lapack_atlas" != "x"; then
    lapack_support="ATLAS"
  elif test "x$lapack_openblas" != "x"; then
    lapack_support="OpenBLAS"
  else
    lapack_support="Other"
  fi

  LIBS="$LAPACK_LIBS $LIBS"
  AC_MSG_RESULT([$lapack_support])
  
  AC_CHECK_FUNCS([dptsv_ dpotrf_ dpotri_ dpotrs_ dposv_ dsytrf_ dsytrs_ dsytri_ dsysvx_ dsysvxx_ dsyevr_ dsyevd_ dgeev_ dgeevx_ dggglm_ dgeqrf_ dgerqf_ dgeqlf_ dgelqf_ dtrsv_],[],[])

  if test "x$lapack_mkl" != "x"; then
    AC_CHECK_HEADERS([mkl_lapacke.h],
      [AC_CHECK_FUNC([LAPACKE_dptsv],
        [AC_DEFINE([HAVE_LAPACKE],[1], [have lapacke support])
         AC_DEFINE([HAVE_MKL_LAPACKE_H],[1], [have mkl_lapacke.h header])
         have_lapacke="yes"][])])
    AC_CHECK_HEADERS([mkl_lapack.h])
  fi

  if test "x$have_lapacke" != "xyes"; then
    AC_CHECK_HEADERS([lapacke.h],
      [AC_SEARCH_LIBS([LAPACKE_dptsv],[lapacke],  
        [AC_DEFINE([HAVE_LAPACKE],[1], [have lapacke support])
         AC_DEFINE([HAVE_LAPACKE_H],[1], [have lapacke.h header])
         have_lapacke="yes"],[])])
  fi

  if test "x$have_lapacke" != "xyes"; then
    AC_CHECK_HEADERS([clapack.h],
      [AC_SEARCH_LIBS([clapack_dpotri],[clapack],
        [AC_DEFINE([HAVE_CLAPACK],[1], [have clapack support])
         AC_DEFINE([HAVE_CLAPACK_H],[1], [have clapack.h header])
         have_clapack="yes"],[])])
  fi
  
  if test "x$have_lapacke" = "xyes"; then
    lapack_support="$lapack_support+LAPACKe"
  fi

  if test "x$have_clapack" = "xyes"; then
    lapack_support="$lapack_support+cLAPACK"
  fi

else
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([Lapack not found.])
fi

dnl ***************************************************************************
dnl Check for cfitsio
dnl ***************************************************************************

PKG_CHECK_MODULES(CFITSIO, [cfitsio >= $CFITSIO_REQUIRED_VERSION],[
  AC_DEFINE([HAVE_CFITSIO],[1], [Have the cfitsio package])
  have_cfitsio_support="#define NUMCOSMO_HAVE_CFITSIO 1"])
AC_SUBST(CFITSIO_CFLAGS)
AC_SUBST(CFITSIO_LIBS)
AC_SUBST(have_cfitsio_support)

SAVED_LIBS="$LIBS"
LIBS="$CFITSIO_LIBS $LIBS"
SAVED_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $CFITSIO_CFLAGS -Werror"

AC_MSG_CHECKING([fits_free_memory presence])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <fitsio.h>
]],[[
fits_free_memory ((void *)0, (void *)0);
]])],[AC_MSG_RESULT([yes])],[
AC_DEFINE([fits_free_memory(p,s)],[fffree(p,s)], [fits_free_memory])
AC_MSG_RESULT([no])
])

AC_MSG_CHECKING([fffree presence])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <fitsio.h>
]],[[
fffree ((void *)0, (void *)0);
]])],[AC_MSG_RESULT([yes])],[
AC_DEFINE([fffree(p,s)],[g_free(p)], [fffree])
AC_MSG_RESULT([no])
])

LIBS="$SAVED_LIBS"
CFLAGS="$SAVED_CFLAGS"

dnl ***************************************************************************
dnl Check for chealpix
dnl ***************************************************************************
dnl 
dnl AC_CHECK_LIB([chealpix], [read_healpix_map],,[chealpix_missing="yes"],[])
dnl if test "x$chealpix_missing" = "xyes"; then
dnl   have_chealpix_support=""
dnl else
dnl   have_chealpix_support="#define NUMCOSMO_HAVE_CHEALPIX 1"
dnl fi
dnl AC_SUBST(have_chealpix_support)
dnl 
dnl ***************************************************************************
dnl Check for libcuba
dnl ***************************************************************************

AC_CHECK_LIB([cuba], [Cuhre],,[cuba_missing="yes"],[])
AM_CONDITIONAL([USE_INTERNAL_LIBCUBA], [test "x$cuba_missing" = "xyes"])

AC_MSG_CHECKING([libcuba])
if test "x$cuba_missing" = "xyes"; then
  have_cuba_support=""
  AC_MSG_RESULT([internal v4.2])
  AC_DEFINE([HAVE_LIBCUBA_4_0],[1], [Have the libcuba >= 4.0])
  AC_DEFINE([HAVE_LIBCUBA],[1], [Have any version of libcuba])
else
  AC_MSG_RESULT([external])
  AC_MSG_CHECKING([libcuba version])
  AC_COMPILE_IFELSE(
              [AC_LANG_PROGRAM([[#include <cuba.h>]],
                 [[Cuhre (0, 0, (void *)0, (void *)0, 0.0, 0.0, 0, 0, 0, 0, (void *)0 /*AQUI*/,  (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0);]])],
              [AC_DEFINE([HAVE_LIBCUBA_3_1],[1], [Have the libcuba >= 3.1])
               AC_DEFINE([HAVE_LIBCUBA],[1], [Have any version of libcuba])
               cuba_3_1="<= 3.3"],
              [cuba_3_1="error"
  AC_COMPILE_IFELSE(
              [AC_LANG_PROGRAM([[#include <cuba.h>]],
                 [[Cuhre (0, 0, (void *)0, (void *)0, 0, 0.0, 0.0, 0, 0, 0, 0, (void *)0 /*AQUI*/,  (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0);]])],
              [AC_DEFINE([HAVE_LIBCUBA_3_3],[1], [Have the libcuba >= 3.3])
               AC_DEFINE([HAVE_LIBCUBA],[1], [Have any version of libcuba])
               cuba_3_1="<= 3.3"],
              [cuba_3_1="error"
  AC_COMPILE_IFELSE(
              [AC_LANG_PROGRAM([[#include <cuba.h>]],
                 [[Cuhre (0, 0, (void *)0, (void *)0, 0, 0.0, 0.0, 0, 0, 0, 0, (void *)0 /*AQUI*/,  /*AQUI*/ (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0);]])],
              [AC_DEFINE([HAVE_LIBCUBA_4_0],[1], [Have the libcuba >= 4.0])
               AC_DEFINE([HAVE_LIBCUBA],[1], [Have any version of libcuba])
               cuba_3_1=">= 4.0"],
              [cuba_3_1="error"
              ])])])
  AC_MSG_RESULT([$cuba_3_1])
fi

dnl ***************************************************************************
dnl Check for nlopt
dnl ***************************************************************************

PKG_CHECK_MODULES(NLOPT, [nlopt >= 1.2],[
  AC_DEFINE([HAVE_NLOPT],[1], [Have the NLopt package])
  have_nlopt_support="#define NUMCOSMO_HAVE_NLOPT 1"],[
  have_nlopt_support=""])
AC_SUBST(NLOPT_CFLAGS)
AC_SUBST(NLOPT_LIBS)
AC_SUBST(have_nlopt_support)

AM_CONDITIONAL([HAVE_LIBNLOPT], [test "x$have_nlopt_support" != x])

TEST_SHARED_LIBS="$TEST_SHARED_LIBS $NLOPT_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $NLOPT_CFLAGS"

if test "x$have_nlopt_support" != "x"; then

PKG_CHECK_EXISTS([nlopt >= 2.2], [
AC_DEFINE([HAVE_NLOPT_2_2],[1], [If nlopt version is bigger than 2.2])
AC_CHECK_LIB([m],[nlopt_create],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with c++ support, recompile the C version.])],[$NLOPT_LIBS])
],[
AC_CHECK_LIB([m],[nlopt_minimize_constrained],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with C++ support, recompile the C version.])],[$NLOPT_LIBS])
])

fi

dnl ***************************************************************************
dnl Check for sundials
dnl ***************************************************************************

dnl AC_ARG_VAR([SUNDIALS_LIBS], [C compiler flags for SUNDIALS])
dnl AC_ARG_VAR([SUNDIALS_CFLAGS], [linker flags for SUNDIALS])
dnl AC_SUBST(SUNDIALS_LIBS)
dnl AC_SUBST(SUNDIALS_CFLAGS)

dnl AC_CHECK_LIB([sundials_nvecserial],[N_VNew_Serial],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_nvecserial"],[AC_MSG_ERROR([Cannot find library sundials_nvecserial, use SUNDIALS_LIBS="-L<path_to_SUNDIALS_libs>".])],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_nvecpthreads],[N_VNew_Pthreads],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_nvecpthreads"],[],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_nvecopenmp],[N_VNew_OpenMP],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_nvecopenmp"],[],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_sunmatrixdense],[SUNDenseMatrix],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_sunmatrixdense"],[],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_sunlinsoldense],[SUNDenseLinearSolver],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_sunlinsoldense"],[],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_cvode],[CVodeCreate],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_cvode"],[AC_MSG_ERROR([Cannot find library sundials_cvode, use SUNDIALS_LIBS="-L<path_to_SUNDIALS_libs>".])],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_cvodes],[CVodeCreate],[SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_cvodes"],[AC_MSG_ERROR([Cannot find library sundials_cvodes, use SUNDIALS_LIBS="-L<path_to_SUNDIALS_libs>".])],[$SUNDIALS_LIBS])
dnl AC_CHECK_LIB([sundials_arkode],[ARKodeCreate],[
dnl 	SUNDIALS_LIBS="$SUNDIALS_LIBS -lsundials_arkode"
dnl 	AC_DEFINE([HAVE_SUNDIALS_ARKODE],[1],[SUNDIALS has ARKode])
dnl 	],[],[$SUNDIALS_LIBS])

dnl SAVED_CFLAGS="$CFLAGS"
dnl CFLAGS="$CFLAGS $SUNDIALS_CFLAGS"

dnl AC_CHECK_HEADER([sundials/sundials_config.h],[],[
dnl AC_MSG_ERROR([Cannot find sundials/sundials_config.h, use SUNDIALS_CFLAGS="-I<path_to_SUNDIALS_headers>".])
dnl ],[])

dnl AC_MSG_CHECKING([SUNDIALS version])

dnl SUNDIALS_VERSION_FOUND=no
dnl m4_foreach([SUNDIALS_VERSION_V], [[2.5.0, 2, 5, 0], [2.6.0, 2, 6, 0], [2.6.1, 2, 6 ,1], [2.6.2, 2, 6 ,2], [2.7.0, 2, 7, 0]],
dnl   [
dnl     SUNDIALS_VERSION=m4_argn([1],SUNDIALS_VERSION_V)
dnl     SUNDIALS_VERSION_REGEX=$(echo $SUNDIALS_VERSION | $SED 's/\./\\./g')
dnl     AC_EGREP_CPP([[\"]?$SUNDIALS_VERSION_REGEX[\"]?],
dnl       [
dnl         #include <sundials/sundials_config.h>
dnl         SUNDIALS_PACKAGE_VERSION;
dnl         SUNDIALS_VERSION;
dnl       ],
dnl       [
dnl         if test "x$SUNDIALS_VERSION_FOUND" != "xno"; then
dnl           AC_MSG_ERROR([Cannot identify SUNDIALS version (2.5.x|2.6.x|2.7.x|3.0.x|3.1.x) double match])
dnl         fi
dnl         HAVE_SUNDIALS_MAJOR=m4_argn([2],SUNDIALS_VERSION_V)
dnl         HAVE_SUNDIALS_MINOR=m4_argn([3],SUNDIALS_VERSION_V)
dnl         HAVE_SUNDIALS_PATCH=m4_argn([4],SUNDIALS_VERSION_V)
dnl         SUNDIALS_VERSION_FOUND=yes
dnl       ])
dnl   ])dnl

dnl if test "x$SUNDIALS_VERSION_FOUND" = "xno"; then
dnl   AC_COMPUTE_INT([HAVE_SUNDIALS_MAJOR],[SUNDIALS_VERSION_MAJOR],[[#include <sundials/sundials_config.h>]],
dnl     AC_MSG_FAILURE([sundials/sundials_config.h: could not determine 'SUNDIALS_VERSION_MAJOR' value]))
dnl   AC_COMPUTE_INT([HAVE_SUNDIALS_MINOR],[SUNDIALS_VERSION_MINOR],[[#include <sundials/sundials_config.h>]],
dnl     AC_MSG_FAILURE([sundials/sundials_config.h: could not determine 'SUNDIALS_VERSION_MINOR' value]))
dnl   AC_COMPUTE_INT([HAVE_SUNDIALS_PATCH],[SUNDIALS_VERSION_PATCH],[[#include <sundials/sundials_config.h>]],
dnl     AC_MSG_FAILURE([sundials/sundials_config.h: could not determine 'SUNDIALS_VERSION_PATCH' value]))
dnl fi

dnl SUNDIALS_VERSION="$HAVE_SUNDIALS_MAJOR.$HAVE_SUNDIALS_MINOR.$HAVE_SUNDIALS_PATCH"
dnl AC_MSG_RESULT([$SUNDIALS_VERSION])

dnl AC_DEFINE_UNQUOTED([HAVE_SUNDIALS_MAJOR],[$HAVE_SUNDIALS_MAJOR],[SUNDIALS version - major])
dnl AC_DEFINE_UNQUOTED([HAVE_SUNDIALS_MINOR],[$HAVE_SUNDIALS_MINOR],[SUNDIALS version - minor])
dnl AC_DEFINE_UNQUOTED([HAVE_SUNDIALS_PATCH],[$HAVE_SUNDIALS_PATCH],[SUNDIALS version - patch])

AC_MSG_CHECKING([SUNDIALS version])
AC_MSG_RESULT([using internal version])

dnl ***************************************************************************
dnl PLC
dnl ***************************************************************************

dnl ***************************************************************************
dnl Using -DLAPACK_CLIK as default (avoid problems with some broken clapack.h)
dnl ***************************************************************************

dnl if test "x$lapack_mkl" = "xyes"; then
dnl   PLC_DEFINES="$PLC_DEFINES -DHAS_MKL"
dnl elif test "x$have_clapack" != "xyes"; then
  PLC_DEFINES="$PLC_DEFINES -DHAS_LAPACK -DLAPACK_CLIK"
dnl fi

PLC_DEFINES="$PLC_DEFINES -DNOHEALPIX -DCLIK_LENSING -DCAMSPEC_V1"
PLC_DEFINES="$PLC_DEFINES -D 'CLIKSVNVERSION=\"0031e3604f0c\"'"

if test "x$has_RTLD_DEFAULT" = "xyes"; then
  PLC_DEFINES="$PLC_DEFINES -DHAS_RTLD_DEFAULT"
fi

AC_SUBST(PLC_DEFINES)

dnl ***************************************************************************
dnl CCL
dnl ***************************************************************************

AC_ARG_VAR([CCL_LDFLAGS], [C linker flags for CCL])
AC_ARG_VAR([CCL_CFLAGS], [C compiler flags for CCL])
AC_SUBST(CCL_LDFLAGS)
AC_SUBST(CCL_CFLAGS)

AC_CHECK_LIB([ccl],[ccl_cosmology_create],[ccl_found=yes],[ccl_found=no],[$CCL_LDFLAGS])

SAVED_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $CCL_CFLAGS"
AC_CHECK_HEADERS([ccl.h],[ccl_header_found=yes],[ccl_header_found=no])
CPPFLAGS="$SAVED_CPPFLAGS"

numcosmo_check_ccl=no
AC_MSG_CHECKING([CCL usability])
if test "x$ccl_found" = "xyes"; then
  CCL_LDFLAGS="$CCL_LDFLAGS -lccl"
  if test "x$ccl_header_found" = "xyes"; then
    CCL_LDFLAGS="$CCL_LDFLAGS -lccl"
    numcosmo_check_ccl=yes
    AC_DEFINE(HAVE_CCL,1,[Define if you have the CCL library.])
    AC_MSG_RESULT([yes])    
  else
    AC_MSG_RESULT([no])
  fi
else
  AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([NUMCOSMO_CHECK_CCL], [test "x$numcosmo_check_ccl" != "xno"])

dnl ***************************************************************************
dnl Appending libs and setting up TEST_SHARED
dnl ***************************************************************************

TEST_SHARED_LIBS="$TEST_SHARED_LIBS $LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS"

dnl ***************************************************************************
dnl Check if we can build a shared library
dnl ***************************************************************************
dnl SAVE_LIBS="$LIBS"
dnl SAVE_CFLAGS="$CFLAGS"
dnl SAVE_LDFLAGS="$LDFLAGS"
dnl LIBS="$TEST_SHARED_LIBS"
dnl CFLAGS="$TEST_SHARED_CFLAGS"
dnl LDFLAGS="$LDFLAGS -shared"
dnl AC_LINK_IFELSE([int main(){}],[],[])
dnl LIBS="$SAVE_LIBS"
dnl CFLAGS="$SAVE_CFLAGS"
dnl LDFLAGS="SAVE_LDFLAGS"
dnl CFLAGS="$CFLAGS $NC_CFLAGS"

dnl ***************************************************************************
dnl Setting up
dnl ***************************************************************************

AC_CONFIG_FILES([
numcosmo.pc
Makefile
numcosmo/Makefile
numcosmo/levmar/Makefile
numcosmo/libcuba/Makefile
numcosmo/libcuba/cuhre/Makefile
numcosmo/libcuba/divonne/Makefile
numcosmo/libcuba/suave/Makefile
numcosmo/libcuba/vegas/Makefile
numcosmo/plc/Makefile
numcosmo/class/Makefile
numcosmo/sundials/Makefile
numcosmo/toeplitz/Makefile
tools/Makefile
tests/Makefile
docs/Makefile
docs/manual_head.xml
examples/Makefile
examples/model_mcmc/Makefile
examples/pydata_simple/Makefile
notebooks/Makefile
data/Makefile
numcosmo/build_cfg.h
])

AC_OUTPUT

echo " "
echo "--------------------------------------------------------------------"
echo " Optional packages:"
echo "--------------------------------------------------------------------"
if [ test "x$have_fftw3_support" != x ]; then
       echo "Building with FFTW3 support: ....................................YES"
else
       echo "Building with FFTW3 support: .....................................NO"
       echo "        Requires FFTW3 (>= $FFTW3_REQUIRED_VERSION)"
fi
       echo "Building with BLAS support: .....................................$blas_support"
if [ test "x$have_lapack" != "xno" ]; then
       echo "Building with Lapack support: ...................................$lapack_support"
else
       echo "Building with Lapack support: ....................................NO"
       echo "        Requires Lapack Linear Algebra PACKage"
fi
if [ test "x$have_cfitsio_support" != x ]; then
       echo "Building with cfitsio support: ..................................YES"
else
       echo "Building with cfitsio support: ...................................NO"
       echo "        Requires cfitsio - A FITS File Subroutine Library"
fi
if [ test "x$have_nlopt_support" != x ]; then
       echo "Building with nlopt support: ....................................YES"
else
       echo "Building with nlopt support: .....................................NO"
       echo "        Requires nlopt - NonLinear OPTimization library"
fi
if [ test "x$with_hdf5" = xyes ]; then
       echo "Building with HDF5 support: .....................................YES"
else
       echo "Building with HDF5 support: ......................................NO"
fi
if [ test "x$HAVE_ARB" = xyes ]; then
       echo "Building with ARB support: ......................................YES"
else
       echo "Building with ARB support: .......................................NO"
fi
echo "--------------------------------------------------------------------"
echo "AM_CFLAGS = $AM_CFLAGS"
echo "CFLAGS = $CFLAGS"
echo "FCFLAGS = $FCFLAGS"
echo "OPENMP_CFLAGS = $OPENMP_CFLAGS"
echo "OPENMP_FCFLAGS = $OPENMP_FCFLAGS"
echo "LIBS = $LIBS"
echo "PLC_DEFINES = $PLC_DEFINES"
echo "PLC_NOWARN_FLAGS = $PLC_NOWARN_FLAGS"
echo "PLC_FC_NOWARN_FLAGS = $PLC_FC_NOWARN_FLAGS"

if [ test "x$HAVE_ARB" = xyes ]; then
  echo "HDF5_CPPFLAGS = $HDF5_CPPFLAGS"
  echo "HDF5_CFLAGS = $HDF5_CFLAGS"
  echo "HDF5_LDFLAGS = $HDF5_LDFLAGS"
  echo "HDF5_LIBS = $HDF5_LIBS"
fi
