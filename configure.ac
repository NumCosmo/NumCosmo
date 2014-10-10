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
m4_define([numcosmo_minor_version], [12])
m4_define([numcosmo_micro_version], [0])
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
          
AC_INIT([numcosmo],[numcosmo_version], https://savannah.nongnu.org/bugs/?func=additem&group=numcosmo, numcosmo, http://www.nongnu.org/numcosmo/)
AC_CONFIG_SRCDIR([numcosmo/numcosmo.h])
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

AM_MAINTAINER_MODE

dnl ***************************************************************************
dnl add the correct m4 directory location.
dnl ***************************************************************************
AC_CONFIG_MACRO_DIR([m4])
  
dnl ***************************************************************************
dnl Gobject Introspection
dnl ***************************************************************************
GOBJECT_INTROSPECTION_CHECK([0.6.7])

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
dnl Checks for programs.
dnl ***************************************************************************
AC_LANG(C)
AC_PROG_CC
AC_PROG_F77
AC_PROG_INSTALL   
AC_PROG_LN_S
AC_PROG_CC_STDC
AC_PROG_GREP
AC_PROG_SED
AC_PROG_FC
AM_PROG_VALAC
dnl LT_INIT([disable-shared])
dnl AC_DISABLE_SHARED
LT_INIT

dnl ***************************************************************************
dnl Makefile.decl variable dependency (here glib_native_win32 is null)
dnl ***************************************************************************
AM_CONDITIONAL(OS_UNIX, [test "$glib_native_win32" != "yes"])

dnl ***************************************************************************
dnl Check for vapigen
dnl ***************************************************************************
AC_ARG_ENABLE(vapi,
              [AC_HELP_STRING([--enable-vapi],
                              [Generate vala bindings, experimental [default=no]])],enable_vapi=yes,
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
GTK_DOC_CHECK([1.14],[--flavour no-tmpl])

AC_ARG_ENABLE(man,
              [AC_HELP_STRING([--enable-man],
                              [regenerate man pages from Docbook [default=no]])],enable_man=yes,
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
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${ac_default_prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
else
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
fi

dnl ***************************************************************************
dnl Set PACKAGE_SOURCE_DIR in config.h.
dnl ***************************************************************************
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}",[PACKAGE_SOURCE_DIR])

dnl ***************************************************************************
dnl Dependencies
dnl ***************************************************************************
GLIB_REQUIRED_VERSION=2.24.0
GSL_REQUIRED_VERSION=1.13
GSL_GLF_REQUIRED_VERSION=1.15
GSL_ODEIV2_REQUIRED_VERSION=1.15
SQLITE3_REQUIRED_VERSION=3.6.10
FFTW3_REQUIRED_VERSION=3.1.2
FFTW3_ALLOC_VERSION=3.3.2

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
    AC_SUBST(have_inline_support)
  fi
else
  AC_DEFINE([G_IMPLEMENT_INLINES], [1], [force glib to disable inline])
  AC_SUBST(G_IMPLEMENT_INLINES)
fi

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
AC_DEFINE_UNQUOTED([NCM_NCORES],[$CPU_COUNT], [Number of cores of the current machine])

AC_ARG_WITH([thread-pool-max],
            AS_HELP_STRING([--with-thread-pool-max=N],[Number of threads in the pool]), 
            [thread_pool_max='yes'], [thread_pool_max='no'])

if test "x$thread_pool_max" = "xyes"; then
  AC_DEFINE_UNQUOTED([NCM_THREAD_POOL_MAX],[$with_thread_pool_max], [Number of threads in the pool])
else
  AC_DEFINE_UNQUOTED([NCM_THREAD_POOL_MAX],[$CPU_COUNT], [Number of threads in the pool])  
fi

dnl ***************************************************************************
dnl Check for fortran lib.
dnl ***************************************************************************
AC_F77_LIBRARY_LDFLAGS
AC_FC_LIBRARY_LDFLAGS
if test "x$FCLIBS" != x; then
  FCLIBS_WL="-Wl"
  for fc_opt in $FCLIBS; do
    FCLIBS_WL="$FCLIBS_WL,${fc_opt}"
  done
  FLIBS_WL="-Wl"
  for fc_opt in $FLIBS; do
    FLIBS_WL="$FLIBS_WL,${fc_opt}"
  done

  FLIBS="$FLIBS_WL"
  FCLIBS="$FCLIBS_WL"
  LIBS="$FCLIBS $LIBS"
fi

dnl ***************************************************************************
dnl Check for math.h
dnl ***************************************************************************

AC_CHECK_HEADERS(math.h,, AC_MSG_ERROR([math.h header files missing]))
AC_CHECK_FUNC([cos],[],[AC_CHECK_LIB([m],[cos])])

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

PKG_CHECK_EXISTS([gsl >= $GSL_GLF_REQUIRED_VERSION], [AC_DEFINE([HAVE_GSL_GLF],[1], [gsl support fixed gauss legendre rules])], [])
PKG_CHECK_EXISTS([gsl >= $GSL_ODEIV2_REQUIRED_VERSION], [AC_DEFINE([HAVE_GSL_ODEIV2],[1], [gsl support odeiv2 ode algorightms])], [])

dnl ***************************************************************************
dnl Check for gmp and mpfr
dnl ***************************************************************************

AC_CHECK_FUNC([__gmpz_init],[],[AC_CHECK_LIB([gmp], [__gmpz_init],,[AC_MSG_ERROR([GMP not found. See README.])])])
AC_CHECK_FUNC([mpfr_init],[],[AC_CHECK_LIB([mpfr], [mpfr_init],,[AC_MSG_ERROR([MPFR not found. See README.])])])

dnl ***************************************************************************
dnl Check for sqlite3
dnl ***************************************************************************

PKG_CHECK_MODULES(SQLITE3, [sqlite3 >= $SQLITE3_REQUIRED_VERSION],[
  AC_DEFINE([HAVE_SQLITE3],[1], [Have the sqlite database])
  have_sqlite3_support="#define NUMCOSMO_HAVE_SQLITE3 1"],[
  have_sqlite3_support=""])
AC_SUBST(SQLITE3_CFLAGS)
AC_SUBST(SQLITE3_LIBS)
AC_SUBST(have_sqlite3_support)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $SQLITE3_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $SQLITE3_CFLAGS"

dnl ***************************************************************************
dnl Check for fftw3
dnl ***************************************************************************

PKG_CHECK_MODULES(FFTW3, [fftw3 >= $FFTW3_REQUIRED_VERSION fftw3l >= $FFTW3_REQUIRED_VERSION],[
  AC_DEFINE([HAVE_FFTW3],[1], [Have the fast fourier package])
  have_fftw3_support="#define NUMCOSMO_HAVE_FFTW3 1"],[
  have_fftw3_support=""])
AC_SUBST(FFTW3_CFLAGS)
AC_SUBST(FFTW3_LIBS)
AC_SUBST(have_fftw3_support)
TEST_SHARED_LIBS="$TEST_SHARED_LIBS $FFTW3_LIBS"
TEST_SHARED_CFLAGS="$TEST_SHARED_CFLAGS $FFTW3_CFLAGS"

PKG_CHECK_EXISTS([fftw3 >= $FFTW3_ALLOC_VERSION], [AC_DEFINE([HAVE_FFTW3_ALLOC],[1], [fftw has alloc functions])], [])

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
    AC_CHECK_HEADERS([mkl_blas.h])
    AC_CHECK_HEADERS([mkl_cblas.h])
  else
    AC_CHECK_HEADERS([cblas.h])
  fi
  LIBS="$BLAS_LIBS $LIBS"

  AC_CHECK_FUNC([cblas_dnrm2], [have_cblas="yes"], [])
  AC_MSG_CHECKING([if cblas is already supported])

  if test "x$have_cblas" != "xyes"; then
    AC_MSG_RESULT([no])
    AC_SEARCH_LIBS([cblas_dgemm],[cblas],[have_cblas="yes"],[])  
  else
    AC_MSG_RESULT([yes])
  fi

  AC_MSG_CHECKING([which cblas to use])
  if test "x$have_cblas" = "xyes"; then
    AC_MSG_RESULT([cblas, removing gslcblas])
    GSL_LIBS=`echo $GSL_LIBS | sed 's/\-lgslcblas//g'`
  else
    AC_MSG_RESULT([none found, using gslcblas])
  fi
else
  AC_MSG_RESULT([GSL's BLAS])
  AC_DEFINE([HAVE_GSL_CBLAS_H], [1], [use lapack functions])
fi

dnl ***************************************************************************
dnl Check for LAPACK
dnl ***************************************************************************

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
  AC_MSG_RESULT([$lapack_support])

  if test "x$lapack_mkl" != "x"; then
    AC_CHECK_HEADERS([mkl_lapacke.h],
      [AC_DEFINE([HAVE_LAPACKE],[1], [have lapacke support])
       AC_DEFINE([HAVE_MKL_LAPACKE_H],[1], [have mkl_lapacke.h header])])
    AC_CHECK_HEADERS([mkl_lapack.h])
  else
    AC_CHECK_HEADERS([lapacke.h],
      [AC_DEFINE([HAVE_LAPACKE],[1], [have lapacke support])
       AC_DEFINE([HAVE_LAPACKE_H],[1], [have lapacke.h header])])
    AC_CHECK_HEADERS([clapack.h],
      [AC_DEFINE([HAVE_CLAPACK],[1], [have clapack support])
       AC_DEFINE([HAVE_CLAPACK_H],[1], [have clapack.h header])])
  fi
  LIBS="$LAPACK_LIBS $LIBS"
else
  AC_MSG_RESULT([no])
fi

dnl ***************************************************************************
dnl Check for levmar
dnl ***************************************************************************

AC_MSG_CHECKING([--with-levmar-path])
AC_ARG_WITH(levmar-path,
  [AS_HELP_STRING([--with-levmar-path=PATH],[Levmar library location])],
  [LEVMAR_PATH="-L$withval"], [])
AC_MSG_RESULT($LEVMAR_PATH)

if test "x$LEVMAR_PATH" != x; then
  LIBS="$LEVMAR_PATH $LIBS"
fi

AC_MSG_CHECKING([--with-levmar-extra-libs])
AC_ARG_WITH(levmar-extra-libs,
  [AS_HELP_STRING([--with-levmar-extra-libs=LIBS],[Levmar extra libraries, default is none])],
  [LEVMAR_EXTRA_LIBS="$withval"], [])
AC_MSG_RESULT($LEVMAR_EXTRA_LIBS)

AC_CHECK_LIB([levmar], [dlevmar_der],[LIBS="-llevmar $LEVMAR_EXTRA_LIBS $LIBS"],[levmar_missing="yes"],[$LEVMAR_EXTRA_LIBS])

AC_CHECK_HEADER(levmar.h, , [AC_CHECK_HEADER(levmar/levmar.h, [levmar_need_prefix="yes"], [levmar_missing="yes"])])
if test "x$levmar_missing" = "xyes"; then
  have_levmar_support=""
else
  have_levmar_support="#define NUMCOSMO_HAVE_LEVMAR 1"
fi
AC_SUBST(have_levmar_support)

AM_CONDITIONAL([HAVE_LIBLEVMAR], [test "x$have_levmar_support" != x])

if test "x$levmar_need_prefix" = "xyes"; then
  AC_DEFINE([NC_LEVMAR_NEED_PREFIX],[1], [Duh!])
fi

dnl ***************************************************************************
dnl Check for cfitsio
dnl ***************************************************************************

PKG_CHECK_MODULES(CFITSIO, [cfitsio >= 3.25],[
  AC_DEFINE([HAVE_CFITSIO],[1], [Have the cfitsio package])
  have_cfitsio_support="#define NUMCOSMO_HAVE_CFITSIO 1"],[
  have_cfitsio_support=""])
AC_SUBST(CFITSIO_CFLAGS)
AC_SUBST(CFITSIO_LIBS)
AC_SUBST(have_cfitsio_support)

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
if test "x$cuba_missing" = "xyes"; then
  have_cuba_support=""
else
  have_cuba_support="#define NUMCOSMO_HAVE_CUBA 1"
fi
AC_SUBST(have_cuba_support)

AC_MSG_CHECKING([if libcuba version is >= 3.1])
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([[#include <cuba.h>]],
               [[Cuhre (0, 0, (void *)0, (void *)0, 0.0, 0.0, 0, 0, 0, 0, (void *)0 /*AQUI*/,  (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0);]])],
            [AC_DEFINE([HAVE_LIBCUBA_3_1],[1], [Have the libcuba >= 3.1])
             cuba_3_1="yes"],
            [cuba_3_1="no"
AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([[#include <cuba.h>]],
               [[Cuhre (0, 0, (void *)0, (void *)0, 0, 0.0, 0.0, 0, 0, 0, 0, (void *)0 /*AQUI*/,  (void *)0, (void *)0, (void *)0, (void *)0, (void *)0, (void *)0);]])],
            [AC_DEFINE([HAVE_LIBCUBA_3_3],[1], [Have the libcuba >= 3.3])
             cuba_3_1="yes"],
            [cuba_3_1="no"
            ])
            ])
AC_MSG_RESULT([$cuba_3_1])

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
AC_CHECK_LIB([nlopt],[nlopt_create],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with c++ support, recompile the C version.])],[])
],[
AC_CHECK_LIB([nlopt],[nlopt_minimize_constrained],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with C++ support, recompile the C version.])],[])
])

fi

dnl ***************************************************************************
dnl Check for sundials
dnl ***************************************************************************

AC_CHECK_PROG(has_sundials, sundials-config, yes, no)
if test "x$has_sundials" = "xyes"; then
  SUNDIALS_LIBS=`sundials-config -s libs -m cvodes -t s -l c`
  SUNDIALS_CFLAGS=`sundials-config -s cppflags -m cvodes -t s -l c`
  AC_SUBST(SUNDIALS_LIBS)
  AC_SUBST(SUNDIALS_CFLAGS)
else
  AC_MSG_ERROR(Cannot find sundials-config.)
fi

SAVED_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS $SUNDIALS_CFLAGS"
AC_CHECK_TYPE([CVDlsDenseJacFn],[],[AC_MSG_ERROR([sundials too old, please install version >= 2.4.0])], [#include <cvodes/cvodes_direct.h>])
CFLAGS="$SAVED_CFLAGS"

AC_MSG_CHECKING([whenever sundials uses long int])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <cvodes/cvodes_direct.h>
typedef int (*CVDlsDenseJacFn)(long int N, realtype t,
                               N_Vector y, N_Vector fy,
                               DlsMat Jac, void *user_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
]],[[1;]])],[AC_DEFINE([SUNDIALS_USES_LONG_INT],[1], [If sundials version >= 2.5.0]) AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no])])

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

dnl ***************************************************************************
dnl Setting up
dnl ***************************************************************************

AC_CONFIG_FILES([
numcosmo.pc
Makefile
numcosmo/Makefile
tools/Makefile
tests/Makefile
docs/Makefile
examples/Makefile
data/Makefile
darkenergy/Makefile
numcosmo/build_cfg.h
])

AC_OUTPUT

echo " "
echo "--------------------------------------------------------------------"
echo " Optional packages:"
echo "--------------------------------------------------------------------"
if [ test "x$have_sqlite3_support" != x ]; then
       echo "Building with SQLite3 support: ..................................YES"
else
       echo "Building  with SQLite3 support: ..................................NO"
       echo "        Requires SQLite3 (>= $SQLITE3_REQUIRED_VERSION)"
fi
if [ test "x$have_fftw3_support" != x ]; then
       echo "Building with FFTW3 support: ....................................YES"
else
       echo "Building  with FFTW3 support: ....................................NO"
       echo "        Requires FFTW3 (>= $FFTW3_REQUIRED_VERSION)"
fi
       echo "Building with BLAS support: .....................................$blas_support"
if [ test "x$have_lapack" != "xno" ]; then
       echo "Building with Lapack support: ...................................$lapack_support"
else
       echo "Building  with Lapack support: ...................................NO"
       echo "        Requires Lapack Linear Algebra PACKage"
fi
if [ test "x$have_levmar_support" != x ]; then
       echo "Building with levmar support: ...................................YES"
else
       echo "Building  with levmar support: ...................................NO"
       echo "        Requires levmar Levenberg-Marquardt least squares library"
fi
if [ test "x$have_cfitsio_support" != x ]; then
       echo "Building with cfitsio support: ..................................YES"
else
       echo "Building  with cfitsio support: ..................................NO"
       echo "        Requires cfitsio - A FITS File Subroutine Library"
fi
if [ test "x$have_nlopt_support" != x ]; then
       echo "Building with nlopt support: ....................................YES"
else
       echo "Building  with nlopt support: ....................................NO"
       echo "        Requires nlopt - NonLinear OPTimization library"
fi
if [ test "x$have_cuba_support" != x ]; then
       echo "Building with cuba support: .....................................YES"
else
       echo "Building  with libcuba support: ..................................NO"
       echo "        Requires libcuba - "
       echo "                 library for multidimensional numerical integration"
fi
echo "--------------------------------------------------------------------"
echo "AM_CFLAGS = $AM_CFLAGS"
echo "CFLAGS = $CFLAGS"
echo "LIBS = $LIBS"
