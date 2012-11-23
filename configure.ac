dnl Process this file with autoconf to produce a configure script.
dnl Created by Anjuta application wizard.
dnl testa

m4_define([numcosmo_major_version], [0])
m4_define([numcosmo_minor_version], [8])
m4_define([numcosmo_micro_version], [0])
m4_define([numcosmo_version],
          [numcosmo_major_version.numcosmo_minor_version.numcosmo_micro_version])
          
AC_INIT([numcosmo],[numcosmo_version], https://savannah.nongnu.org/bugs/?func=additem&group=numcosmo, numcosmo, http://www.nongnu.org/numcosmo/)
AC_CONFIG_SRCDIR([numcosmo/numcosmo.h])
AC_CONFIG_MACRO_DIR([m4])

AC_PREREQ([2.69])
AM_INIT_AUTOMAKE([gnu -Wno-portability])

PACKAGE=numcosmo
NUMCOSMO_MAJOR_VERSION=numcosmo_major_version
NUMCOSMO_MINOR_VERSION=numcosmo_minor_version
NUMCOSMO_MICRO_VERSION=numcosmo_micro_version
NUMCOSMO_VERSION=$NUMCOSMO_MAJOR_VERSION.$NUMCOSMO_MINOR_VERSION.$NUMCOSMO_MICRO_VERSION
NUMCOSMO_LT_VERSION="$NUMCOSMO_MAJOR_VERSION:$NUMCOSMO_MINOR_VERSION:$NUMCOSMO_MICRO_VERSION"

AC_SUBST(NUMCOSMO_MAJOR_VERSION)
AC_SUBST(NUMCOSMO_MINOR_VERSION)
AC_SUBST(NUMCOSMO_MICRO_VERSION)
AC_SUBST(NUMCOSMO_VERSION)
AC_SUBST(NUMCOSMO_LT_VERSION)

AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE

dnl
dnl Gobject Introspection
dnl
dnl add the correct m4 directory location.
AC_CONFIG_MACRO_DIR([m4])
  
dnl any typos here, and you will end up with a message about HAVE_INTROSPECTION not being defined
GOBJECT_INTROSPECTION_CHECK([0.6.7])

dnl things required by automake
dnl AC_ARG_PROGRAM
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_LANG(C)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL   
AC_PROG_LN_S
AC_PROG_CC_STDC
dnl LT_INIT([disable-shared])
dnl AC_DISABLE_SHARED
LT_INIT

dnl check for gtk-doc
GTK_DOC_CHECK([1.14],[--flavour no-tmpl])

AC_ARG_ENABLE(man,
              [AC_HELP_STRING([--enable-man],
                              [regenerate man pages from Docbook [default=no]])],enable_man=yes,
              enable_man=no)

if test "${enable_man}" != no; then
dnl
dnl Check for xsltproc
dnl
AC_PATH_PROG([XSLTPROC], [xsltproc])
  if test -z "$XSLTPROC"; then
    enable_man=no
  fi
fi

AM_CONDITIONAL(ENABLE_MAN, test x$enable_man != xno)

dnl Set PACKAGE_DATA_DIR in config.h.
if test "x${prefix}" = "xNONE"; then
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${ac_default_prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
else
  AC_DEFINE_UNQUOTED(PACKAGE_DATA_DIR, "${prefix}/share/${PACKAGE}-${VERSION}", [PACKAGE_DATA_DIR])
fi

dnl Set PACKAGE_SOURCE_DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}",[PACKAGE_SOURCE_DIR])

dnl Dependencies
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

if test "$enable_inline" = "yes"; then
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

AC_ARG_WITH([thread-pool-max],
            AS_HELP_STRING([--with-thread-pool-max=N],[Number of threads in the pool]), 
            [thread_pool_max='yes'], [thread_pool_max='no'])

if test "x$thread_pool_max" = "xyes"; then
  AC_DEFINE_UNQUOTED([NCM_THREAD_POOL_MAX],[$with_thread_pool_max], [Number of threads in the pool])
fi

dnl ***************************************************************************
dnl Check for math.h
dnl ***************************************************************************

AC_CHECK_HEADERS(math.h,, AC_MSG_ERROR([math.h header files missing]))
AC_CHECK_LIB([m],[cos])

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

AC_CHECK_LIB([gmp], [__gmpz_init],,[AC_MSG_ERROR([GMP not found. See README.])],[$GLIB_LIBS $GSL_LIBS])
AC_CHECK_LIB([mpfr], [mpfr_init],,[AC_MSG_ERROR([MPFR not found. See README.])],[$GLIB_LIBS $GSL_LIBS])

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

PKG_CHECK_MODULES(FFTW3, [fftw3 >= $FFTW3_REQUIRED_VERSION],[
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
dnl Check for atlas + lapack
dnl ***************************************************************************

AC_MSG_CHECKING([--with-atlas-path])
AC_ARG_WITH(atlas-path,
  [AS_HELP_STRING([--with-atlas-path=PATH],[Full Path to Atlas.])],
  [ATLAS_PATH="-L$withval"])
AC_MSG_RESULT($ATLAS_PATH)

AC_MSG_CHECKING([--with-atlas])
AC_ARG_WITH(atlas,
  [AS_HELP_STRING([--with-atlas=atlas],[Atlas library name, default is 'atlas'])],
  [ATLAS="$withval"], [ATLAS="atlas"])
AC_MSG_RESULT($ATLAS)

AC_MSG_CHECKING([--with-atlas-cblas])
AC_ARG_WITH(atlas-cblas,
  [AS_HELP_STRING([--with-atlas-cblas=cblas],[Atlas cblas library name, default is 'cblas'])],
  [ATLAS_CBLAS="$withval"], [ATLAS_CBLAS="cblas"])
AC_MSG_RESULT($ATLAS_CBLAS)

AC_MSG_CHECKING([--with-atlas-lapack])
AC_ARG_WITH(atlas-lapack,
  [AS_HELP_STRING([--with-atlas-lapack=lapack],[Atlas lapack library name, default is 'lapack'])],
  [ATLAS_LAPACK="$withval"], [ATLAS_LAPACK="lapack"])
AC_MSG_RESULT($ATLAS_LAPACK)

AC_MSG_CHECKING([--with-atlas-extra-libs])
AC_ARG_WITH(atlas-extra-libs,
  [AS_HELP_STRING([--with-atlas-extra-libs=LIBS],[Atlas extra libraries, default is empty])],
  [ATLAS_EXTRA_LIBS="$withval"], [ATLAS_EXTRA_LIBS=""])
AC_MSG_RESULT($ATLAS_EXTRA_LIBS)

AC_CHECK_LIB([$ATLAS_CBLAS], [cblas_dgemv], [have_atlas_cblas='yes'],[],[$ATLAS_PATH -l$ATLAS $ATLAS_EXTRA_LIBS])
AC_CHECK_LIB([$ATLAS_LAPACK], [sgesv_], [have_lapack='yes'],[],[$ATLAS_PATH -l$ATLAS_CBLAS -l$ATLAS $ATLAS_EXTRA_LIBS])

if test "x$have_atlas_cblas" = "xyes"; then
  GSL_LIBS=`echo $GSL_LIBS | sed 's/\-lgslcblas//g'`
  APPEND_LIBS="$APPEND_LIBS $ATLAS_PATH -l$ATLAS_CBLAS -l$ATLAS $ATLAS_EXTRA_LIBS"
  have_atlas_support="#define NUMCOSMO_HAVE_ATLAS 1"
else
  have_atlas_support=""
fi
AC_SUBST(have_atlas_support)

if test "x$have_lapack" = "xyes"; then
  GSL_LIBS=`echo $GSL_LIBS | sed 's/\-lgslcblas//g'`
  APPEND_LIBS="$APPEND_LIBS -l$ATLAS_LAPACK"
  have_lapack_support="#define NUMCOSMO_HAVE_LAPACK 1"
  AC_DEFINE([HAVE_LAPACK], [1], [use lapack functions])
  AC_SUBST(HAVE_LAPACK)  
else
  have_lapack_support=""
fi
AC_SUBST(have_lapack_support)

dnl ***************************************************************************
dnl Check for levmar
dnl ***************************************************************************

AC_MSG_CHECKING([--with-levmar-path])
AC_ARG_WITH(levmar-path,
  [AS_HELP_STRING([--with-levmar-path=PATH],[Levmar library location])],
  [LEVMAR_PATH="-L$withval"], [LEVMAR_PATH=""])
AC_MSG_RESULT($LEVMAR_PATH)

AC_MSG_CHECKING([--with-levmar-extra-libs])
AC_ARG_WITH(levmar-extra-libs,
  [AS_HELP_STRING([--with-levmar-extra-libs=LIBS],[Levmar extra libraries, default is '-lf77blas -lgfortran'])],
  [LEVMAR_EXTRA_LIBS="$withval"], [LEVMAR_EXTRA_LIBS="-lf77blas -lgfortran"])
AC_MSG_RESULT($LEVMAR_EXTRA_LIBS)

AC_CHECK_LIB([levmar], [dlevmar_der],,[levmar_missing="yes"],[$GLIB_LIBS $GSL_LIBS $APPEND_LIBS $LEVMAR_PATH $LEVMAR_EXTRA_LIBS])
AC_CHECK_HEADER(levmar.h, , [AC_CHECK_HEADER(levmar/levmar.h, [levmar_need_prefix="yes"], [levmar_missing="yes"])])
if test "x$levmar_missing" = "xyes"; then
  have_levmar_support=""
else
  have_levmar_support="#define NUMCOSMO_HAVE_LEVMAR 1"
  APPEND_LIBS="$APPEND_LIBS $LEVMAR_PATH $LEVMAR_EXTRA_LIBS"
fi
AC_SUBST(have_levmar_support)

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
dnl AC_CHECK_LIB([chealpix], [read_healpix_map],,[chealpix_missing="yes"],[$GLIB_LIBS $APPEND_LIBS])
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

AC_CHECK_LIB([cuba], [Cuhre],,[cuba_missing="yes"],[$GLIB_LIBS $APPEND_LIBS])
if test "x$cuba_missing" = "xyes"; then
  have_cuba_support=""
else
  have_cuba_support="#define NUMCOSMO_HAVE_CUBA 1"
fi
AC_SUBST(have_cuba_support)

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

PKG_CHECK_EXISTS([nlopt >= 2.2], [
AC_DEFINE([HAVE_NLOPT_2_2],[1], [If nlopt version is bigger than 2.2])
AC_CHECK_LIB([nlopt],[nlopt_create],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with c++ support, recompile the C version.])],[$APPEND_LIBS])
],[
AC_CHECK_LIB([nlopt],[nlopt_minimize_constrained],,[AC_MSG_ERROR([Can't link to -lnlopt. The nlopt library is probably compiled with C++ support, recompile the C version.])],[$APPEND_LIBS])
])

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
AC_CHECK_TYPE([CVDlsDenseJacFn],[],[AC_MSG_ERROR([sundials too old, please install version >= 2.4.0])], [#include <cvodes/cvodes_direct.h>])
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

LIBS="$LIBS $APPEND_LIBS"
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
if [ test "x$have_atlas_support" != x ]; then
       echo "Building with Atlas support: ....................................YES"
else
       echo "Building  with Atlas support: ....................................NO"
       echo "        Requires ATLAS - Automatically Tuned Linear Algebra Software"
fi
if [ test "x$have_lapack_support" != x ]; then
       echo "Building with Lapack support: ...................................YES"
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
