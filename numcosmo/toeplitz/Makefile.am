## Process this file with automake to produce Makefile.in

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS) \
	 $(SUNDIALS_CFLAGS) \
	 $(GSL_CFLAGS) \
	 $(FFTW3_CFLAGS) \
	 $(CFITSIO_CFLAGS) \
	 $(NLOPT_CFLAGS) \
	 -I$(top_srcdir) \
	 -I$(builddir)/..

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

SUBDIRS = .

noinst_LTLIBRARIES = libtoeplitz.la

libtoeplitz_la_LIBADD = \
	$(SUNDIALS_LIBS) \
	$(FFTW3_LIBS) \
	$(GLIB_LIBS) \
	$(NLOPT_LIBS) \
	$(GSL_LIBS) \
	$(CFITSIO_LIBS)

libtoeplitz_la_SOURCES = \
	solvers/dcommonsub.c  \
	solvers/dlevcl.c      \
	solvers/dschur.c      \
	solvers/dschurlev.c   \
	solvers/polytool.c    \
	solvers/dlevin.c      \
	solvers/drefine.c     \
	solvers/dschurin.c    \
	solvers/dschurlevin.c \
	solvers/toetool.c     \
	solvers/dlev.c        \
	solvers/dopt_k.c      \
	solvers/dschurinup.c  \
	solvers/dschurup.c    \
	solvers/dcommonsub.h  \
	solvers/dopt_k.h      \
	solvers/drefine.h     \
	solvers/fblas.h       \
	solvers/ftypes.h      \
	solvers/polytool.h    \
	solvers/toeplitz.h    \
	solvers/toetool.h

EXTRA_DIST = \
	COPYING \
	README

CLEANFILES = 

