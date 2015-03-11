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
	 -I$(srcdir)/../common

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

SUBDIRS = .

noinst_LTLIBRARIES = libcuba_cuhre.la

libcuba_cuhre_la_SOURCES = \
	decl.h       \
	Cuhre.c

EXTRA_DIST = \
	common.c    \
	Cuhre.c     \
	decl.h      \
	Integrate.c \
	Rule.c

Cuhre.c: \
	common.c    \
	Integrate.c \
	Rule.c      \
	../common/stddecl.h   \
	../common/ChiSquare.c \
	../common/Erf.c       \
	../common/CSample.c   \
	../common/Parallel.c  \
	../common/sock.h

CLEANFILES = 
