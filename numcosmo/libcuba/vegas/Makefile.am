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

noinst_LTLIBRARIES = libcuba_vegas.la

libcuba_vegas_la_LIBADD = 

libcuba_vegas_la_SOURCES = \
	decl.h             \
	Vegas.c

EXTRA_DIST =        \
	common.c    \
	Grid.c      \
	Integrate.c

Divonne.c: \
	common.c    \
	Grid.c      \
	Integrate.c \
	../common/stddecl.h   \
	../common/ChiSquare.c \
	../common/Erf.c       \
	../common/CSample.c   \
	../common/Parallel.c  \
	../common/sock.h

CLEANFILES = 

