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

noinst_LTLIBRARIES = libcuba_divonne.la

libcuba_divonne_la_LIBADD = 

libcuba_divonne_la_SOURCES = \
	decl.h       \
	Divonne.c

EXTRA_DIST =                 \
	FindMinimum.c        \
	Sample.c             \
	common.c             \
	Rule.c               \
	Split.c              \
	Explore.c            \
	KorobovCoeff.c       \
	Iterate.c            \
	Integrate.c          \
	KorobovCoeff.c-5879  \
	KorobovCoeff.c-3319  \
	KorobovCoeff.c-9689  \
	KorobovCoeff.c-16033 \
	KorobovCoeff.c-2423

Divonne.c: \
	KorobovCoeff.c        \
	KorobovCoeff.c-5879   \
	KorobovCoeff.c-3319   \
	KorobovCoeff.c-9689   \
	KorobovCoeff.c-16033  \
	KorobovCoeff.c-2423   \
	Sample.c              \
	common.c              \
	Rule.c                \
	Split.c               \
	Explore.c             \
	Iterate.c             \
	Integrate.c           \
	../common/stddecl.h   \
	../common/ChiSquare.c \
	../common/Erf.c       \
	../common/CSample.c   \
	../common/Parallel.c  \
	../common/sock.h

CLEANFILES = 

