## Process this file with automake to produce Makefile.in

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS)              \
	 $(SUNDIALS_CFLAGS)          \
	 $(GSL_CFLAGS)               \
	 $(FFTW3_CFLAGS)             \
	 $(CFITSIO_CFLAGS)           \
	 $(NLOPT_CFLAGS)             \
	 -I$(top_srcdir)             \
	 $(VISIBILITY_HIDDEN_FLAG)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

AM_LDFLAGS = $(VISIBILITY_HIDDEN_FLAG)

SUBDIRS = cuhre divonne suave vegas

noinst_LTLIBRARIES = libcuba.la

libcuba_la_LIBADD =                \
	cuhre/libcuba_cuhre.la     \
	divonne/libcuba_divonne.la \
	suave/libcuba_suave.la     \
	vegas/libcuba_vegas.la

libcuba_la_SOURCES =    \
	cuba.h          \
	common/Fork.c   \
	common/Global.c \
	common/Data.c

EXTRA_DIST = \
	COPYING \
	common/ChiSquare.c \
	common/CSample.c   \
	common/Erf.c       \
	common/MSample.c   \
	common/Parallel.c  \
	common/Random.c    \
	common/sock.h      \
	common/stddecl.h

CLEANFILES = 

