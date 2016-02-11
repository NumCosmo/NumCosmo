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

noinst_LTLIBRARIES = liblevmar.la

liblevmar_la_LIBADD = \
	$(SUNDIALS_LIBS) \
	$(FFTW3_LIBS) \
	$(GLIB_LIBS) \
	$(NLOPT_LIBS) \
	$(GSL_LIBS) \
	$(CFITSIO_LIBS)

misc.c: misc_core.c
lmbc.c: lmbc_core.c
lmblec.c: lmblec_core.c
lmbleic.c: lmbleic_core.c
lm.c: lm_core.c
lmlec.c: lmlec_core.c

liblevmar_la_SOURCES = \
	Axb.c          \
	compiler.h     \
	levmar.h       \
	lm.h           \
	misc.c         \
	misc.h         \
	levmar.h       \
	lmbc.c         \
	lmblec.c       \
	lmbleic.c      \
	lm.c           \
	lmlec.c

EXTRA_DIST =           \
	LICENSE        \
	Axb_core.c     \
	lmbc_core.c    \
	lmblec_core.c  \
	lmbleic_core.c \
	lm_core.c      \
	lmlec_core.c   \
	misc_core.c

CLEANFILES = 

