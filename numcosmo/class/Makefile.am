## Process this file with automake to produce Makefile.in

CLASS_DIR = $(top_srcdir)/numcosmo/class

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS)     \
	 $(SUNDIALS_CFLAGS) \
	 $(GSL_CFLAGS)      \
	 $(FFTW3_CFLAGS)    \
	 $(CFITSIO_CFLAGS)  \
	 $(NLOPT_CFLAGS)    \
	 -I$(top_srcdir)    \
	 -I$(CLASS_DIR)     \
	 -I$(CLASS_DIR)/include \
	 -I$(CLASS_DIR)/hyrec   \
	 $(CLASS_DEFINES)       \
	 $(OPENMP_CFLAGS)       \
	 -DHYREC                \
	 -D__CLASSDIR__=\"$(PACKAGE_DATA_DIR)\/class_data\" \
	 $(VISIBILITY_HIDDEN_FLAG)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)
	 
AM_LDFLAGS = $(VISIBILITY_HIDDEN_FLAG)

SUBDIRS = .

noinst_LTLIBRARIES = libclass.la

libclass_la_LIBADD = \
	$(OPENMP_CFLAGS)

## CLik modules sources

class_tool_sources =           \
	tools/growTable.c      \
	tools/dei_rkck.c       \
	tools/sparse.c         \
	tools/evolver_rkck.c   \
	tools/evolver_ndf15.c  \
	tools/arrays.c         \
	tools/parser.c         \
	tools/quadrature.c     \
	tools/hyperspherical.c \
	tools/common.c         \
	tools/hermite3_interpolation_csource.h \
	tools/hermite4_interpolation_csource.h \
	tools/hermite6_interpolation_csource.h \
	include/common.h         \
	include/hyperspherical.h \
	include/evolver_rkck.h   \
	include/evolver_ndf15.h  \
	include/growTable.h      \
	include/quadrature.h     \
	include/parser.h         \
	include/dei_rkck.h       \
	include/sparse.h         \
	include/arrays.h

class_hyrec_sources = \
	hyrec/hyrectools.c \
	hyrec/helium.c     \
	hyrec/hydrogen.c   \
	hyrec/history.c    \
	hyrec/hyrec.h      \
	hyrec/history.h    \
	hyrec/helium.h     \
	hyrec/hydrogen.h   \
	hyrec/hyrectools.h

class_sources = \
	source/background.c      \
	source/thermodynamics.c  \
	source/perturbations.c   \
	source/primordial.c      \
	source/nonlinear.c       \
	source/transfer.c        \
	source/spectra.c         \
	source/lensing.c         \
	include/nonlinear.h      \
	include/background.h     \
	include/thermodynamics.h \
	include/primordial.h     \
	include/lensing.h        \
	include/perturbations.h  \
	include/spectra.h        \
	include/transfer.h       \
	include/class.h          \
	include/svnversion.h
	

libclass_la_SOURCES = \
	$(class_tool_sources)  \
	$(class_hyrec_sources) \
	$(class_sources)

EXTRA_DIST =

CLEANFILES =
