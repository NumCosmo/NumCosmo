## Process this file with automake to produce Makefile.in

AM_CFLAGS = \
	-Wall \
	$(GLIB_CFLAGS) \
	$(SUNDIALS_CFLAGS) \
	$(GSL_CFLAGS) \
	-I$(top_srcdir) \
	-I$(top_srcdir)/numcosmo

bin_PROGRAMS = darkenergy

darkenergy_SOURCES = \
	darkenergy.c\
	data_cluster.c \
	data_cluster.h \
	de_options.c \
	de_options.h \
	savedata.c \
	savedata.h

darkenergy_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la

