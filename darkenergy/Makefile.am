## Process this file with automake to produce Makefile.in

AM_CFLAGS = \
	-Wall \
	$(GLIB_CFLAGS) \
	$(SUNDIALS_CFLAGS) \
	$(GSL_CFLAGS) \
	-I$(top_srcdir) \
	-I$(top_srcdir)/numcosmo

AM_VALAFLAGS = \
	--pkg numcosmomath-1.0 \
	--pkg numcosmo-1.0     \
	--vapidir=$(abs_top_builddir)/numcosmo

bin_PROGRAMS = darkenergy \
	gdarkenergy

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

gdarkenergy_SOURCES =  \
	gdarkenergy.vala

gdarkenergy_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la
