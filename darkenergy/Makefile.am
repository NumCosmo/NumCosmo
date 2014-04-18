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

EXTRA_DIST = 

AM_VALAFLAGS = 

if HAVE_GTK

AM_CFLAGS += $(GTK_CFLAGS) -Wl,--export-dynamic

AM_VALAFLAGS += \
	--pkg numcosmomath-1.0 \
	--pkg numcosmo-1.0     \
	--pkg gtk+-3.0         \
	--pkg gmodule-2.0      \
	--vapidir=$(abs_top_builddir)/numcosmo

bin_PROGRAMS += gdarkenergy

gdarkenergy_SOURCES =  \
	gdarkenergy.vala

gdarkenergy_DATA =  \
	gdarkenergy.glade

EXTRA_DIST += $(gdarkenergy_DATA)

gdarkenergydir = \
        $(pkgdatadir)-$(VERSION)/gdarkenergy

gdarkenergy_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la $(GTK_LIBS)

endif
