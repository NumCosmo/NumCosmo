## Process this file with automake to produce Makefile.in

MY_CFLAGS = \
         -DG_LOG_DOMAIN=\"NUMCOSMO\" \
         $(GLIB_CFLAGS)     \
         $(SUNDIALS_CFLAGS) \
         $(GSL_CFLAGS)      \
         $(FFTW3_CFLAGS)    \
         $(FFTW3F_CFLAGS)   \
         $(CFITSIO_CFLAGS)  \
         $(NLOPT_CFLAGS)    \
         -I$(top_srcdir)

AM_CFLAGS = \
	-Wall \
	$(MY_CFLAGS)

noinst_PROGRAMS =  \
	cmb_maps   \
	gobj_itest

cmb_maps_SOURCES = \
	cmb_maps.c

cmb_maps_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

gobj_itest_SOURCES = \
	gobj_itest.c

gobj_itest_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

mcat_analyze_SOURCES = \
	mcat_analyze.c

mcat_analyze_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la  \
	$(GLIB_LIBS)

mset_gen_SOURCES =  \
	mset_gen.c

mset_gen_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

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

AM_CPPFLAGS = 

bin_PROGRAMS = \
	darkenergy   \
	mcat_analyze \
	mset_gen

