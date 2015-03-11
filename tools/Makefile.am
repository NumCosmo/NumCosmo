## Process this file with automake to produce Makefile.in

MY_CFLAGS = \
         -DG_LOG_DOMAIN=\"NUMCOSMO\" \
         $(GLIB_CFLAGS) \
         $(SUNDIALS_CFLAGS) \
         $(GSL_CFLAGS) \
         $(FFTW3_CFLAGS) \
         $(CFITSIO_CFLAGS) \
         $(NLOPT_CFLAGS) \
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

AM_CPPFLAGS = 

bin_PROGRAMS = 

## File created by the gnome-build tools
