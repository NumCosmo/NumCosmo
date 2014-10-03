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

noinst_PROGRAMS = \
	cmb_maps \
	q-piecewise \
	q-linear \
	q-constant \
	recomb \
	quantum-gravity \
	gobj_itest

quantum_gravity_SOURCES = \
	quantum-gravity.c

quantum_gravity_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la \
        $(GLIB_LIBS)

recomb_SOURCES = \
	recomb.c

recomb_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

q_constant_SOURCES = \
	q-constant.c

q_constant_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

q_linear_SOURCES = \
	q-linear.c

q_linear_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

q_piecewise_SOURCES = \
	q-piecewise.c

q_piecewise_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS)

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
