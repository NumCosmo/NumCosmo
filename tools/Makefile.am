
INCLUDES = 

AM_CFLAGS = \
	-Wall \
	$(GLIB_CFLAGS) \
	$(SUNDIALS_CFLAGS) \
	$(GSL_CFLAGS) \
	-I$(top_srcdir) \
	-I$(top_srcdir)/numcosmo

noinst_PROGRAMS = \
	cmb_maps \
	q-piecewise \
	q-linear \
	q-constant \
	recomb \
	quantum-gravity

quantum_gravity_SOURCES = \
	quantum-gravity.c

quantum_gravity_LDADD = \
        $(top_builddir)/numcosmo/libnumcosmo.la

recomb_SOURCES = \
	recomb.c

recomb_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

q_constant_SOURCES = \
	q-constant.c

q_constant_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

q_linear_SOURCES = \
	q-linear.c

q_linear_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

q_piecewise_SOURCES = \
	q-piecewise.c

q_piecewise_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

cmb_maps_SOURCES = \
	cmb_maps.c

cmb_maps_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

AM_CPPFLAGS = 

bin_PROGRAMS = \
	gobj_itest

gobj_itest_SOURCES = \
	gobj_itest.c

gobj_itest_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la

## File created by the gnome-build tools

