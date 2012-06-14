## Process this file with automake to produce Makefile.in

AM_CFLAGS = -Wall \
	$(GLIB_CFLAGS) \
	$(GSL_CFLAGS) \
	$(SUNDIALS_CFLAGS) \
	$(NUMCOSMO_CFLAGS) \
	-I$(top_srcdir) 

test_ncm_sparam_SOURCES =  \
	test_ncm_sparam.c

test_ncm_vector_SOURCES =  \
	test_ncm_vector.c

test_ncm_matrix_SOURCES =  \
	test_ncm_matrix.c
	
test_ncm_spline_SOURCES =  \
        test_ncm_spline.c

test_ncm_spline2d_SOURCES =  \
        test_ncm_spline2d.c

test_ncm_model_SOURCES =  \
        test_ncm_model.c \
	ncm_model_test.c \
	ncm_model_test.h
        
test_nc_window_SOURCES =  \
        test_nc_window.c
        
test_nc_transfer_func_SOURCES =  \
        test_nc_transfer_func.c        

test_nc_galaxy_acf_SOURCES =  \
	test_nc_galaxy_acf.c

check_PROGRAMS =  \
	test_ncm_sparam    \
	test_ncm_vector    \
	test_ncm_matrix    \
	test_ncm_spline    \
	test_ncm_spline2d  \
	test_ncm_model     \
	test_nc_window     \
	test_nc_transfer_func \
	test_nc_galaxy_acf

test_ncm_sparam_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_vector_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_matrix_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_spline_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_spline2d_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_model_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_window_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_transfer_func_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_galaxy_acf_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

TESTS = $(check_PROGRAMS)

