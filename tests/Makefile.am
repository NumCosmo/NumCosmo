## Process this file with automake to produce Makefile.in
include $(top_srcdir)/Makefile.decl

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

test_ncm_sf_sbessel_SOURCES =  \
	test_ncm_sf_sbessel.c

test_ncm_model_SOURCES =  \
        test_ncm_model.c \
	ncm_model_test.c \
	ncm_model_test.h

test_ncm_object_serialization_SOURCES = \
	test_ncm_object_serialization.c

test_ncm_data_gauss_cov_SOURCES = \
	test_ncm_data_gauss_cov.c \
	ncm_data_gauss_cov_test.c \
	ncm_data_gauss_cov_test.h

test_nc_window_SOURCES =  \
        test_nc_window.c
        
test_nc_transfer_func_SOURCES =  \
        test_nc_transfer_func.c        

test_nc_galaxy_acf_SOURCES =  \
	test_nc_galaxy_acf.c

test_nc_recomb_SOURCES =  \
	test_nc_recomb.c

check_PROGRAMS =  \
	test_ncm_sparam               \
	test_ncm_vector               \
	test_ncm_matrix               \
	test_ncm_spline               \
	test_ncm_spline2d             \
	test_ncm_sf_sbessel           \
	test_ncm_model                \
	test_ncm_object_serialization \
	test_ncm_data_gauss_cov       \
	test_nc_window                \
	test_nc_transfer_func         \
	test_nc_galaxy_acf            \
	test_nc_recomb

# TEST_PROGS += $(check_PROGRAMS)

test_ncm_sparam_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_vector_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_matrix_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_ncm_spline_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_spline2d_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_sf_sbessel_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_model_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_object_serialization_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_ncm_data_gauss_cov_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_window_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_transfer_func_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la

test_nc_galaxy_acf_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

test_nc_recomb_LDADD = $(top_builddir)/numcosmo/libnumcosmo.la 

TESTS = $(check_PROGRAMS)
