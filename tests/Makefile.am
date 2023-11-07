## Process this file with automake to produce Makefile.in

# gcov support
include $(top_srcdir)/aminclude_static.am

include $(top_srcdir)/glib-tap.mk

TEST_EXTENSIONS = .py

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS)     \
	 $(SUNDIALS_CFLAGS) \
	 $(GSL_CFLAGS)      \
	 $(FFTW3_CFLAGS)    \
	 $(CFITSIO_CFLAGS)  \
	 $(NLOPT_CFLAGS)    \
	 $(OPENMP_CFLAGS)   \
	 -I$(top_srcdir)    \
	 $(CODE_COVERAGE_CFLAGS)

AM_CFLAGS = \
	 -Wall \
	 $(MY_CFLAGS)

COVLIBS = \
	$(CODE_COVERAGE_LIBS)

test_ncm_cfg_SOURCES =  \
	test_ncm_cfg.c

test_ncm_sparam_SOURCES =  \
	test_ncm_sparam.c

test_ncm_diff_SOURCES =  \
	test_ncm_diff.c

test_ncm_ode_SOURCES =  \
	test_ncm_ode.c

test_ncm_fftlog_SOURCES =  \
	test_ncm_fftlog.c

test_ncm_vector_SOURCES =  \
	test_ncm_vector.c

test_ncm_matrix_SOURCES =  \
	test_ncm_matrix.c

test_ncm_stats_vec_SOURCES =  \
	test_ncm_stats_vec.c
	
test_ncm_stats_dist1d_epdf_SOURCES =  \
	test_ncm_stats_dist1d_epdf.c
	
test_ncm_stats_dist_SOURCES =  \
	test_ncm_stats_dist.c

test_ncm_stats_dist_kernel_SOURCES = \
	test_ncm_stats_dist_kernel.c
	
test_ncm_spline_SOURCES =  \
	test_ncm_spline.c

test_ncm_spline2d_SOURCES =  \
	test_ncm_spline2d.c

test_ncm_integral1d_SOURCES =  \
	test_ncm_integral1d.c

test_ncm_integral_nd_SOURCES =  \
	test_ncm_integral_nd.c

test_ncm_sf_sbessel_SOURCES =  \
	test_ncm_sf_sbessel.c

test_ncm_mpsf_0F1_SOURCES =  \
	test_ncm_mpsf_0F1.c

test_ncm_mpsf_trig_int_SOURCES =  \
	test_ncm_mpsf_trig_int.c

test_ncm_model_SOURCES =  \
        test_ncm_model.c \
	ncm_model_test.c \
	ncm_model_test.h

test_ncm_model_ctrl_SOURCES =  \
        test_ncm_model_ctrl.c

test_ncm_mset_SOURCES = \
	test_ncm_mset.c

test_ncm_mset_catalog_SOURCES = \
	test_ncm_mset_catalog.c

test_ncm_obj_array_SOURCES = \
	test_ncm_obj_array.c

test_ncm_serialize_SOURCES = \
	test_ncm_serialize.c

test_ncm_data_gauss_cov_SOURCES = \
	test_ncm_data_gauss_cov.c \
	ncm_data_gauss_cov_test.c \
	ncm_data_gauss_cov_test.h

test_ncm_fit_SOURCES = \
	test_ncm_fit.c

test_ncm_fit_state_SOURCES = \
	test_ncm_fit_state.c

test_ncm_fit_esmcmc_SOURCES = \
	test_ncm_fit_esmcmc.c

test_ncm_func_eval_SOURCES =  \
	test_ncm_func_eval.c

test_ncm_sf_spherical_harmonics_SOURCES =  \
	test_ncm_sf_spherical_harmonics.c

test_ncm_sphere_map_SOURCES =  \
	test_ncm_sphere_map.c
	
test_nc_halo_bias_SOURCES =  \
	test_nc_halo_bias.c

test_nc_hicosmo_de_SOURCES =  \
	test_nc_hicosmo_de.c

test_nc_window_SOURCES =  \
	test_nc_window.c
        
test_nc_transfer_func_SOURCES =  \
	test_nc_transfer_func.c        

test_nc_galaxy_acf_SOURCES =  \
	test_nc_galaxy_acf.c

test_nc_recomb_SOURCES =  \
	test_nc_recomb.c

test_nc_cbe_SOURCES =  \
	test_nc_cbe.c

test_nc_hiqg_1d_SOURCES =  \
	test_nc_hiqg_1d.c

test_nc_data_bao_rdv_SOURCES =  \
	test_nc_data_bao_rdv.c

test_nc_data_bao_dvdv_SOURCES =  \
	test_nc_data_bao_dvdv.c
        
test_nc_data_cluster_ncount_SOURCES =  \
	test_nc_data_cluster_ncount.c
        
test_nc_cluster_pseudo_counts_SOURCES =  \
	test_nc_cluster_pseudo_counts.c

test_nc_cluster_abundance_SOURCES =  \
	test_nc_cluster_abundance.c

test_nc_halo_density_profile_SOURCES =  \
	test_nc_halo_density_profile.c

test_nc_wl_surface_mass_density_SOURCES =  \
	test_nc_wl_surface_mass_density.c

test_nc_powspec_SOURCES =  \
	test_nc_powspec.c

test_nc_distance_SOURCES =  \
	test_nc_distance.c

test_nc_galaxy_wl_ellipticity_kde_SOURCES =  \
	test_nc_galaxy_wl_ellipticity_kde.c

test_nc_galaxy_wl_ellipticity_binned_SOURCES =  \
	test_nc_galaxy_wl_ellipticity_binned.c

test_programs =  \
	test_ncm_cfg                             \
	test_ncm_vector                          \
	test_ncm_matrix                          \
	test_ncm_stats_vec                       \
	test_ncm_stats_dist1d_epdf               \
	test_ncm_stats_dist                      \
	test_ncm_stats_dist_kernel               \
	test_ncm_spline                          \
	test_ncm_spline2d                        \
	test_ncm_integral1d                      \
	test_ncm_integral_nd                     \
	test_ncm_sf_sbessel                      \
	test_ncm_mpsf_0F1                        \
	test_ncm_mpsf_trig_int                   \
	test_ncm_func_eval                       \
	test_ncm_sparam                          \
	test_ncm_diff                            \
	test_ncm_ode                             \
	test_ncm_fftlog                          \
	test_ncm_model                           \
	test_ncm_model_ctrl                      \
	test_ncm_serialize                       \
	test_ncm_mset                            \
	test_ncm_mset_catalog                    \
	test_ncm_obj_array                       \
	test_ncm_data_gauss_cov                  \
	test_ncm_fit                             \
	test_ncm_fit_state                       \
	test_ncm_fit_esmcmc                      \
	test_ncm_sf_spherical_harmonics          \
	test_ncm_sphere_map                      \
	test_nc_halo_bias                        \
	test_nc_hiqg_1d                          \
	test_nc_hicosmo_de                       \
	test_nc_window                           \
	test_nc_transfer_func                    \
	test_nc_galaxy_acf                       \
	test_nc_recomb                           \
	test_nc_cbe                              \
	test_nc_data_bao_rdv                     \
	test_nc_data_bao_dvdv                    \
	test_nc_data_cluster_ncount              \
	test_nc_cluster_pseudo_counts            \
	test_nc_cluster_abundance                \
	test_nc_halo_density_profile             \
	test_nc_wl_surface_mass_density          \
	test_nc_powspec                          \
	test_nc_distance                         \
	test_nc_galaxy_wl_ellipticity_kde        \
	test_nc_galaxy_wl_ellipticity_binned

if HAVE_PYTEST_TAP

PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = -m pytest --tap-stream
AM_TESTS_ENVIRONMENT += \
	GI_TYPELIB_PATH="$(abs_top_builddir)/numcosmo$${GI_TYPELIB_PATH:+:$${GI_TYPELIB_PATH}}" \
	LD_LIBRARY_PATH="$(abs_top_builddir)/numcosmo/.libs$${LD_LIBRARY_PATH:+:$${LD_LIBRARY_PATH}}" \
	PYTHONPATH="$(abs_top_srcdir)$${PYTHONPATH:+:$${PYTHONPATH}}"

dist_test_scripts = \
	test_py_sanity.py                     \
	test_py_integralnd.py                 \
	test_py_diff.py                       \
	test_py_cluster_abundance_sampling.py

endif # HAVE_PYTEST_TAP

test_ncm_cfg_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_sparam_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_diff_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_ode_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_fftlog_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_vector_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_matrix_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_stats_vec_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_stats_dist1d_epdf_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_stats_dist_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)
	
test_ncm_stats_dist_kernel_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)
	
test_ncm_spline_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_spline2d_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_integral1d_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_integral_nd_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_sf_sbessel_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_mpsf_0F1_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_mpsf_trig_int_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_model_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_model_ctrl_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_mset_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_mset_catalog_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_obj_array_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_serialize_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_data_gauss_cov_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_fit_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_fit_state_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_fit_esmcmc_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_func_eval_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_sf_spherical_harmonics_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_ncm_sphere_map_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)
	
test_nc_halo_bias_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_hicosmo_de_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_window_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_transfer_func_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_galaxy_acf_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_recomb_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_cbe_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_hiqg_1d_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_data_bao_rdv_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_data_bao_dvdv_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_data_cluster_ncount_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_cluster_pseudo_counts_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_cluster_abundance_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_halo_density_profile_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)	

test_nc_powspec_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)	

test_nc_wl_surface_mass_density_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_distance_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_galaxy_wl_ellipticity_kde_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

test_nc_galaxy_wl_ellipticity_binned_LDADD = \
	$(top_builddir)/numcosmo/libnumcosmo.la \
	$(GLIB_LIBS) \
	$(GSL_LIBS) \
	$(COVLIBS)

EXTRA_DIST = \
        test_ncm_template.c

