## Process this file with automake to produce Makefile.in

@CODE_COVERAGE_RULES@

if USE_INTERNAL_LIBCUBA
  MAYBE_LIBCUBA  = libcuba
  LIBCUBA_PLACE  = libcuba/libcuba.la
  LIBCUBA_INCDIR = -I$(top_srcdir)/numcosmo/libcuba
else

endif

PLC_DIR = $(top_srcdir)/numcosmo/plc
SUNDIALS_INCL = $(top_srcdir)/numcosmo/sundials/include

MY_CFLAGS = \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS)       \
	 $(SUNDIALS_CFLAGS)   \
	 $(GSL_CFLAGS)        \
	 $(FFTW3_CFLAGS)      \
	 $(FFTW3F_CFLAGS)     \
	 $(CFITSIO_CFLAGS)    \
	 $(NLOPT_CFLAGS)      \
	 $(LIBCUBA_INCDIR)    \
	 $(OPENMP_CFLAGS)     \
	 $(HDF5_CPPFLAGS)     \
	 $(HDF5_CFLAGS)       \
	 -I$(PLC_DIR)         \
	 -I$(PLC_DIR)/cldf    \
	 -I$(PLC_DIR)/minipmc \
	 -I$(PLC_DIR)/plik    \
	 -I$(top_srcdir)      \
	 -I$(SUNDIALS_INCL)   \
	 $(CODE_COVERAGE_CFLAGS)

AM_CFLAGS =    \
	 -Wall \
	 $(MY_CFLAGS)

if NUMCOSMO_CHECK_CCL
  AM_CFLAGS += $(CCL_CFLAGS)
endif

SUBDIRS = $(MAYBE_LIBCUBA) levmar plc class toeplitz sundials .

lib_LTLIBRARIES = libnumcosmo.la

libnumcosmo_la_LDFLAGS = \
	-version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE) $(NOUNDEFINED_FLAG)

libnumcosmo_la_LIBADD =         \
	$(SUNDIALS_LIBS)        \
	$(FFTW3_LIBS)           \
	$(FFTW3F_LIBS)          \
	$(GLIB_LIBS)            \
	$(NLOPT_LIBS)           \
	$(GSL_LIBS)             \
	$(CFITSIO_LIBS)         \
	$(LIBCUBA_PLACE)        \
	$(HDF5_LDFLAGS)         \
	$(HDF5_LIBS)            \
	$(CODE_COVERAGE_LIBS)   \
	plc/libclik.la          \
	levmar/liblevmar.la     \
	class/libclass.la       \
	sundials/libsundials.la \
	toeplitz/libtoeplitz.la

if NUMCOSMO_CHECK_CCL
  libnumcosmo_la_LIBADD += $(CCL_LDFLAGS)
endif

ncm_sources = \
	ncm_build_inline.c                   \
	math/ncm_c.c                         \
	math/ncm_cfg.c                       \
	math/ncm_mpi_job.c                   \
	math/ncm_mpi_job_test.c              \
	math/ncm_mpi_job_fit.c               \
	math/ncm_mpi_job_mcmc.c              \
	math/ncm_util.c                      \
	math/ncm_diff.c                      \
	math/ncm_ode.c                       \
	math/ncm_ode_eval.c                  \
	math/ncm_timer.c                     \
	math/ncm_lapack.c                    \
	math/ncm_vector.c                    \
	math/ncm_matrix.c                    \
	math/ncm_serialize.c                 \
	math/ncm_obj_array.c                 \
	math/ncm_integral1d.c                \
	math/ncm_integral1d_ptr.c            \
	math/ncm_rng.c                       \
	math/ncm_stats_vec.c                 \
	math/ncm_stats_dist1d.c              \
	math/ncm_stats_dist1d_spline.c       \
	math/ncm_stats_dist1d_epdf.c         \
	math/ncm_stats_dist2d.c              \
        math/ncm_stats_dist2d_spline.c       \
	math/ncm_stats_dist_nd.c             \
	math/ncm_stats_dist_nd_kde_gauss.c   \
	math/ncm_bootstrap.c                 \
	math/ncm_spline.c                    \
	math/ncm_spline_gsl.c                \
	math/ncm_spline_cubic.c              \
	math/ncm_spline_cubic_notaknot.c     \
	math/ncm_spline_rbf.c                \
	math/ncm_spline_func.c               \
	math/ncm_spline2d.c                  \
	math/ncm_spline2d_spline.c           \
	math/ncm_spline2d_gsl.c              \
	math/ncm_spline2d_bicubic.c          \
	math/ncm_powspec.c                   \
	math/ncm_powspec_filter.c            \
	math/ncm_powspec_corr3d.c            \
	math/ncm_hoaa.c                      \
	math/ncm_csq1d.c                     \
	math/ncm_func_eval.c                 \
	math/ncm_mpsf_trig_int.c             \
	math/ncm_mpsf_sbessel.c              \
	math/ncm_sf_sbessel.c                \
	math/ncm_sf_sbessel_int.c            \
	math/ncm_sf_spherical_harmonics.c    \
	math/ncm_mpsf_sbessel_int.c          \
	math/ncm_mpsf_0F1.c                  \
	math/ncm_fftlog.c                    \
	math/ncm_fftlog_sbessel_j.c          \
	math/ncm_fftlog_tophatwin2.c         \
	math/ncm_fftlog_gausswin2.c          \
	math/ncm_sparam.c                    \
	math/ncm_vparam.c                    \
	math/ncm_model.c                     \
	math/ncm_model_ctrl.c                \
	math/ncm_model_builder.c             \
	math/ncm_model_mvnd.c                \
	math/ncm_mset.c                      \
	math/ncm_mset_func.c                 \
	math/ncm_mset_func1.c                \
	math/ncm_mset_func_list.c            \
	math/ncm_calc.c                      \
	math/ncm_ode_spline.c                \
	math/ncm_reparam.c                   \
	math/ncm_reparam_linear.c            \
	math/ncm_data.c                      \
	math/ncm_data_dist1d.c               \
	math/ncm_data_dist2d.c               \
	math/ncm_data_gauss.c                \
	math/ncm_data_gauss_cov.c            \
	math/ncm_data_gauss_diag.c           \
	math/ncm_data_poisson.c              \
	math/ncm_data_gauss_cov_mvnd.c       \
	math/ncm_dataset.c                   \
	math/ncm_likelihood.c                \
	math/ncm_prior.c                     \
	math/ncm_prior_gauss.c               \
	math/ncm_prior_gauss_param.c         \
	math/ncm_prior_gauss_func.c          \
	math/ncm_prior_flat.c                \
	math/ncm_prior_flat_param.c          \
	math/ncm_prior_flat_func.c           \
	math/ncm_mset_trans_kern.c           \
	math/ncm_mset_trans_kern_flat.c      \
	math/ncm_mset_trans_kern_gauss.c     \
	math/ncm_mset_trans_kern_cat.c       \
	math/ncm_fit_state.c                 \
	math/ncm_fit.c                       \
	math/ncm_fit_levmar.c                \
	math/ncm_fit_gsl_ls.c                \
	math/ncm_fit_gsl_mm.c                \
	math/ncm_fit_gsl_mms.c               \
	math/ncm_mset_catalog.c              \
	math/ncm_fit_mc.c                    \
	math/ncm_fit_mcbs.c                  \
	math/ncm_fit_mcmc.c                  \
	math/ncm_fit_esmcmc.c                \
	math/ncm_fit_esmcmc_walker.c         \
	math/ncm_fit_esmcmc_walker_stretch.c \
	math/ncm_fit_esmcmc_walker_walk.c    \
	math/ncm_fit_esmcmc_walker_aps.c     \
	math/ncm_lh_ratio1d.c                \
	math/ncm_lh_ratio2d.c                \
	math/ncm_abc.c                       \
	math/ncm_quaternion.c                \
	math/ncm_sphere_map.c                \
	math/ncm_function_cache.c            \
	math/poly.c               \
	math/integral.c           \
	math/quadrature.c         \
	math/matrix_exp.c         \
	math/magnus_iserles_ode.c \
	math/binsplit.c           \
	math/mpq_tree.c           \
	math/grid_one.c           \
	math/dividedifference.c   \
	math/ncm_memory_pool.c

ncm_sources_nogi = \
	math/Faddeeva.c  \
	math/rcm.c       \
	math/LowRankQP.c \
	gslextras/cqp/cqpminimizer.c  \
	gslextras/cqp/initial_point.c \
	gslextras/cqp/mg_pdip.c

ncm_headers =                                \
	math/ncm_c.h                         \
	math/ncm_cfg.h                       \
	math/ncm_mpi_job.h                   \
	math/ncm_mpi_job_test.h              \
	math/ncm_mpi_job_fit.h               \
	math/ncm_mpi_job_mcmc.h              \
	math/ncm_util.h                      \
	math/ncm_diff.h                      \
	math/ncm_ode.h                       \
	math/ncm_ode_eval.h                  \
	math/ncm_timer.h                     \
	math/ncm_lapack.h                    \
	math/ncm_vector.h                    \
	math/ncm_matrix.h                    \
	math/ncm_serialize.h                 \
	math/ncm_obj_array.h                 \
	math/ncm_integral1d.h                \
	math/ncm_integral1d_ptr.h            \
	math/ncm_rng.h                       \
	math/ncm_stats_vec.h                 \
	math/ncm_stats_dist1d.h              \
	math/ncm_stats_dist1d_spline.h       \
	math/ncm_stats_dist1d_epdf.h         \
	math/ncm_stats_dist2d.h              \
        math/ncm_stats_dist2d_spline.h       \
	math/ncm_stats_dist_nd.h             \
	math/ncm_stats_dist_nd_kde_gauss.h   \
	math/ncm_bootstrap.h                 \
	math/ncm_spline.h                    \
	math/ncm_spline_gsl.h                \
	math/ncm_spline_cubic.h              \
	math/ncm_spline_cubic_notaknot.h     \
	math/ncm_spline_rbf.h                \
	math/ncm_spline_func.h               \
	math/ncm_spline2d.h                  \
	math/ncm_spline2d_spline.h           \
	math/ncm_spline2d_gsl.h              \
	math/ncm_spline2d_bicubic.h          \
	math/ncm_powspec.h                   \
	math/ncm_powspec_filter.h            \
	math/ncm_powspec_corr3d.h            \
	math/ncm_hoaa.h                      \
	math/ncm_csq1d.h                     \
	math/ncm_func_eval.h                 \
	math/ncm_mpsf_trig_int.h             \
	math/ncm_mpsf_sbessel.h              \
	math/ncm_sf_sbessel.h                \
	math/ncm_sf_sbessel_int.h            \
	math/ncm_sf_spherical_harmonics.h    \
	math/ncm_mpsf_sbessel_int.h          \
	math/ncm_mpsf_0F1.h                  \
	math/ncm_fftlog.h                    \
	math/ncm_fftlog_sbessel_j.h          \
	math/ncm_fftlog_tophatwin2.h         \
	math/ncm_fftlog_gausswin2.h          \
	math/ncm_sparam.h                    \
	math/ncm_vparam.h                    \
	math/ncm_model.h                     \
	math/ncm_model_ctrl.h                \
	math/ncm_model_builder.h             \
	math/ncm_model_mvnd.h                \
	math/ncm_mset.h                      \
	math/ncm_mset_func.h                 \
	math/ncm_mset_func1.h                \
	math/ncm_mset_func_list.h            \
	math/ncm_calc.h                      \
	math/ncm_ode_spline.h                \
	math/ncm_reparam.h                   \
	math/ncm_reparam_linear.h            \
	math/ncm_data.h                      \
	math/ncm_data_dist1d.h               \
	math/ncm_data_dist2d.h               \
	math/ncm_data_gauss.h                \
	math/ncm_data_gauss_cov.h            \
	math/ncm_data_gauss_diag.h           \
	math/ncm_data_poisson.h              \
	math/ncm_data_gauss_cov_mvnd.h       \
	math/ncm_dataset.h                   \
	math/ncm_likelihood.h                \
	math/ncm_mset_trans_kern.h           \
	math/ncm_mset_trans_kern_flat.h      \
	math/ncm_mset_trans_kern_gauss.h     \
	math/ncm_mset_trans_kern_cat.h       \
	math/ncm_prior.h                     \
	math/ncm_prior_gauss.h               \
	math/ncm_prior_gauss_param.h         \
	math/ncm_prior_gauss_func.h          \
	math/ncm_prior_flat.h                \
	math/ncm_prior_flat_param.h          \
	math/ncm_prior_flat_func.h           \
	math/ncm_fit_state.h                 \
	math/ncm_fit.h                       \
	math/ncm_fit_levmar.h                \
	math/ncm_fit_gsl_ls.h                \
	math/ncm_fit_gsl_mm.h                \
	math/ncm_fit_gsl_mms.h               \
	math/ncm_mset_catalog.h              \
	math/ncm_fit_mc.h                    \
	math/ncm_fit_mcbs.h                  \
	math/ncm_fit_mcmc.h                  \
	math/ncm_fit_esmcmc.h                \
	math/ncm_fit_esmcmc_walker.h         \
	math/ncm_fit_esmcmc_walker_stretch.h \
	math/ncm_fit_esmcmc_walker_walk.h    \
	math/ncm_fit_esmcmc_walker_aps.h     \
	math/ncm_lh_ratio1d.h                \
	math/ncm_lh_ratio2d.h                \
	math/ncm_abc.h                       \
	math/ncm_quaternion.h                \
	math/ncm_sphere_map.h                \
	math/ncm_function_cache.h            \
	math/integral.h           \
	math/magnus_iserles_ode.h \
	math/matrix_exp.h         \
	math/poly.h               \
	math/quadrature.h         \
	math/binsplit.h           \
	math/mpq_tree.h           \
	math/grid_one.h           \
	math/dividedifference.h   \
	math/ncm_memory_pool.h

ncm_headers_nogi = \
	math/ncm_gsl_blas_types.h

ncm_headers_internal = \
	gslextras/cqp/gsl_cqp.h           \
	gslextras/cqp/initial_point.h     \
	math/rcm.h                        \
	math/ncm_flapack.h                \
	math/Faddeeva.h                   \
	perturbations/nc_hipert_private.h

nc_headers = \
	nc_hicosmo.h                  \
	nc_hiprim.h                   \
	nc_distance.h                 \
	nc_scalefactor.h              \
	nc_hicosmo_priors.h           \
	nc_recomb.h                   \
	nc_recomb_seager.h            \
	nc_recomb_cbe.h               \
	nc_hireion.h                  \
	nc_hireion_camb.h             \
	nc_hireion_camb_reparam_tau.h \
	nc_powspec_ml.h               \
	nc_powspec_ml_fix_spline.h    \
	nc_powspec_ml_transfer.h      \
	nc_powspec_ml_cbe.h           \
	nc_powspec_mnl.h              \
	nc_powspec_mnl_halofit.h      \
	nc_snia_dist_cov.h            \
	nc_planck_fi.h                \
	nc_planck_fi_cor_tt.h         \
	nc_planck_fi_cor_ttteee.h     \
	nc_cbe_precision.h            \
	nc_cbe.h                      \
	nc_hiqg_1d.h                  \
	model/nc_hicosmo_idem2.h      \
	model/nc_hicosmo_gcg.h        \
	model/nc_hicosmo_de.h         \
	model/nc_hicosmo_de_reparam_ok.h \
	model/nc_hicosmo_de_reparam_cmb.h \
	model/nc_hicosmo_de_cpl.h  \
	model/nc_hicosmo_de_jbp.h     \
	model/nc_hicosmo_de_xcdm.h    \
	model/nc_hicosmo_lcdm.h       \
	model/nc_hicosmo_qconst.h     \
	model/nc_hicosmo_qlinear.h    \
	model/nc_hicosmo_qspline.h    \
	model/nc_hicosmo_qrbf.h       \
	model/nc_hicosmo_qgrw.h       \
	model/nc_hicosmo_Vexp.h       \
	model/nc_hiprim_power_law.h   \
	model/nc_hiprim_atan.h        \
	model/nc_hiprim_expc.h        \
	model/nc_hiprim_bpl.h         \
	model/nc_hiprim_sbpl.h        \
	lss/nc_window.h               \
	lss/nc_window_tophat.h        \
	lss/nc_window_gaussian.h      \
	lss/nc_transfer_func.h        \
	lss/nc_transfer_func_bbks.h   \
	lss/nc_transfer_func_eh.h     \
	lss/nc_transfer_func_camb.h   \
	lss/nc_growth_func.h          \
	lss/nc_density_profile.h      \
	lss/nc_density_profile_nfw.h  \
	lss/nc_density_profile_einasto.h \
	lss/nc_density_profile_dk14.h \
	lss/nc_galaxy_acf.h           \
	lss/nc_multiplicity_func.h    \
	lss/nc_multiplicity_func_ps.h           \
	lss/nc_multiplicity_func_st.h           \
	lss/nc_multiplicity_func_jenkins.h      \
	lss/nc_multiplicity_func_warren.h       \
	lss/nc_multiplicity_func_tinker.h       \
	lss/nc_multiplicity_func_tinker_mean.h  \
	lss/nc_multiplicity_func_tinker_crit.h  \
	lss/nc_multiplicity_func_tinker_mean_normalized.h \
	lss/nc_multiplicity_func_crocce.h      \
	lss/nc_halo_mass_function.h            \
        lss/nc_halo_bias_type.h                \
        lss/nc_halo_bias_type_ps.h             \
        lss/nc_halo_bias_type_st_spher.h       \
        lss/nc_halo_bias_type_st_ellip.h       \
        lss/nc_halo_bias_type_tinker.h         \
        lss/nc_halo_bias_func.h                \
        lss/nc_cluster_redshift.h              \
        lss/nc_cluster_redshift_nodist.h       \
        lss/nc_cluster_photoz_gauss.h          \
        lss/nc_cluster_photoz_gauss_global.h   \
        lss/nc_cluster_mass.h                  \
        lss/nc_cluster_mass_nodist.h           \
        lss/nc_cluster_mass_lnnormal.h         \
        lss/nc_cluster_mass_vanderlinde.h      \
        lss/nc_cluster_mass_benson.h           \
        lss/nc_cluster_mass_benson_xray.h      \
        lss/nc_cluster_mass_plcl.h             \
        lss/nc_cluster_mass_ascaso.h           \
        lss/nc_cluster_abundance.h             \
        lss/nc_cluster_pseudo_counts.h         \
        lss/nc_galaxy_redshift.h               \
        lss/nc_galaxy_redshift_spec.h          \
        lss/nc_galaxy_redshift_spline.h        \
        lss/nc_cor_cluster_cmb_lens_limber.h   \
        lss/nc_wl_surface_mass_density.h       \
        lss/nc_reduced_shear_cluster_mass.h    \
        lss/nc_reduced_shear_calib.h           \
        lss/nc_reduced_shear_calib_wtg.h       \
        lss/nc_galaxy_selfunc.h                \
	data/nc_data_snia.h                    \
	data/nc_data_dist_mu.h                 \
	data/nc_data_snia_cov.h                \
	data/nc_data_hubble.h                  \
	data/nc_data_hubble_bao.h              \
	data/nc_data_bao_rdv.h                 \
	data/nc_data_bao_dvdv.h                \
	data/nc_data_bao_a.h                   \
	data/nc_data_bao_dv.h                  \
	data/nc_data_bao_empirical_fit.h       \
	data/nc_data_bao_empirical_fit_2d.h    \
	data/nc_data_bao_dhr_dar.h             \
	data/nc_data_bao_dmr_hr.h              \
	data/nc_data_bao.h                     \
	data/nc_data_cmb_shift_param.h         \
	data/nc_data_cmb_dist_priors.h         \
	data/nc_data_cmb.h                     \
	data/nc_data_cluster_ncount.h          \
	data/nc_data_cluster_poisson.h         \
	data/nc_data_cluster_counts_box_poisson.h \
	data/nc_data_cluster_pseudo_counts.h      \
	data/nc_data_reduced_shear_cluster_mass.h \
	data/nc_data_planck_lkl.h                 \
	data/nc_data_xcor.h                       \
	abc/nc_abc_cluster_ncount.h               \
	perturbations/nc_hipert.h                 \
	perturbations/nc_hipert_bg_var.h          \
	perturbations/nc_hipert_wkb.h             \
	perturbations/nc_hipert_itwo_fluids.h     \
	perturbations/nc_hipert_adiab.h           \
	perturbations/nc_hipert_gw.h              \
	perturbations/nc_hipert_two_fluids.h      \
	perturbations/nc_hipert_boltzmann.h       \
	perturbations/nc_hipert_boltzmann_std.h   \
	perturbations/nc_hipert_boltzmann_cbe.h   \
	perturbations/nc_hipert_first_order.h     \
	perturbations/nc_hipert_grav.h            \
	perturbations/nc_hipert_grav_einstein.h   \
	perturbations/nc_hipert_comp.h            \
	perturbations/nc_hipert_comp_pb.h         \
	xcor/nc_xcor.h                            \
	xcor/nc_xcor_AB.h                         \
	xcor/nc_xcor_limber_kernel.h              \
	xcor/nc_xcor_limber_kernel_gal.h          \
	xcor/nc_xcor_limber_kernel_CMB_lensing.h  \
	xcor/nc_xcor_limber_kernel_weak_lensing.h

nc_sources = \
	nc_hicosmo.c                  \
	nc_hiprim.c                   \
	nc_distance.c                 \
	nc_scalefactor.c              \
	nc_hicosmo_priors.c           \
	nc_recomb.c                   \
	nc_recomb_seager.c            \
	nc_recomb_cbe.c               \
	nc_hireion.c                  \
	nc_hireion_camb.c             \
	nc_hireion_camb_reparam_tau.c \
	nc_powspec_ml.c               \
	nc_powspec_ml_fix_spline.c    \
	nc_powspec_ml_transfer.c      \
	nc_powspec_ml_cbe.c           \
	nc_powspec_mnl.c              \
	nc_powspec_mnl_halofit.c      \
	nc_snia_dist_cov.c            \
	nc_planck_fi.c                \
	nc_planck_fi_cor_tt.c         \
	nc_planck_fi_cor_ttteee.c     \
	nc_cbe_precision.c            \
	nc_cbe.c                      \
	nc_hiqg_1d.c                  \
	model/nc_hicosmo_idem2.c      \
	model/nc_hicosmo_gcg.c        \
	model/nc_hicosmo_de.c         \
	model/nc_hicosmo_de_reparam_ok.c \
	model/nc_hicosmo_de_reparam_cmb.c \
	model/nc_hicosmo_de_cpl.c  \
	model/nc_hicosmo_de_jbp.c     \
	model/nc_hicosmo_de_xcdm.c    \
	model/nc_hicosmo_lcdm.c       \
	model/nc_hicosmo_qconst.c     \
	model/nc_hicosmo_qlinear.c    \
	model/nc_hicosmo_qspline.c    \
	model/nc_hicosmo_qrbf.c       \
	model/nc_hicosmo_qgrw.c       \
	model/nc_hicosmo_Vexp.c       \
	model/nc_hiprim_power_law.c   \
	model/nc_hiprim_atan.c        \
	model/nc_hiprim_expc.c        \
	model/nc_hiprim_bpl.c         \
	model/nc_hiprim_sbpl.c        \
	lss/nc_window.c               \
	lss/nc_window_tophat.c        \
	lss/nc_window_gaussian.c      \
	lss/nc_transfer_func.c        \
	lss/nc_transfer_func_bbks.c   \
	lss/nc_transfer_func_eh.c     \
	lss/nc_transfer_func_camb.c   \
	lss/nc_growth_func.c          \
	lss/nc_density_profile.c      \
	lss/nc_density_profile_nfw.c  \
	lss/nc_density_profile_einasto.c \
	lss/nc_density_profile_dk14.c \
	lss/nc_galaxy_acf.c           \
	lss/nc_multiplicity_func.c    \
	lss/nc_multiplicity_func_ps.c           \
	lss/nc_multiplicity_func_st.c           \
	lss/nc_multiplicity_func_jenkins.c      \
	lss/nc_multiplicity_func_warren.c       \
	lss/nc_multiplicity_func_tinker.c       \
	lss/nc_multiplicity_func_tinker_mean.c  \
	lss/nc_multiplicity_func_tinker_crit.c  \
	lss/nc_multiplicity_func_tinker_mean_normalized.c \
	lss/nc_multiplicity_func_crocce.c      \
	lss/nc_halo_mass_function.c            \
	lss/nc_halo_bias_type.c                \
	lss/nc_halo_bias_type_ps.c             \
	lss/nc_halo_bias_type_st_spher.c       \
	lss/nc_halo_bias_type_st_ellip.c       \
	lss/nc_halo_bias_type_tinker.c         \
	lss/nc_halo_bias_func.c                \
	lss/nc_cluster_redshift.c              \
	lss/nc_cluster_redshift_nodist.c       \
	lss/nc_cluster_photoz_gauss.c          \
	lss/nc_cluster_photoz_gauss_global.c   \
	lss/nc_cluster_mass.c                  \
	lss/nc_cluster_mass_nodist.c           \
	lss/nc_cluster_mass_lnnormal.c         \
	lss/nc_cluster_mass_vanderlinde.c      \
	lss/nc_cluster_mass_benson.c           \
	lss/nc_cluster_mass_benson_xray.c      \
	lss/nc_cluster_mass_plcl.c             \
        lss/nc_cluster_mass_ascaso.c           \
	lss/nc_cluster_abundance.c             \
	lss/nc_cluster_pseudo_counts.c         \
	lss/nc_galaxy_redshift.c               \
	lss/nc_galaxy_redshift_spec.c          \
	lss/nc_galaxy_redshift_spline.c        \
	lss/nc_cor_cluster_cmb_lens_limber.c   \
	lss/nc_wl_surface_mass_density.c       \
	lss/nc_reduced_shear_cluster_mass.c    \
	lss/nc_reduced_shear_calib.c           \
	lss/nc_reduced_shear_calib_wtg.c       \
        lss/nc_galaxy_selfunc.c                \
	data/nc_data_snia.c                    \
	data/nc_data_dist_mu.c                 \
	data/nc_data_snia_cov.c                \
	data/nc_data_hubble.c                  \
	data/nc_data_hubble_bao.c              \
	data/nc_data_bao_rdv.c                 \
	data/nc_data_bao_dvdv.c                \
	data/nc_data_bao_a.c                   \
	data/nc_data_bao_dv.c                  \
	data/nc_data_bao_empirical_fit.c       \
	data/nc_data_bao_empirical_fit_2d.c    \
	data/nc_data_bao_dhr_dar.c             \
	data/nc_data_bao_dmr_hr.c              \
	data/nc_data_bao.c                     \
	data/nc_data_cmb_shift_param.c         \
	data/nc_data_cmb_dist_priors.c         \
	data/nc_data_cmb.c                     \
	data/nc_data_cluster_ncount.c          \
	data/nc_data_cluster_poisson.c         \
	data/nc_data_cluster_counts_box_poisson.c \
	data/nc_data_cluster_pseudo_counts.c      \
	data/nc_data_reduced_shear_cluster_mass.c \
	data/nc_data_planck_lkl.c                 \
	data/nc_data_xcor.c                       \
	abc/nc_abc_cluster_ncount.c               \
	perturbations/nc_hipert.c                 \
	perturbations/nc_hipert_bg_var.c          \
	perturbations/nc_hipert_wkb.c             \
	perturbations/nc_hipert_itwo_fluids.c     \
	perturbations/nc_hipert_adiab.c           \
	perturbations/nc_hipert_gw.c              \
	perturbations/nc_hipert_two_fluids.c      \
	perturbations/nc_hipert_boltzmann.c       \
	perturbations/nc_hipert_boltzmann_std.c   \
	perturbations/nc_hipert_boltzmann_cbe.c   \
	perturbations/nc_hipert_first_order.c     \
	perturbations/nc_hipert_grav.c            \
	perturbations/nc_hipert_grav_einstein.c   \
	perturbations/nc_hipert_comp.c            \
	perturbations/nc_hipert_comp_pb.c         \
	xcor/nc_xcor.c                            \
	xcor/nc_xcor_AB.c                         \
	xcor/nc_xcor_limber_kernel.c              \
	xcor/nc_xcor_limber_kernel_gal.c          \
	xcor/nc_xcor_limber_kernel_CMB_lensing.c  \
	xcor/nc_xcor_limber_kernel_weak_lensing.c

libnumcosmo_la_includedir = \
	$(pkgincludedir)

ncm_enum_types_built_headers =

ncm_enum_types_headers = \
	ncm_enum_types.h

nc_enum_types_headers = \
	nc_enum_types.h

ncm_enum_types_sources = \
	ncm_enum_types.c

nc_enum_types_sources = \
	nc_enum_types.c

ncm_cond_sources =

ncm_cond_headers =

if HAVE_LIBNLOPT

ncm_fit_nlopt_enum_meta.h: ncm_fit_nlopt_enum_meta.c
	$(AM_V_GEN)$(CPP) $(AM_CFLAGS) $(srcdir)/ncm_fit_nlopt_enum_meta.c -o ncm_fit_nlopt_enum_meta.h.tmp
	$(AM_V_GEN)$(GREP) -v '^#' ncm_fit_nlopt_enum_meta.h.tmp > ncm_fit_nlopt_enum_meta.h.tmp1
	$(AM_V_GEN)$(SED) 's/nlopt_algorithm/NcmFitNloptAlgorithm/' ncm_fit_nlopt_enum_meta.h.tmp1 > ncm_fit_nlopt_enum_meta.h.tmp2
	$(AM_V_GEN)$(SED) 's/nlopt_result/NcmFitNloptResult/' ncm_fit_nlopt_enum_meta.h.tmp2 > ncm_fit_nlopt_enum_meta.h
	$(AM_V_GEN)rm ncm_fit_nlopt_enum_meta.h.tmp
	$(AM_V_GEN)rm ncm_fit_nlopt_enum_meta.h.tmp1
	$(AM_V_GEN)rm ncm_fit_nlopt_enum_meta.h.tmp2

ncm_fit_nlopt_enum.h: ncm_fit_nlopt_enum_meta.h $(srcdir)/ncm_fit_nlopt_enum.h.template
	$(AM_V_GEN) ( $(GLIB_MKENUMS) --template $(srcdir)/ncm_fit_nlopt_enum.h.template ncm_fit_nlopt_enum_meta.h ) > \
	  ncm_fit_nlopt_enum.h.tmp && mv ncm_fit_nlopt_enum.h.tmp ncm_fit_nlopt_enum.h

ncm_fit_nlopt_enum.c: ncm_fit_nlopt_enum_meta.h $(srcdir)/ncm_fit_nlopt_enum.c.template
	$(AM_V_GEN) ( $(GLIB_MKENUMS) --template $(srcdir)/ncm_fit_nlopt_enum.c.template ncm_fit_nlopt_enum_meta.h ) > \
	  ncm_fit_nlopt_enum.c.tmp && mv ncm_fit_nlopt_enum.c.tmp ncm_fit_nlopt_enum.c

ncm_enum_types_built_headers += \
	ncm_fit_nlopt_enum_meta.h

ncm_enum_types_sources += \
	ncm_fit_nlopt_enum.c

ncm_cond_sources += \
	math/ncm_fit_nlopt.c

ncm_enum_types_headers += \
	ncm_fit_nlopt_enum.h

ncm_cond_headers += \
	math/ncm_fit_nlopt.h

endif

BUILT_SOURCES = \
	$(ncm_enum_types_headers) \
	$(ncm_enum_types_sources) \
	$(nc_enum_types_headers)  \
	$(nc_enum_types_sources)  \
	$(ncm_enum_types_built_headers)

ncm_enum_types.h: $(ncm_headers) ncm_enum_types.h.template
	$(AM_V_GEN) ( cd $(srcdir) && $(GLIB_MKENUMS) --template ncm_enum_types.h.template $(ncm_headers) $(ncm_cond_headers) ) > \
	  ncm_enum_types.h.tmp && mv ncm_enum_types.h.tmp ncm_enum_types.h

ncm_enum_types.c: $(ncm_headers) ncm_enum_types.c.template
	$(AM_V_GEN) ( cd $(srcdir) && $(GLIB_MKENUMS) --template ncm_enum_types.c.template $(ncm_headers) $(ncm_cond_headers) ) > \
	  ncm_enum_types.c.tmp && mv ncm_enum_types.c.tmp ncm_enum_types.c

nc_enum_types.h: $(nc_headers) nc_enum_types.h.template
	$(AM_V_GEN) ( cd $(srcdir) && $(GLIB_MKENUMS) --template nc_enum_types.h.template $(nc_headers) ) > \
	  nc_enum_types.h.tmp && mv nc_enum_types.h.tmp nc_enum_types.h

nc_enum_types.c: $(nc_headers) nc_enum_types.c.template
	$(AM_V_GEN) ( cd $(srcdir) && $(GLIB_MKENUMS) --template nc_enum_types.c.template $(nc_headers) ) > \
	  nc_enum_types.c.tmp && mv nc_enum_types.c.tmp nc_enum_types.c

libnumcosmo_la_SOURCES = \
	$(ncm_sources)          \
	$(ncm_headers)          \
	$(ncm_sources_nogi)     \
	$(ncm_headers_nogi)     \
	$(ncm_headers_internal) \
	$(nc_sources)           \
	$(nc_headers)           \
	$(ncm_cond_sources)     \
	$(ncm_cond_headers)     \
	config_extra.h          \
	numcosmo-math.h         \
	numcosmo.h

nodist_libnumcosmo_la_SOURCES = \
	build_cfg.h               \
	$(ncm_enum_types_headers) \
	$(ncm_enum_types_sources) \
	$(nc_enum_types_headers)  \
	$(nc_enum_types_sources)

nobase_libnumcosmo_la_include_HEADERS = \
	$(ncm_headers)            \
	$(ncm_headers_nogi)       \
	$(nc_headers)             \
	$(ncm_cond_headers)       \
	math/binsplit_eval.c      \
	numcosmo-math.h           \
	numcosmo.h

nodist_libnumcosmo_la_include_HEADERS = \
	build_cfg.h               \
	$(ncm_enum_types_headers) \
	$(nc_enum_types_headers)

EXTRA_DIST = \
	$(nobase_libnumcosmo_la_include_HEADERS) \
	build_cfg.h.in                  \
	ncm_enum_types.h.template       \
	ncm_enum_types.c.template       \
	ncm_fit_nlopt_enum.h.template   \
	ncm_fit_nlopt_enum.c.template   \
	nc_enum_types.h.template        \
	nc_enum_types.c.template        \
	ncm_gir_scan.h                  \
	ncm_fit_nlopt_enum_meta.c       \
	math/ncm_sphere_map_block.c     \
	math/rquantile.c                \
	math/gsl_rstat.h

CLEANFILES = $(BUILT_SOURCES)

-include $(INTROSPECTION_MAKEFILE)
INTROSPECTION_SCANNER_ENV = CC="$(CC)"
INTROSPECTION_GIRS =
INTROSPECTION_SCANNER_ARGS = --add-include-path=$(srcdir) -I$(srcdir)/../ -I$(top_builddir) --warn-all

INTROSPECTION_COMPILER_ARGS = --includedir=$(srcdir)

if HAVE_INTROSPECTION

ncm_introspection_sources = \
	ncm_gir_scan.h            \
	$(ncm_sources)            \
	$(ncm_headers)            \
	$(ncm_cond_sources)       \
	$(ncm_cond_headers)       \
	$(ncm_enum_types_headers) \
	$(ncm_enum_types_sources) \
	numcosmo-math.h

nc_introspection_sources = \
	ncm_gir_scan.h           \
	$(nc_sources)            \
	$(nc_headers)            \
	$(nc_enum_types_headers) \
	$(nc_enum_types_sources) \
	numcosmo.h

NumCosmoMath-1.0.gir: libnumcosmo.la $(INTROSPECTION_SCANNER)
NumCosmoMath_1_0_gir_VERSION      = 1.0
NumCosmoMath_1_0_gir_INCLUDES     = GObject-2.0
NumCosmoMath_1_0_gir_CFLAGS       = $(MY_CFLAGS) -DNUMCOSMO_GIR_SCAN
NumCosmoMath_1_0_gir_LIBS         = libnumcosmo.la
NumCosmoMath_1_0_gir_FILES        = $(ncm_introspection_sources)
NumCosmoMath_1_0_gir_SCANNERFLAGS = \
	--c-include="numcosmo/numcosmo-math.h"  \
	-n NumCosmoMath --identifier-prefix=Ncm

INTROSPECTION_GIRS += NumCosmoMath-1.0.gir

NumCosmo-1.0.gir: libnumcosmo.la $(INTROSPECTION_SCANNER) NumCosmoMath-1.0.gir
NumCosmo_1_0_gir_VERSION      = 1.0
NumCosmo_1_0_gir_INCLUDES     = GObject-2.0
NumCosmo_1_0_gir_CFLAGS       = $(MY_CFLAGS) -DNUMCOSMO_GIR_SCAN
NumCosmo_1_0_gir_LIBS         = libnumcosmo.la
NumCosmo_1_0_gir_FILES        = $(nc_introspection_sources)
NumCosmo_1_0_gir_SCANNERFLAGS = \
	--include-uninstalled $(builddir)/NumCosmoMath-1.0.gir \
	--c-include="numcosmo/numcosmo.h" \
	-n NumCosmo --identifier-prefix=Nc

INTROSPECTION_GIRS += NumCosmo-1.0.gir

girdir = $(datadir)/gir-1.0
gir_DATA = $(INTROSPECTION_GIRS)

typelibdir = $(libdir)/girepository-1.0
typelib_DATA = $(INTROSPECTION_GIRS:.gir=.typelib)

CLEANFILES += $(gir_DATA) $(typelib_DATA)
EXTRA_DIST += $(gir_DATA) $(typelib_DATA)

if HAVE_VAPIGEN

vapidir = $(datadir)/vala/vapi
vapi_DATA = numcosmo-1.0.vapi numcosmomath-1.0.vapi

numcosmomath_vapiopt = \
	--vapidir=$(builddir)      \
	--girdir=$(builddir)       \
	--pkg GObject-2.0          \
	--library numcosmomath-1.0

numcosmomath-1.0.vapi: NumCosmoMath-1.0.gir
	$(AM_V_GEN)$(VAPIGEN) $(numcosmomath_vapiopt) $<
	@touch $@

numcosmo_vapiopt = \
	--vapidir=$(builddir)  \
	--girdir=$(builddir)   \
	--pkg numcosmomath-1.0 \
	--library numcosmo-1.0

numcosmo-1.0.vapi: NumCosmo-1.0.gir numcosmomath-1.0.vapi
	$(AM_V_GEN)$(VAPIGEN) $(numcosmo_vapiopt) $<
	@touch $@

CLEANFILES += $(vapi_DATA)
EXTRA_DIST += $(vapi_DATA)

endif # HAVE_INTROSPECTION

endif # HAVE_VAPIGEN
