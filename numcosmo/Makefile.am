## Process this file with automake to produce Makefile.in
include $(top_srcdir)/Makefile.decl

AM_CFLAGS = \
	 -Wall \
	 -DG_LOG_DOMAIN=\"NUMCOSMO\" \
	 $(GLIB_CFLAGS) \
	 $(SUNDIALS_CFLAGS) \
	 $(GSL_CFLAGS) \
	 $(FFTW3_CFLAGS) \
	 $(CFITSIO_CFLAGS) \
	 $(NLOPT_CFLAGS) \
	 -I$(top_srcdir)

SUBDIRS = 

lib_LTLIBRARIES = libnumcosmo.la

libnumcosmo_la_LIBADD = \
	$(SUNDIALS_LIBS) \
	$(FFTW3_LIBS) \
	$(GLIB_LIBS) \
	$(NLOPT_LIBS) \
	$(GSL_LIBS) \
	$(SQLITE3_LIBS) \
	$(CFITSIO_LIBS)

ncm_sources = \
	ncm_build_inline.c                 \
	math/ncm_c.c                       \
	math/ncm_cfg.c                     \
	math/ncm_lapack.c                  \
	math/ncm_vector.c                  \
	math/ncm_matrix.c                  \
	math/ncm_spline.c                  \
	math/ncm_spline_gsl.c              \
	math/ncm_spline_cubic.c            \
	math/ncm_spline_cubic_notaknot.c   \
	math/ncm_spline_func.c             \
	math/ncm_spline2d.c                \
	math/ncm_spline2d_spline.c         \
	math/ncm_spline2d_gsl.c            \
	math/ncm_spline2d_bicubic.c        \
	math/ncm_func_eval.c               \
	math/ncm_mpsf_trig_int.c           \
	math/ncm_mpsf_sbessel.c            \
	math/ncm_sf_sbessel.c              \
	math/ncm_sf_sbessel_int.c          \
	math/ncm_mpsf_sbessel_int.c        \
	math/ncm_mpsf_0F1.c                \
	math/ncm_fftlog.c                  \
	math/ncm_sparam.c                  \
	math/ncm_vparam.c                  \
	math/ncm_model.c                   \
	math/ncm_model_ctrl.c              \
	math/ncm_mset.c                    \
	math/ncm_mset_func.c               \
	math/ncm_ode_spline.c              \
	math/ncm_reparam.c                 \
	math/ncm_reparam_linear.c          \
	math/ncm_data.c                    \
	math/ncm_data_dist1d.c             \
	math/ncm_data_gauss.c              \
	math/ncm_data_gauss_cov.c          \
	math/ncm_data_gauss_diag.c         \
	math/ncm_data_poisson.c            \
	math/ncm_dataset.c                 \
	math/ncm_likelihood.c              \
	math/ncm_priors.c                  \
	math/ncm_fit_state.c               \
	math/ncm_fit.c                     \
	math/ncm_fit_gsl_ls.c              \
	math/ncm_fit_gsl_mm.c              \
	math/ncm_fit_gsl_mms.c             \
	math/ncm_fit_mc.c                  \
	math/ncm_lh_ratio1d.c              \
	math/ncm_lh_ratio2d.c              \
	math/poly.c               \
	math/integral.c           \
	math/function_cache.c     \
	math/quaternion.c         \
	math/cvode_util.c         \
	math/quadrature.c         \
	math/matrix_exp.c         \
	math/magnus_iserles_ode.c \
	math/binsplit.c           \
	math/mpq_tree.c           \
	math/grid_one.c           \
	math/dividedifference.c   \
	math/memory_pool.c        \
	math/util.c               \
	sphere/map.c              \
	sphere/healpix.c

ncm_headers = \
	math/ncm_c.h                     \
	math/ncm_cfg.h                   \
	math/ncm_lapack.h                \
	math/ncm_vector.h                \
	math/ncm_matrix.h                \
	math/ncm_spline.h                \
	math/ncm_spline_gsl.h            \
	math/ncm_spline_cubic.h          \
	math/ncm_spline_cubic_notaknot.h \
	math/ncm_spline_func.h           \
	math/ncm_spline2d.h              \
	math/ncm_spline2d_spline.h       \
	math/ncm_spline2d_gsl.h          \
	math/ncm_spline2d_bicubic.h      \
	math/ncm_func_eval.h             \
	math/ncm_mpsf_trig_int.h         \
	math/ncm_mpsf_sbessel.h          \
	math/ncm_sf_sbessel.h            \
	math/ncm_sf_sbessel_int.h        \
	math/ncm_mpsf_sbessel_int.h      \
	math/ncm_mpsf_0F1.h              \
	math/ncm_fftlog.h                \
	math/ncm_sparam.h                \
	math/ncm_vparam.h                \
	math/ncm_model.h                 \
	math/ncm_model_ctrl.h            \
	math/ncm_mset.h                  \
	math/ncm_mset_func.h             \
	math/ncm_ode_spline.h            \
	math/ncm_reparam.h               \
	math/ncm_reparam_linear.h        \
	math/ncm_data.h                  \
	math/ncm_data_dist1d.h           \
	math/ncm_data_gauss.h            \
	math/ncm_data_gauss_cov.h        \
	math/ncm_data_gauss_diag.h       \
	math/ncm_data_poisson.h          \
	math/ncm_dataset.h               \
	math/ncm_likelihood.h            \
	math/ncm_priors.h                \
	math/ncm_fit_state.h             \
	math/ncm_fit.h                   \
	math/ncm_fit_gsl_ls.h            \
	math/ncm_fit_gsl_mm.h            \
	math/ncm_fit_gsl_mms.h           \
	math/ncm_fit_mc.h                \
	math/ncm_lh_ratio1d.h            \
	math/ncm_lh_ratio2d.h            \
	math/cvode_util.h         \
	math/function_cache.h     \
	math/integral.h           \
	math/magnus_iserles_ode.h \
	math/matrix_exp.h         \
	math/poly.h               \
	math/quadrature.h         \
	math/quaternion.h         \
	math/binsplit.h           \
	math/mpq_tree.h           \
	math/grid_one.h           \
	math/dividedifference.h   \
	math/util.h               \
	math/memory_pool.h        \
	sphere/map.h              \
	sphere/healpix.h

ncm_headers_internal = 

nc_headers = \
	nc_hicosmo.h                  \
	nc_distance.h                 \
	nc_hicosmo_priors.h           \
	nc_recomb.h                   \
	nc_recomb_seager.h            \
	nc_snia_dist_cov.h            \
	model/nc_hicosmo_de.h         \
	model/nc_hicosmo_de_linder.h  \
	model/nc_hicosmo_de_pad.h     \
	model/nc_hicosmo_de_qe.h      \
	model/nc_hicosmo_de_xcdm.h    \
	model/nc_hicosmo_lcdm.h       \
	model/nc_hicosmo_qconst.h     \
	model/nc_hicosmo_qlinear.h    \
	model/nc_hicosmo_qpw.h        \
	model/nc_hicosmo_qspline.h    \
	lss/nc_window.h               \
	lss/nc_window_tophat.h        \
	lss/nc_window_gaussian.h      \
	lss/nc_transfer_func.h        \
	lss/nc_transfer_func_bbks.h   \
	lss/nc_transfer_func_eh.h     \
	lss/nc_transfer_func_camb.h   \
	lss/nc_transfer_func_pert.h   \
	lss/nc_growth_func.h          \
	lss/nc_matter_var.h           \
	lss/nc_galaxy_acf.h           \
	lss/nc_multiplicity_func.h    \
	lss/nc_multiplicity_func_ps.h \
	lss/nc_multiplicity_func_st.h \
	lss/nc_multiplicity_func_jenkins.h \
	lss/nc_multiplicity_func_warren.h \
	lss/nc_multiplicity_func_tinker.h \
	lss/nc_multiplicity_func_tinker_mean.h \
	lss/nc_multiplicity_func_tinker_crit.h \
	lss/nc_mass_function.h           \
        lss/nc_halo_bias_type.h          \
        lss/nc_halo_bias_type_ps.h       \
        lss/nc_halo_bias_type_st_spher.h \
        lss/nc_halo_bias_type_st_ellip.h \
        lss/nc_halo_bias_type_tinker.h   \
        lss/nc_halo_bias_func.h          \
        lss/nc_cluster_redshift.h            \
        lss/nc_cluster_redshift_nodist.h     \
        lss/nc_cluster_photoz_gauss.h        \
        lss/nc_cluster_photoz_gauss_global.h \
        lss/nc_cluster_mass.h                \
        lss/nc_cluster_mass_nodist.h         \
        lss/nc_cluster_mass_lnnormal.h       \
        lss/nc_cluster_mass_vanderlinde.h    \
        lss/nc_cluster_mass_benson.h         \
        lss/nc_cluster_mass_benson_xray.h    \
        lss/nc_cluster_abundance.h      \
	nc_data_snia.h                  \
	nc_data_dist_mu.h               \
	nc_data_snia_cov.h              \
	nc_data_hubble.h                \
	nc_data_hubble_bao.h            \
	nc_data_bao_rdv.h               \
	nc_data_bao_dvdv.h              \
	nc_data_bao_a.h                 \
	nc_data_bao_dv.h                \
	nc_data_bao.h                   \
	nc_data_cmb_shift_param.h       \
	nc_data_cmb_dist_priors.h       \
	nc_data_cmb.h                   \
	nc_data_cluster_ncount.h        \
	nc_data_cluster_poisson.h       \
	scalefactor.h                   \
	model/quantum_gravity.h         \
	perturbations/linear.h          \
	perturbations/covariance.h      \
	perturbations/linear_internal.h \
	perturbations/hydrodyn_adiabatic.h

nc_sources = \
	nc_hicosmo.c                 \
	nc_distance.c                \
	nc_hicosmo_priors.c          \
	nc_recomb.c                  \
	nc_recomb_seager.c           \
	nc_snia_dist_cov.c           \
	model/nc_hicosmo_de.c        \
	model/nc_hicosmo_de_linder.c \
	model/nc_hicosmo_de_pad.c    \
	model/nc_hicosmo_de_qe.c     \
	model/nc_hicosmo_de_xcdm.c   \
	model/nc_hicosmo_lcdm.c      \
	model/nc_hicosmo_qconst.c    \
	model/nc_hicosmo_qlinear.c   \
	model/nc_hicosmo_qpw.c       \
	model/nc_hicosmo_qspline.c    \
	lss/nc_window.c              \
	lss/nc_window_tophat.c       \
	lss/nc_window_gaussian.c     \
	lss/nc_transfer_func.c       \
	lss/nc_transfer_func_bbks.c  \
	lss/nc_transfer_func_eh.c    \
	lss/nc_transfer_func_camb.c  \
	lss/nc_transfer_func_pert.c  \
	lss/nc_growth_func.c         \
	lss/nc_matter_var.c          \
	lss/nc_galaxy_acf.c          \
	lss/nc_multiplicity_func.c   \
	lss/nc_multiplicity_func_ps.c          \
	lss/nc_multiplicity_func_st.c          \
	lss/nc_multiplicity_func_jenkins.c     \
	lss/nc_multiplicity_func_warren.c      \
	lss/nc_multiplicity_func_tinker.c      \
	lss/nc_multiplicity_func_tinker_mean.c \
	lss/nc_multiplicity_func_tinker_crit.c \
	lss/nc_mass_function.c                 \
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
	lss/nc_cluster_abundance.c             \
	nc_data_snia.c               \
	nc_data_dist_mu.c            \
	nc_data_snia_cov.c           \
	nc_data_hubble.c             \
	nc_data_hubble_bao.c         \
	nc_data_bao_rdv.c            \
	nc_data_bao_dvdv.c           \
	nc_data_bao_a.c              \
	nc_data_bao_dv.c             \
	nc_data_bao.c                \
	nc_data_cmb_shift_param.c    \
	nc_data_cmb_dist_priors.c    \
	nc_data_cmb.c                \
	nc_data_cluster_ncount.c     \
	nc_data_cluster_poisson.c    \
	scalefactor.c                   \
	model/quantum_gravity.c         \
	perturbations/linear.c          \
	perturbations/covariance.c      \
	perturbations/linear_cvodes.c   \
	perturbations/linear_gsl_ode.c  \
	perturbations/linear_gsl_ode2.c \
	perturbations/hydrodyn_adiabatic.c

libnumcosmo_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)

cosmo_includedir = \
	$(pkgincludedir)

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
	$(CPP) $(AM_CFLAGS) $(srcdir)/ncm_fit_nlopt_enum_meta.c -o ncm_fit_nlopt_enum_meta.h.tmp
	$(GREP) -v '^#' ncm_fit_nlopt_enum_meta.h.tmp > ncm_fit_nlopt_enum_meta.h.tmp1
	$(SED) 's/nlopt_algorithm/NcmFitNloptAlgorithm/' ncm_fit_nlopt_enum_meta.h.tmp1 > ncm_fit_nlopt_enum_meta.h.tmp2
	$(SED) 's/nlopt_result/NcmFitNloptResult/' ncm_fit_nlopt_enum_meta.h.tmp2 > ncm_fit_nlopt_enum_meta.h

ncm_fit_nlopt_enum.h: ncm_fit_nlopt_enum_meta.h
	$(AM_V_GEN) ( $(GLIB_MKENUMS) --template $(srcdir)/ncm_fit_nlopt_enum.h.template ncm_fit_nlopt_enum_meta.h ) > \
	  ncm_fit_nlopt_enum.h.tmp && mv ncm_fit_nlopt_enum.h.tmp ncm_fit_nlopt_enum.h

ncm_fit_nlopt_enum.c: ncm_fit_nlopt_enum_meta.h
	$(AM_V_GEN) ( $(GLIB_MKENUMS) --template $(srcdir)/ncm_fit_nlopt_enum.c.template ncm_fit_nlopt_enum_meta.h ) > \
	  ncm_fit_nlopt_enum.c.tmp && mv ncm_fit_nlopt_enum.c.tmp ncm_fit_nlopt_enum.c

ncm_enum_types_sources += \
	ncm_fit_nlopt_enum.c

ncm_cond_sources += \
	math/ncm_fit_nlopt.c

ncm_enum_types_headers += \
	ncm_fit_nlopt_enum.h

ncm_cond_headers += \
	math/ncm_fit_nlopt.h

endif

if HAVE_LIBLEVMAR

ncm_cond_headers += \
	math/ncm_fit_levmar.h

ncm_cond_sources += \
	math/ncm_fit_levmar.c

endif

BUILT_SOURCES = \
	$(ncm_enum_types_headers) \
	$(ncm_enum_types_sources) \
	$(nc_enum_types_headers)  \
	$(nc_enum_types_sources)

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
	$(ncm_sources)      \
	$(ncm_headers)      \
	$(nc_sources)       \
	$(nc_headers)       \
	$(ncm_cond_sources) \
	$(ncm_cond_headers) \
	numcosmo-math.h     \
	numcosmo.h

nodist_libnumcosmo_la_SOURCES = \
	build_cfg.h               \
	$(ncm_enum_types_headers) \
	$(ncm_enum_types_sources) \
	$(nc_enum_types_headers)  \
	$(nc_enum_types_sources)

nobase_cosmo_include_HEADERS = \
	$(ncm_headers)                 \
	$(ncm_headers_internal)        \
	$(nc_headers)                  \
	$(ncm_cond_headers)            \
	math/binsplit_eval.c           \
	perturbations/linear_generic.c \
	numcosmo-math.h                \
	numcosmo.h

nodist_cosmo_include_HEADERS = \
	build_cfg.h               \
	$(ncm_enum_types_headers) \
	$(nc_enum_types_headers)

EXTRA_DIST = \
	$(nobase_cosmo_include_HEADERS) \
	build_cfg.h.in                  \
	ncm_enum_types.h.template       \
	ncm_enum_types.c.template       \
	ncm_fit_nlopt_enum.h.template   \
	ncm_fit_nlopt_enum.c.template   \
	nc_enum_types.h.template        \
	nc_enum_types.c.template        \
	ncm_gir_scan.h                  \
	ncm_fit_nlopt_enum_meta.c

CLEANFILES = $(BUILT_SOURCES)

-include $(INTROSPECTION_MAKEFILE)
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
NumCosmoMath_1_0_gir_CFLAGS       = $(INCLUDES) -DNUMCOSMO_GIR_SCAN
NumCosmoMath_1_0_gir_LIBS         = libnumcosmo.la
NumCosmoMath_1_0_gir_FILES        = $(ncm_introspection_sources)
NumCosmoMath_1_0_gir_SCANNERFLAGS = \
	--c-include="numcosmo/numcosmo-math.h"  \
	-n NumCosmoMath --identifier-prefix=Ncm

INTROSPECTION_GIRS += NumCosmoMath-1.0.gir

NumCosmo-1.0.gir: libnumcosmo.la $(INTROSPECTION_SCANNER) NumCosmoMath-1.0.gir
NumCosmo_1_0_gir_VERSION      = 1.0
NumCosmo_1_0_gir_INCLUDES     = GObject-2.0
NumCosmo_1_0_gir_CFLAGS       = $(INCLUDES) -DNUMCOSMO_GIR_SCAN
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

endif

endif
