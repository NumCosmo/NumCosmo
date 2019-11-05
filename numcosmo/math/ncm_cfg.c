/***************************************************************************
 *            ncm_cfg.c
 *
 *  Wed Aug 13 20:59:22 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
 * numcosmo is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * numcosmo is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * SECTION:ncm_cfg
 * @title: NcmCfg
 * @short_description: Library configuration and helper functions.
 *
 * FIXME
 *
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "math/ncm_cfg.h"
#include "math/ncm_rng.h"
#include "math/ncm_memory_pool.h"
#include "math/ncm_mpi_job.h"
#include "math/ncm_mpi_job_test.h"
#include "math/ncm_mpi_job_fit.h"
#include "math/ncm_mpi_job_mcmc.h"
#include "math/ncm_vector.h"
#include "math/ncm_spline_gsl.h"
#include "math/ncm_spline_cubic.h"
#include "math/ncm_spline_cubic_notaknot.h"
#include "math/ncm_spline2d_bicubic.h"
#include "math/ncm_spline2d_gsl.h"
#include "math/ncm_spline2d_spline.h"
#include "math/ncm_powspec.h"
#include "math/ncm_powspec_filter.h"
#include "math/ncm_powspec_corr3d.h"
#include "math/ncm_model.h"
#include "math/ncm_model_ctrl.h"
#include "math/ncm_model_builder.h"
#include "math/ncm_model_mvnd.h"
#include "math/ncm_reparam_linear.h"
#include "math/ncm_data.h"
#include "math/ncm_data_gauss_cov_mvnd.h"
#include "math/ncm_stats_vec.h"
#include "math/ncm_fit_esmcmc_walker_stretch.h"
#include "math/ncm_data.h"
#include "math/ncm_stats_dist1d_epdf.h"
#include "math/ncm_stats_dist1d_spline.h"
#include "math/ncm_dataset.h"
#include "math/ncm_fit.h"
#include "math/ncm_fit_nlopt.h"
#include "math/ncm_prior_gauss_param.h"
#include "math/ncm_prior_gauss_func.h"
#include "nc_hicosmo.h"
#include "nc_cbe_precision.h"
#include "model/nc_hicosmo_qconst.h"
#include "model/nc_hicosmo_qlinear.h"
#include "model/nc_hicosmo_qspline.h"
#include "model/nc_hicosmo_qrbf.h"
#include "model/nc_hicosmo_lcdm.h"
#include "model/nc_hicosmo_gcg.h"
#include "model/nc_hicosmo_idem2.h"
#include "model/nc_hicosmo_de_xcdm.h"
#include "model/nc_hicosmo_de_cpl.h"
#include "model/nc_hicosmo_de_jbp.h"
#include "model/nc_hicosmo_qgrw.h"
#include "model/nc_hicosmo_Vexp.h"
#include "model/nc_hicosmo_de_reparam_ok.h"
#include "model/nc_hicosmo_de_reparam_cmb.h"
#include "model/nc_hiprim_power_law.h"
#include "model/nc_hiprim_atan.h"
#include "model/nc_hiprim_expc.h"
#include "model/nc_hiprim_bpl.h"
#include "model/nc_hiprim_sbpl.h"
#include "lss/nc_window_tophat.h"
#include "lss/nc_window_gaussian.h"
#include "lss/nc_growth_func.h"
#include "lss/nc_transfer_func.h"
#include "lss/nc_transfer_func_bbks.h"
#include "lss/nc_transfer_func_eh.h"
#include "lss/nc_transfer_func_camb.h"
#include "lss/nc_density_profile.h"
#include "lss/nc_density_profile_nfw.h"
#include "lss/nc_density_profile_einasto.h"
#include "lss/nc_density_profile_dk14.h"
#include "lss/nc_multiplicity_func.h"
#include "lss/nc_multiplicity_func_st.h"
#include "lss/nc_multiplicity_func_ps.h"
#include "lss/nc_multiplicity_func_jenkins.h"
#include "lss/nc_multiplicity_func_warren.h"
#include "lss/nc_multiplicity_func_tinker.h"
#include "lss/nc_multiplicity_func_tinker_mean.h"
#include "lss/nc_multiplicity_func_tinker_crit.h"
#include "lss/nc_multiplicity_func_tinker_mean_normalized.h"
#include "lss/nc_multiplicity_func_crocce.h"
#include "lss/nc_halo_mass_function.h"
#include "lss/nc_galaxy_acf.h"
#include "lss/nc_galaxy_redshift_spec.h"
#include "lss/nc_galaxy_redshift_spline.h"
#include "lss/nc_galaxy_acf.h"
#include "lss/nc_cluster_mass.h"
#include "lss/nc_cluster_mass_nodist.h"
#include "lss/nc_cluster_mass_lnnormal.h"
#include "lss/nc_cluster_mass_vanderlinde.h"
#include "lss/nc_cluster_mass_benson.h"
#include "lss/nc_cluster_mass_benson_xray.h"
#include "lss/nc_cluster_mass_plcl.h"
#include "lss/nc_cluster_mass_ascaso.h"
#include "lss/nc_cluster_redshift.h"
#include "lss/nc_cluster_redshift_nodist.h"
#include "lss/nc_cluster_photoz_gauss_global.h"
#include "lss/nc_cluster_photoz_gauss.h"
#include "lss/nc_halo_bias_func.h"
#include "lss/nc_halo_bias_type_ps.h"
#include "lss/nc_halo_bias_type_st_ellip.h"
#include "lss/nc_halo_bias_type_st_spher.h"
#include "lss/nc_halo_bias_type_tinker.h"
#include "lss/nc_cluster_abundance.h"
#include "lss/nc_cluster_pseudo_counts.h"
#include "lss/nc_cor_cluster_cmb_lens_limber.h"
#include "lss/nc_wl_surface_mass_density.h"
#include "lss/nc_reduced_shear_cluster_mass.h"
#include "lss/nc_reduced_shear_calib.h"
#include "lss/nc_reduced_shear_calib_wtg.h"
#include "nc_distance.h"
#include "nc_recomb.h"
#include "nc_recomb_cbe.h"
#include "nc_recomb_seager.h"
#include "nc_hireion.h"
#include "nc_hireion_camb.h"
#include "nc_powspec_ml.h"
#include "nc_powspec_ml_transfer.h"
#include "nc_powspec_ml_cbe.h"
#include "nc_powspec_mnl.h"
#include "nc_powspec_mnl_halofit.h"
#include "nc_snia_dist_cov.h"
#include "nc_planck_fi.h"
#include "nc_planck_fi_cor_tt.h"
#include "nc_planck_fi_cor_ttteee.h"
#include "perturbations/nc_hipert_boltzmann_cbe.h"
#include "data/nc_data_bao_a.h"
#include "data/nc_data_bao_dv.h"
#include "data/nc_data_bao_dvdv.h"
#include "data/nc_data_bao_rdv.h"
#include "data/nc_data_bao_empirical_fit.h"
#include "data/nc_data_bao_empirical_fit_2d.h"
#include "data/nc_data_bao_dhr_dar.h"
#include "data/nc_data_bao_dmr_hr.h"
#include "data/nc_data_dist_mu.h"
#include "data/nc_data_cluster_pseudo_counts.h"
#include "data/nc_data_cluster_counts_box_poisson.h"
#include "data/nc_data_reduced_shear_cluster_mass.h"
#include "data/nc_data_cmb_shift_param.h"
#include "data/nc_data_cmb_dist_priors.h"
#include "data/nc_data_hubble.h"
#include "data/nc_data_snia_cov.h"
#include "data/nc_data_xcor.h"
#include "data/nc_data_planck_lkl.h"
#include "xcor/nc_xcor.h"
#include "xcor/nc_xcor_AB.h"
#include "xcor/nc_xcor_limber_kernel.h"
#include "xcor/nc_xcor_limber_kernel_gal.h"
#include "xcor/nc_xcor_limber_kernel_CMB_lensing.h"
#include "xcor/nc_xcor_limber_kernel_weak_lensing.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <stdlib.h>
#include <gio/gio.h>
#ifdef NUMCOSMO_HAVE_FFTW3
#include <fftw3.h>
#endif /* NUMCOSMO_HAVE_FFTW3 */
#include <cuba.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */

#ifndef G_VALUE_INIT
#define G_VALUE_INIT {0}
#endif

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif /* HAVE_EXECINFO_H */
#endif /* NUMCOSMO_GIR_SCAN */

static gchar *numcosmo_path = NULL;
static gboolean numcosmo_init = FALSE;
static FILE *_log_stream = NULL;
static FILE *_log_stream_err = NULL;
static guint _log_msg_id = 0;
static guint _log_err_id = 0;
static gboolean _enable_msg = TRUE;
static gboolean _enable_msg_flush = TRUE;
static gsl_error_handler_t *gsl_err = NULL;

static void
_ncm_cfg_log_message (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  NCM_UNUSED (log_domain);
  NCM_UNUSED (log_level);
  NCM_UNUSED (user_data);
  if (_enable_msg && _log_stream)
  {
    fprintf (_log_stream, "%s", message);
    if (_enable_msg_flush)
      fflush (_log_stream);
  }
}

static void
_ncm_cfg_log_error (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  const gchar *pname = g_get_prgname ();
  NCM_UNUSED (log_domain);
  NCM_UNUSED (log_level);
  NCM_UNUSED (user_data);

  fprintf (_log_stream_err, "# (%s): %s-ERROR: %s\n", pname, log_domain, message);
#if defined (HAVE_BACKTRACE) && defined (HAVE_BACKTRACE_SYMBOLS)
  {
    gpointer tarray[30];
    gsize size = backtrace (tarray, 30);
    gchar **trace = backtrace_symbols (tarray, size);
    gsize i;
    // print out all the frames to stderr
    for (i = 0; i < size; i++)
    {
      fprintf (_log_stream_err, "# (%s): %s-BACKTRACE:[%02zd] %s\n", pname, log_domain, i, trace[i]);
    }
    g_free (trace);
  }
#endif
  fflush (_log_stream_err);

  abort ();
}

void clencurt_gen (int M);

#ifdef HAVE_OPENBLAS_SET_NUM_THREADS
  void openblas_set_num_threads (gint);
#endif /* HAVE_OPENBLAS_SET_NUM_THREADS */

#ifdef HAVE_MKL_SET_NUM_THREADS
  void MKL_Set_Num_Threads (gint);
#endif /* HAVE_MKL_SET_NUM_THREADS */

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

void _nc_hicosmo_register_functions (void);
void _nc_hicosmo_de_register_functions (void);
void _nc_hiprim_register_functions (void);
void _nc_hireion_register_functions (void);
void _nc_distance_register_functions (void);
void _nc_planck_fi_cor_tt_register_functions (void);

#ifdef HAVE_MPI
static void _ncm_cfg_mpi_main_loop (void);
#endif /* HAVE_MPI */

NcmMPIJobCtrl _mpi_ctrl;

void 
ncm_cfg_mpi_kill_all_slaves (void)
{
#ifdef HAVE_MPI
	if (_mpi_ctrl.rank != NCM_MPI_CTRL_MASTER_ID)
		return;

	if (_mpi_ctrl.size > 1)
	{
		gint i;

		for (i = 0; i < _mpi_ctrl.nslaves; i++)
		{
			gint slave_id = i + 1;
			gint cmd      = NCM_MPI_CTRL_SLAVE_KILL;
			
			MPI_Send (&cmd, 1, MPI_INT, slave_id, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD);	
		}

		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All slaves killed!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
	}
#else
#endif /* HAVE_MPI */
}

static void
_ncm_cfg_exit (void)
{
#ifdef HAVE_MPI
	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Dying [%d]!\n", _mpi_ctrl.size, _mpi_ctrl.rank, _mpi_ctrl.initialized);
	if (_mpi_ctrl.initialized)
	{
		ncm_cfg_mpi_kill_all_slaves ();
		MPI_Barrier (MPI_COMM_WORLD);
		MPI_Finalize ();
	}
#endif /* HAVE_MPI */
}

/**
 * ncm_cfg_init:
 *
 * Main library configuration function. Must be called before any other function of NumCosmo.
 *
 * Initializes internal variables and sets all other library number of threads to one.
 *
 * See also: ncm_cfg_init_full() ncm_cfg_init_full_ptr().
 */
void
ncm_cfg_init (void)
{
	ncm_cfg_init_full (0, NULL);
}

static gchar **
_ncm_cfg_make_strv (gint argc, gchar **argv)
{
	if ((argc == 0) || (argv == NULL))
		return NULL;
	else
	{
		gchar **argv_dup = g_new (gchar *, argc + 1);
		gint i;
		for (i = 0; i < argc; i++)
		{
			argv_dup[i] = g_strdup (argv[i]);
		}
		argv_dup[i] = NULL;

		return argv_dup;
	}
}

/**
 * ncm_cfg_init_full:
 * @argc: a pointer to argc
 * @argv: (array length=argc): a pointer to argv
 *
 * Main library configuration function. Must be called before any other function of NumCosmo.
 *
 * Initializes internal variables and sets all other library number of threads to one.
 * This function passes the arguments to other libraries, e.g, MPI. If that is not necessary
 * ncm_cfg_init() should be used. This version is compatible with bindings and can be safely
 * called from other languages.
 * 
 * See also: ncm_cfg_init() ncm_cfg_init_full_ptr().
 * 
 * Returns: (transfer full) (array zero-terminated=1): the possibly modified argv
 */
gchar **
ncm_cfg_init_full (gint argc, gchar **argv)
{
	gchar **argv1 = _ncm_cfg_make_strv (argc, argv);
	gchar **argv2 = argv1;
	gchar **argv_ret;
	
	ncm_cfg_init_full_ptr (&argc, &argv1);

	argv_ret = _ncm_cfg_make_strv (argc, argv1);
	g_strfreev (argv2);
	
	return argv_ret;
}

/**
 * ncm_cfg_init_full_ptr:
 * @argc: a pointer to argc
 * @argv: (array length=argc): a pointer to argv
 *
 * Main library configuration function. Must be called before any other function of NumCosmo.
 *
 * Initializes internal variables and sets all other library number of threads to one.
 * This function passes the arguments to other libraries, e.g, MPI. If that is not necessary
 * ncm_cfg_init() should be used. This version should be used from C applications passing
 * @argc and @argv pointers from main.
 * 
 * See also: ncm_cfg_init() ncm_cfg_init_full().
 */
void
ncm_cfg_init_full_ptr (gint *argc, gchar ***argv)
{
  const gchar *home;
  if (numcosmo_init)
    return;

  if (sizeof (NcmComplex) != sizeof (fftw_complex))
    g_warning ("NcmComplex is not binary compatible with complex double, expect problems with it!");

  home = g_get_home_dir ();
  numcosmo_path = g_build_filename (home, ".numcosmo", NULL);
  if (!g_file_test (numcosmo_path, G_FILE_TEST_EXISTS))
    g_mkdir_with_parents (numcosmo_path, 0755);

#ifdef HAVE_OPENBLAS_SET_NUM_THREADS
  openblas_set_num_threads (1);
#endif /* HAVE_OPENBLAS_SET_NUM_THREADS */

#ifdef HAVE_MKL_SET_NUM_THREADS
  MKL_Set_Num_Threads (1);
#endif /* HAVE_MKL_SET_NUM_THREADS */

#ifdef _OPENMP
  omp_set_num_threads (1);
#endif /* _OPENMP */

  g_setenv ("CUBACORES", "0", TRUE);
  g_setenv ("CUBACORESMAX", "0", TRUE);
  g_setenv ("CUBAACCEL", "0", TRUE);
  g_setenv ("CUBAACCELMAX", "0", TRUE);
#ifdef HAVE_LIBCUBA_4_0
  cubaaccel (0, 0);
  cubacores (0, 0);
#endif

  gsl_err = gsl_set_error_handler_off ();
  
#ifdef NUMCOSMO_HAVE_FFTW3
  fftw_set_timelimit (10.0);
#endif /* NUMCOSMO_HAVE_FFTW3 */
#ifdef HAVE_FFTW3F
  fftwf_set_timelimit (10.0);
#endif /* HAVE_FFTW3F */
  
#if !GLIB_CHECK_VERSION(2,36,0)
  g_type_init ();
#endif

  _log_stream = stdout;
  _log_stream_err = stderr;

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message, NULL);
  _log_err_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, _ncm_cfg_log_error, NULL);

  ncm_cfg_register_obj (NCM_TYPE_RNG);

  ncm_cfg_register_obj (NCM_TYPE_VECTOR);
  ncm_cfg_register_obj (NCM_TYPE_MATRIX);

  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_TEST);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_FIT);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_MCMC);

  ncm_cfg_register_obj (NCM_TYPE_SPLINE);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC_NOTAKNOT);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_GSL);

  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_BICUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_GSL);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_SPLINE);

  ncm_cfg_register_obj (NCM_TYPE_POWSPEC);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_FILTER);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_CORR3D);

  ncm_cfg_register_obj (NCM_TYPE_MODEL);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_CTRL);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_BUILDER);

  ncm_cfg_register_obj (NCM_TYPE_MODEL_MVND);

  ncm_cfg_register_obj (NCM_TYPE_REPARAM);
  ncm_cfg_register_obj (NCM_TYPE_REPARAM_LINEAR);

  ncm_cfg_register_obj (NCM_TYPE_BOOTSTRAP);
  ncm_cfg_register_obj (NCM_TYPE_STATS_VEC);

  ncm_cfg_register_obj (NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH);

	ncm_cfg_register_obj (NCM_TYPE_DATA);
	ncm_cfg_register_obj (NCM_TYPE_DATASET);

  ncm_cfg_register_obj (NCM_TYPE_DATA_GAUSS_COV_MVND);
  
	ncm_cfg_register_obj (NCM_TYPE_FIT);
	ncm_cfg_register_obj (NCM_TYPE_FIT_NLOPT);

	ncm_cfg_register_obj (NCM_TYPE_PRIOR_GAUSS_PARAM);
	ncm_cfg_register_obj (NCM_TYPE_PRIOR_GAUSS_FUNC);

  ncm_cfg_register_obj (NCM_TYPE_DATA);

	ncm_cfg_register_obj (NCM_TYPE_STATS_DIST1D_EPDF);
	ncm_cfg_register_obj (NCM_TYPE_STATS_DIST1D_SPLINE);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QCONST);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QLINEAR);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QSPLINE);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QRBF);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_LCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_GCG);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_IDEM2);  
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_XCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_CPL);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_JBP);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QGRW);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_VEXP);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_REPARAM_OK);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_REPARAM_CMB);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_GCG_REPARAM_OK);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_GCG_REPARAM_CMB);
  
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_IDEM2_REPARAM_OK);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_IDEM2_REPARAM_CMB);
    
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_POWER_LAW);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_ATAN);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_EXPC);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_BPL);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_SBPL);

  ncm_cfg_register_obj (NC_TYPE_CBE_PRECISION);

  ncm_cfg_register_obj (NC_TYPE_WINDOW);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_TOPHAT);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_GAUSSIAN);

  ncm_cfg_register_obj (NC_TYPE_GROWTH_FUNC);

  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_BBKS);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_EH);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_CAMB);

  ncm_cfg_register_obj (NC_TYPE_DENSITY_PROFILE);
  ncm_cfg_register_obj (NC_TYPE_DENSITY_PROFILE_NFW);
  ncm_cfg_register_obj (NC_TYPE_DENSITY_PROFILE_EINASTO);
  ncm_cfg_register_obj (NC_TYPE_DENSITY_PROFILE_DK14);

  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_PS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_ST);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_JENKINS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_WARREN);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_CRIT);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_CROCCE);	

  ncm_cfg_register_obj (NC_TYPE_HALO_MASS_FUNCTION);

  ncm_cfg_register_obj (NC_TYPE_GALAXY_ACF);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_SPEC);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_SPLINE);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_LNNORMAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_VANDERLINDE);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON_XRAY);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_PLCL);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS);

  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_FUNC);

  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_PS);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_ST_ELLIP);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_ST_SPHER);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TYPE_TINKER);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_ABUNDANCE);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PSEUDO_COUNTS);

  ncm_cfg_register_obj (NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER);

  ncm_cfg_register_obj (NC_TYPE_WL_SURFACE_MASS_DENSITY);

  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS);
  
  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CALIB);
  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CALIB_WTG);

  ncm_cfg_register_obj (NC_TYPE_DISTANCE);

  ncm_cfg_register_obj (NC_TYPE_RECOMB);
  ncm_cfg_register_obj (NC_TYPE_RECOMB_CBE);
  ncm_cfg_register_obj (NC_TYPE_RECOMB_SEAGER);

  ncm_cfg_register_obj (NC_TYPE_HIREION);
  ncm_cfg_register_obj (NC_TYPE_HIREION_CAMB);

  ncm_cfg_register_obj (NC_TYPE_POWSPEC_ML);
  ncm_cfg_register_obj (NC_TYPE_POWSPEC_ML_TRANSFER);
  ncm_cfg_register_obj (NC_TYPE_POWSPEC_ML_CBE);

  ncm_cfg_register_obj (NC_TYPE_POWSPEC_MNL);
  ncm_cfg_register_obj (NC_TYPE_POWSPEC_MNL_HALOFIT);

  ncm_cfg_register_obj (NC_TYPE_SNIA_DIST_COV);

  ncm_cfg_register_obj (NC_TYPE_PLANCK_FI);
  ncm_cfg_register_obj (NC_TYPE_PLANCK_FI_COR_TT);
  ncm_cfg_register_obj (NC_TYPE_PLANCK_FI_COR_TTTEEE);

	ncm_cfg_register_obj (NC_TYPE_HIPERT_BOLTZMANN_CBE);
	
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_A);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DV);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DVDV);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_RDV);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_EMPIRICAL_FIT);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_EMPIRICAL_FIT_2D);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DHR_DAR);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DMR_HR);

  ncm_cfg_register_obj (NC_TYPE_DATA_DIST_MU);

  ncm_cfg_register_obj (NC_TYPE_DATA_HUBBLE);

	ncm_cfg_register_obj (NC_TYPE_DATA_SNIA_COV);

  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_COUNTS_BOX_POISSON);
  ncm_cfg_register_obj (NC_TYPE_DATA_REDUCED_SHEAR_CLUSTER_MASS);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS);

  ncm_cfg_register_obj (NC_TYPE_DATA_CMB_SHIFT_PARAM);
  ncm_cfg_register_obj (NC_TYPE_DATA_CMB_DIST_PRIORS);

  ncm_cfg_register_obj (NC_TYPE_XCOR);
  ncm_cfg_register_obj (NC_TYPE_XCOR_LIMBER_KERNEL);
  ncm_cfg_register_obj (NC_TYPE_XCOR_LIMBER_KERNEL_GAL);
  ncm_cfg_register_obj (NC_TYPE_XCOR_LIMBER_KERNEL_CMB_LENSING);
  ncm_cfg_register_obj (NC_TYPE_XCOR_LIMBER_KERNEL_WEAK_LENSING);
  ncm_cfg_register_obj (NC_TYPE_DATA_XCOR);
  ncm_cfg_register_obj (NC_TYPE_XCOR_AB);

	ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_LKL);

  _nc_hicosmo_register_functions ();
  _nc_hicosmo_de_register_functions ();
  _nc_hiprim_register_functions ();
  _nc_hireion_register_functions ();
  _nc_distance_register_functions ();
  _nc_planck_fi_cor_tt_register_functions ();

  numcosmo_init = TRUE;

	_mpi_ctrl.initialized    = 0;
	_mpi_ctrl.size           = 1;
	_mpi_ctrl.rank           = 0;
	_mpi_ctrl.nslaves        = 0;
	_mpi_ctrl.working_slaves = 0;

	atexit (_ncm_cfg_exit);

#ifdef HAVE_MPI
	MPI_Initialized (&_mpi_ctrl.initialized);
	if (!_mpi_ctrl.initialized)
	{
		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] MPI not initalized, calling MPI_Init.\n", _mpi_ctrl.size, _mpi_ctrl.rank);
		MPI_Init (argc, argv);
		MPI_Initialized (&_mpi_ctrl.initialized);
	}
	else
	{
		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] MPI was already initalized!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
	}
  {
    gchar mpi_hostname[MPI_MAX_PROCESSOR_NAME];
    gint len = 0;
    
    MPI_Comm_size (MPI_COMM_WORLD, &_mpi_ctrl.size);
    MPI_Comm_rank (MPI_COMM_WORLD, &_mpi_ctrl.rank);
    MPI_Get_processor_name (mpi_hostname, &len);

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] We have %d mpi process!! My rank is %d and I'm running on `%s'.\n", _mpi_ctrl.size, _mpi_ctrl.rank, _mpi_ctrl.size, _mpi_ctrl.rank, mpi_hostname);

    if (_mpi_ctrl.rank != NCM_MPI_CTRL_MASTER_ID)
    {
      _ncm_cfg_mpi_main_loop ();
    }
		else
		{
			_mpi_ctrl.nslaves        = (_mpi_ctrl.size - 1);
			_mpi_ctrl.working_slaves = 0;
		}
  }
#endif /* HAVE_MPI */

  return;
}

#ifdef HAVE_MPI

static gboolean _ncm_cfg_mpi_cmd_handler (gpointer user_data);

static void
_ncm_cfg_mpi_main_loop (void)
{
  GMainLoop *mpi_ml = g_main_loop_new (NULL, FALSE);
  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Starting slave!\n", _mpi_ctrl.size, _mpi_ctrl.rank);

  g_timeout_add (100, &_ncm_cfg_mpi_cmd_handler, mpi_ml);
  
  g_main_loop_run (mpi_ml);
  
  g_main_loop_unref (mpi_ml);

	NCM_MPI_JOB_DEBUG_PRINT ("#[%d %d] Dying slave!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
  exit (0);
}

static gboolean 
_ncm_cfg_mpi_cmd_handler (gpointer user_data)
{
	enum buf_type { input_type, ret_type, msg_type, };
	struct buf_desc { gpointer obj; gpointer buf; enum buf_type t; };
	GMainLoop *mpi_ml        = user_data;
	NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_CLEAN_DUP);
	NcmMPIJob *mpi_job       = NULL;
	gboolean init            = FALSE;
	GArray *input_array      = g_array_new (FALSE, FALSE, sizeof (gdouble));
	GArray *work_ret_request = g_array_new (FALSE, FALSE, sizeof (MPI_Request));
	GArray *work_ret_bufs    = g_array_new (FALSE, TRUE, sizeof (struct buf_desc));
	gpointer input           = NULL;
	gpointer input_buf       = NULL;
	gint input_len           = 0;
	gint input_size          = 0;
	gint return_len          = 0;
	gint return_size         = 0;
	gboolean normal_exit     = TRUE;
	MPI_Datatype input_dtype;
	MPI_Datatype return_dtype;

	while (TRUE)
	{
		gboolean end  = FALSE;
		gint cmd      = 0;
		MPI_Status status;

		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Waiting for command...\n", _mpi_ctrl.size, _mpi_ctrl.rank);

		MPI_Recv (&cmd, 1, MPI_INT, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_CMD, MPI_COMM_WORLD, &status);

		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Received %d\n", _mpi_ctrl.size, _mpi_ctrl.rank, cmd);

		switch (cmd)
		{
			case NCM_MPI_CTRL_SLAVE_INIT:
			{
				GVariant *job_ser = NULL;
				gint job_size     = 0;
				gint job_recv     = 0;
				gchar *job        = NULL;

				if (init)
					g_error ("_ncm_cfg_mpi_cmd_handler: slave %d already initialized.", _mpi_ctrl.rank);

				NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Initializing slave.\n", _mpi_ctrl.size, _mpi_ctrl.rank);

				MPI_Probe (NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &status);
				MPI_Get_count (&status, MPI_BYTE, &job_size);

				NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave object size %d.\n", _mpi_ctrl.size, _mpi_ctrl.rank, job_size);
			
				job = g_new (gchar, job_size);
				MPI_Recv (job, job_size, MPI_BYTE, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_JOB, MPI_COMM_WORLD, &status);
				MPI_Get_count (&status, MPI_BYTE, &job_recv);

				NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave object received size %d.\n", _mpi_ctrl.size, _mpi_ctrl.rank, job_recv);

				g_assert_cmpint (job_recv, ==, job_size);
			
				job_ser = g_variant_new_from_data (G_VARIANT_TYPE (NCM_SERIALIZE_OBJECT_TYPE), job, job_size, TRUE, g_free, job);

				/*NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave object received string `%s'.\n", _mpi_ctrl.size, _mpi_ctrl.rank, g_variant_print (job_ser, TRUE));*/
			
				mpi_job = NCM_MPI_JOB (ncm_serialize_from_variant (ser, job_ser));

				ncm_mpi_job_work_init (mpi_job);

				input_dtype  = ncm_mpi_job_input_datatype  (mpi_job, &input_len,  &input_size);
				return_dtype = ncm_mpi_job_return_datatype (mpi_job, &return_len, &return_size);
				input        = ncm_mpi_job_create_input (mpi_job);
				input_buf    = ncm_mpi_job_get_input_buffer (mpi_job, input);
				
				g_assert (NCM_IS_MPI_JOB (mpi_job));

				g_variant_unref (job_ser);
				
				init = TRUE;
				break;
			}
			case NCM_MPI_CTRL_SLAVE_FREE:
				end = TRUE;
				break;
			case NCM_MPI_CTRL_SLAVE_KILL:
				end         = TRUE;
				normal_exit = FALSE;
				break;
			case NCM_MPI_CTRL_SLAVE_WORK:
			{
				if (!init)
				{
					g_error ("_ncm_cfg_mpi_cmd_handler: uninitialized slave `%d' received work (vector).", _mpi_ctrl.rank);
				}
				else
				{
					struct buf_desc bd = {NULL, NULL, ret_type};
					gint input_recv    = 0;
					MPI_Request wr_request;

					MPI_Recv (input_buf, input_len, input_dtype, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD, &status);
					MPI_Get_count (&status, input_dtype, &input_recv);

					g_assert_cmpint (input_recv, ==, input_len);

					ncm_mpi_job_unpack_input (mpi_job, input_buf, input);

					bd.obj = ncm_mpi_job_create_return (mpi_job);
					
					ncm_mpi_job_run (mpi_job, input, bd.obj);
					bd.buf = ncm_mpi_job_pack_return (mpi_job, bd.obj);
					
					NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Returning!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
					
					MPI_Isend (bd.buf, return_len, return_dtype, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_WORK_RETURN, MPI_COMM_WORLD, &wr_request);
					
					g_array_append_val (work_ret_request, wr_request);
					g_array_append_val (work_ret_bufs,    bd);
				}
				break;
			}
			default:
				g_error ("_ncm_cfg_mpi_cmd_handler: unknown MPI message `%d' to slave %d", cmd, _mpi_ctrl.rank);
				break;
		}

		if (work_ret_request->len > 0)
		{
			gint i;
			NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Testing %d sends:\n", _mpi_ctrl.size, _mpi_ctrl.rank, work_ret_request->len);
			for (i = work_ret_request->len - 1; i >= 0; i--)
			{
				gint done = 0;
				MPI_Test (&g_array_index (work_ret_request, MPI_Request, i), &done, &status);
				NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Send %d is %s!\n", _mpi_ctrl.size, _mpi_ctrl.rank, i, done ? "done" : "not done");
				if (done)
				{
					struct buf_desc bd = g_array_index (work_ret_bufs, struct buf_desc, i);
					ncm_mpi_job_destroy_return_buffer (mpi_job, bd.obj, bd.buf);
					ncm_mpi_job_destroy_return (mpi_job, bd.obj);
					
					g_array_remove_index_fast (work_ret_request, i);
					g_array_remove_index_fast (work_ret_bufs, i);
				}
			}
		}

		NCM_MPI_JOB_DEBUG_PRINT	("#[%3d %3d] Finished command, %d requests left, %d%s\n", _mpi_ctrl.size, _mpi_ctrl.rank, work_ret_request->len, cmd, end ? ", exiting!" : ".");

		if (end)
			break;
	}

	if (work_ret_request->len > 0)
	{
		gint i;
		MPI_Waitall (work_ret_request->len, (MPI_Request *)work_ret_request->data, MPI_STATUSES_IGNORE);

		NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] All sent, freeing %d buffers!\n", _mpi_ctrl.size, _mpi_ctrl.rank, work_ret_request->len);

		for (i = 0; i < work_ret_bufs->len; i++)
		{
			struct buf_desc bd = g_array_index (work_ret_bufs, struct buf_desc, i);
			ncm_mpi_job_destroy_return_buffer (mpi_job, bd.obj, bd.buf);
			ncm_mpi_job_destroy_return (mpi_job, bd.obj);
		}		
	}

	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Freeing arrays!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
	
	g_array_unref (input_array);
	g_array_unref (work_ret_request);
	g_array_unref (work_ret_bufs);

	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Freeing input buffer %p [input %p, mpi_job %p]!\n", _mpi_ctrl.size, _mpi_ctrl.rank, input_buf, input, mpi_job);

	if (input_buf != NULL)
		ncm_mpi_job_destroy_input_buffer (mpi_job, input, input_buf);
	
	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Freeing input %p [mpi_job %p]!\n", _mpi_ctrl.size, _mpi_ctrl.rank, input, mpi_job);

	if (input != NULL)
		ncm_mpi_job_destroy_input (mpi_job, input);
	
	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Clearing MPI job %p!\n", _mpi_ctrl.size, _mpi_ctrl.rank, mpi_job);

	ncm_mpi_job_clear (&mpi_job);

	NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Returning %d!\n", _mpi_ctrl.size, _mpi_ctrl.rank, normal_exit);

	if (!normal_exit)
		g_main_loop_quit (mpi_ml);

	return normal_exit;
}

#endif /* HAVE_MPI */

/** 
 * ncm_cfg_enable_gsl_err_handler:
 *
 * FIXME
 */
void
ncm_cfg_enable_gsl_err_handler (void)
{
  g_assert (numcosmo_init);
  gsl_set_error_handler (gsl_err);
}

static uint nreg_model = 0;

/**
 * ncm_cfg_register_obj:
 * @obj: FIXME
 *
 * FIXME
 */
void
ncm_cfg_register_obj (GType obj)
{
#if GLIB_CHECK_VERSION(2,34,0)
  g_type_ensure (obj);
#endif /* GLIB >= 2.34*/
  gpointer obj_class = g_type_class_ref (obj);
  g_type_class_unref (obj_class);
  nreg_model++;
}

/**
 * ncm_cfg_mpi_nslaves:
 * 
 * Returns: the total number of available slaves.
 */ 
guint 
ncm_cfg_mpi_nslaves (void)
{
	return _mpi_ctrl.nslaves;
}

/**
 * ncm_cfg_set_logfile:
 * @filename: name of the log-file
 *
 * Sets all log information to @filename.
 */
void
ncm_cfg_set_logfile (gchar *filename)
{
  FILE *out = g_fopen (filename, "w");

  if (out != NULL)
    _log_stream = out;
  else
    g_error ("ncm_cfg_set_logfile: Can't open logfile `%s' %s", filename, g_strerror (errno));
}

/**
 * ncm_cfg_set_logstream:
 * @stream: a stream
 *
 * Sets all log information to @stream.
 */
void
ncm_cfg_set_logstream (FILE *stream)
{
  g_assert (stream != NULL);
  _log_stream = stream;
}

static void
_ncm_cfg_log_message_logger (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  NCM_UNUSED (log_domain);
  NCM_UNUSED (log_level);
  NCM_UNUSED (user_data);
  if (_enable_msg && _log_stream)
  {
    void (*logger) (const gchar *msg) = user_data;
    logger (message);
  }
}

/**
 * ncm_cfg_set_log_handler:
 * @logger: (scope notified): a logger function
 *
 * Sets all log information to @stream.
 */
void
ncm_cfg_set_log_handler (NcmCfgLoggerFunc logger)
{
  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message_logger, logger);
}

/**
 * ncm_cfg_set_error_log_handler:
 * @logger: (scope notified): a logger function
 *
 * Sets all log information to @stream.
 */
void
ncm_cfg_set_error_log_handler (NcmCfgLoggerFunc logger)
{
  _log_err_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, _ncm_cfg_log_message_logger, logger);
}

/**
 * ncm_cfg_logfile:
 * @on: FIXME
 *
 * FIXME
 */
void
ncm_cfg_logfile (gboolean on) 
{ 
  _enable_msg = on; 
}

/**
 * ncm_cfg_logfile_flush:
 * @on: FIXME
 *
 * FIXME
 */
void
ncm_cfg_logfile_flush (gboolean on) 
{ 
  _enable_msg_flush = on; 
}

/**
 * ncm_cfg_logfile_flush_now:
 *
 * FIXME
 */
void
ncm_cfg_logfile_flush_now (void) 
{
  fflush (_log_stream); 
}

/**
 * ncm_message:
 * @msg: FIXME
 * @...: FIXME
 *
 * FIXME
 */
void
ncm_message (const gchar *msg, ...)
{
  va_list ap;
  va_start(ap, msg);
  g_logv (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, msg, ap);
  va_end (ap);
}

/**
 * ncm_string_ww:
 * @msg: FIXME
 * @first: FIXME
 * @rest: FIXME
 * @ncols: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): word wraped string @msg.
 */
gchar *
ncm_string_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols)
{
  gchar **msg_split = g_strsplit (msg, " ", 0);
  guint size_print = strlen (msg_split[0]);
  guint first_size = strlen (first);
  guint rest_size = strlen (rest);
  guint msg_len = strlen (msg);
  guint msg_final_len = msg_len + first_size + msg_len / ncols * rest_size;
  GString *msg_ww = g_string_sized_new (msg_final_len);
  guint i = 1;

  g_string_append_printf (msg_ww, "%s%s", first, msg_split[0]);
  while (msg_split[i] != NULL)
  {
    guint lsize = strlen (msg_split[i]);
    if (size_print + lsize + 1 > (ncols - first_size))
      break;
    g_string_append_printf (msg_ww, " %s", msg_split[i]);
    size_print += lsize + 1;
    i++;
  }
  g_string_append_printf (msg_ww, "\n");

  while (msg_split[i] != NULL)
  {
    size_print = 0;
    g_string_append_printf (msg_ww, "%s", rest);
    while (msg_split[i] != NULL)
    {
      guint lsize = strlen (msg_split[i]);
      if (size_print + lsize + 1 > (ncols - rest_size))
        break;
      g_string_append_printf (msg_ww, " %s", msg_split[i]);
      size_print += lsize + 1;
      i++;
    }
    g_string_append_printf (msg_ww, "\n");
  }

  g_strfreev (msg_split);

  return g_string_free (msg_ww, FALSE);
}

/**
 * ncm_message_ww:
 * @msg: FIXME
 * @first: FIXME
 * @rest: FIXME
 * @ncols: FIXME
 *
 * FIXME
 */
void
ncm_message_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols)
{
  gchar *msg_ww = ncm_string_ww (msg, first, rest, ncols);
  g_message ("%s", msg_ww);
  g_free (msg_ww);
}

/**
 * ncm_cfg_msg_sepa:
 *
 * Log a message separator.
 *
 */
void
ncm_cfg_msg_sepa (void)
{
  g_message ("#----------------------------------------------------------------------------------\n");
}

/**
 * ncm_cfg_get_fullpath:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gchar *
ncm_cfg_get_fullpath (const gchar *filename, ...)
{
  gchar *file, *full_filename;

  g_assert (numcosmo_init);
  va_list ap;
  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  full_filename = g_build_filename (numcosmo_path, file, NULL);

  g_free (file);
  return full_filename;
}

/**
 * ncm_cfg_keyfile_to_arg:
 * @kfile: FIXME
 * @group_name: FIXME
 * @entries: FIXME
 * @argv: FIXME
 * @argc: FIXME
 *
 * FIXME
 *
 */
void
ncm_cfg_keyfile_to_arg (GKeyFile *kfile, const gchar *group_name, GOptionEntry *entries, gchar **argv, gint *argc)
{
  if (g_key_file_has_group (kfile, group_name))
  {
    GError *error = NULL;
    gint i;
    for (i = 0; entries[i].long_name != NULL; i++)
    {
      if (g_key_file_has_key (kfile, group_name, entries[i].long_name, &error))
      {
        if (entries[i].arg == G_OPTION_ARG_STRING_ARRAY || entries[i].arg == G_OPTION_ARG_FILENAME_ARRAY)
        {
          guint j;
          gsize length;
          gchar **vals = g_key_file_get_string_list (kfile, group_name, entries[i].long_name, &length, &error);
          if (error != NULL)
            g_error ("ncm_cfg_keyfile_to_arg: Cannot parse key file[%s]", error->message);
          for (j = 0; j < length; j++)
          {
            argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
            argv[argc[0]++] = vals[j];
          }
          g_free (vals);
        }
        else
        {
          gchar *val = g_key_file_get_value (kfile, group_name, entries[i].long_name, &error);
          if (error != NULL)
            g_error ("ncm_cfg_keyfile_to_arg: Cannot parse key file[%s]", error->message);

          if (entries[i].arg == G_OPTION_ARG_NONE)
          {
            if ((g_ascii_strcasecmp (val, "1") == 0) ||
                (g_ascii_strcasecmp (val, "true") == 0))
            {
              argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
            }
            g_free (val);
          }
          else if (strlen (val) > 0)
          {
            argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);
            argv[argc[0]++] = val;
          }
        }
      }
    }
  }
}

/**
 * ncm_cfg_string_to_comment:
 * @str: FIXME
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
gchar *
ncm_cfg_string_to_comment (const gchar *str)
{
  g_assert (str != NULL);
  {
    gchar *desc_ww = ncm_string_ww (str, "  ", "  ", 80);
    gchar *desc = g_strdup_printf ("###############################################################################\n\n%s\n", desc_ww);
    g_free (desc_ww);
    return desc;
  }
}

/**
 * ncm_cfg_entries_to_keyfile:
 * @kfile: FIXME
 * @group_name: FIXME
 * @entries: FIXME
 *
 * FIXME
 *
 */
void
ncm_cfg_entries_to_keyfile (GKeyFile *kfile, const gchar *group_name, GOptionEntry *entries)
{
  GError *error = NULL;
  gint i;
  for (i = 0; entries[i].long_name != NULL; i++)
  {
    gboolean skip_comment = FALSE;
    switch (entries[i].arg)
    {
      case G_OPTION_ARG_NONE:
      {
        gboolean arg_b = ((gboolean *)entries[i].arg_data)[0];
        g_key_file_set_boolean (kfile, group_name, entries[i].long_name, arg_b);
        break;
      }
      case G_OPTION_ARG_STRING:
      case G_OPTION_ARG_FILENAME:
      {
        gchar **arg_s = (gchar **)entries[i].arg_data;
        g_key_file_set_string (kfile, group_name, entries[i].long_name, *arg_s != NULL ? *arg_s : "");
        break;
      }
      case G_OPTION_ARG_STRING_ARRAY:
      case G_OPTION_ARG_FILENAME_ARRAY:
      {
        const gchar ***arg_as = (const gchar ***)entries[i].arg_data;
        gchar ***arg_cas = (gchar ***)entries[i].arg_data;
        if (*arg_cas != NULL)
          g_key_file_set_string_list (kfile, group_name, entries[i].long_name, *arg_as, g_strv_length (*arg_cas));
        else
          g_key_file_set_string_list (kfile, group_name, entries[i].long_name, NULL, 0);
        break;
      }
      case G_OPTION_ARG_INT:
      {
        gint arg_i = ((gint *)entries[i].arg_data)[0];
        g_key_file_set_integer (kfile, group_name, entries[i].long_name, arg_i);
        break;
      }
      case G_OPTION_ARG_INT64:
      {
        gint64 arg_l = ((gint64 *)entries[i].arg_data)[0];
        g_key_file_set_int64 (kfile, group_name, entries[i].long_name, arg_l);
        break;
      }
      case G_OPTION_ARG_DOUBLE:
      {
        gdouble arg_d = ((double *)entries[i].arg_data)[0];
        g_key_file_set_double (kfile, group_name, entries[i].long_name, arg_d);
        break;
      }
      case G_OPTION_ARG_CALLBACK:
      default:
        skip_comment = TRUE;
        //g_error ("ncm_cfg_entries_to_keyfile: cannot convert entry type %d to keyfile", entries[i].arg);
        break;
    }
    if (!skip_comment)
    {
      gchar *desc = ncm_cfg_string_to_comment (entries[i].description);
      if (!g_key_file_set_comment (kfile, group_name, entries[i].long_name, desc, &error))
        g_error ("ncm_cfg_entries_to_keyfile: %s", error->message);
      g_free (desc);
    }
  }
}

/**
 * ncm_cfg_get_enum_by_id_name_nick:
 * @enum_type: FIXME
 * @id_name_nick: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
const GEnumValue *
ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick)
{
  g_assert (id_name_nick != NULL);
  {
    GEnumValue *res = NULL;
    gchar *endptr = NULL;
    gint64 id = g_ascii_strtoll (id_name_nick, &endptr, 10);
    GEnumClass *enum_class = NULL;

    g_assert (G_TYPE_IS_ENUM (enum_type));

    enum_class = g_type_class_ref (enum_type);
    if ((endptr == id_name_nick) || (strlen (endptr) > 0))
    {
      res = g_enum_get_value_by_name (enum_class, id_name_nick);
      if (res == NULL)
        res = g_enum_get_value_by_nick (enum_class, id_name_nick);
    }
    else
      res = g_enum_get_value (enum_class, id);

    g_type_class_unref (enum_class);

    return res;
  }
}

/**
 * ncm_cfg_enum_get_value:
 * @enum_type: FIXME
 * @n: FIXME
 *
 * FIXME
 *
 * Returns: (transfer none)
 */
const GEnumValue *
ncm_cfg_enum_get_value (GType enum_type, guint n)
{
  GEnumClass *enum_class;
  GEnumValue *val;
  g_assert (G_TYPE_IS_ENUM (enum_type));
  enum_class = g_type_class_ref (enum_type);

  val = g_enum_get_value (enum_class, n);

  g_type_class_unref (enum_class);

  return val;
}

/**
 * ncm_cfg_enum_print_all:
 * @enum_type: FIXME
 * @header: FIXME
 *
 * FIXME
 *
 */
void
ncm_cfg_enum_print_all (GType enum_type, const gchar *header)
{
  GEnumClass *enum_class;
  GEnumValue *snia;
  gint i = 0;
  gint name_max_len = 4;
  gint nick_max_len = 4;
  gint pad;

  g_assert (G_TYPE_IS_ENUM (enum_type));

  enum_class = g_type_class_ref (enum_type);

  while ((snia = g_enum_get_value (enum_class, i++)) != NULL)
  {
    name_max_len = GSL_MAX (name_max_len, strlen (snia->value_name));
    nick_max_len = GSL_MAX (nick_max_len, strlen (snia->value_nick));
  }

  printf ("# %s:\n", header);
  pad = 10 + name_max_len + nick_max_len;
  printf ("#");
  while (pad-- != 0)
    printf ("-");
  printf ("#\n");
  printf ("# Id | %-*s | %-*s |\n", name_max_len, "Name", nick_max_len, "Nick");
  i = 0;
  while ((snia = g_enum_get_value (enum_class, i++)) != NULL)
  {
    printf ("# %02d | %-*s | %-*s |\n", snia->value, name_max_len, snia->value_name, nick_max_len, snia->value_nick);
  }
  pad = 10 + name_max_len + nick_max_len;
  printf ("#");
  while (pad-- != 0)
    printf ("-");
  printf ("#\n");

  g_type_class_unref (enum_class);
}

#ifdef NUMCOSMO_HAVE_FFTW3

G_LOCK_DEFINE_STATIC (fftw_saveload_lock);
G_LOCK_DEFINE_STATIC (fftw_plan_lock);

/**
 * ncm_cfg_lock_plan_fftw:
 *
 * FIXME
 *
 */
void
ncm_cfg_lock_plan_fftw (void)
{
  G_LOCK (fftw_plan_lock);
}

/**
 * ncm_cfg_unlock_plan_fftw:
 *
 * FIXME
 *
 */
void
ncm_cfg_unlock_plan_fftw (void)
{
  G_UNLOCK (fftw_plan_lock);
}

/**
 * ncm_cfg_load_fftw_wisdom:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_fftw_wisdom (const gchar *filename, ...)
{
  gchar *file, *file_ext;
  gchar *full_filename;
  va_list ap;
  gboolean ret = FALSE;
  g_assert (numcosmo_init);

  G_LOCK (fftw_saveload_lock);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  g_free (file);
  file = g_strdup ("ncm_cfg_wisdom"); /* overwrite, unifying wisdom */

  file_ext      = g_strdup_printf ("%s.fftw3", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    fftw_import_wisdom_from_filename (full_filename);
    ret = TRUE;
  }

#ifdef HAVE_FFTW3F
  g_free (file_ext);
  g_free (full_filename);

  file_ext      = g_strdup_printf ("%s.fftw3f", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    fftwf_import_wisdom_from_filename (full_filename);
    ret = TRUE;
  }
#endif

  g_free (file);
  g_free (file_ext);
  g_free (full_filename);

  G_UNLOCK (fftw_saveload_lock);

  return ret;
}

/**
 * ncm_cfg_save_fftw_wisdom:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_fftw_wisdom (const gchar *filename, ...)
{
  gchar *file, *file_ext;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  G_LOCK (fftw_saveload_lock);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  g_free (file);
  file = g_strdup ("ncm_cfg_wisdom"); /* overwrite, unifying wisdom */

  file_ext      = g_strdup_printf ("%s.fftw3", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  fftw_export_wisdom_to_filename (full_filename);

#ifdef HAVE_FFTW3F
  g_free (file_ext);
  g_free (full_filename);

  file_ext      = g_strdup_printf ("%s.fftw3f", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  fftwf_export_wisdom_to_filename (full_filename);
#endif

  g_free (file);
  g_free (file_ext);
  g_free (full_filename);

  G_UNLOCK (fftw_saveload_lock);

  return TRUE;
}
#endif /* NUMCOSMO_HAVE_FFTW3 */

/**
 * ncm_cfg_fopen:
 * @filename: FIXME
 * @mode: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
FILE *
ncm_cfg_fopen (const gchar *filename, const gchar *mode, ...)
{
  FILE *F;
  gchar *file;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  va_start (ap, mode);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);
  full_filename = g_build_filename (numcosmo_path, file, NULL);

  F = g_fopen (full_filename, mode);
  if (F == NULL)
  {
    g_error ("ncm_cfg_fopen: cannot open file %s [%s].", full_filename, g_strerror (errno));
  }

  g_free (file);
  g_free (full_filename);
  return F;
}

/**
 * ncm_cfg_vfopen: (skip)
 * @filename: FIXME
 * @mode: FIXME
 * @ap: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
FILE *
ncm_cfg_vfopen (const gchar *filename, const gchar *mode, va_list ap)
{
  FILE *F;
  gchar *file;
  gchar *full_filename;

  g_assert (numcosmo_init);
  file = g_strdup_vprintf (filename, ap);
  full_filename = g_build_filename (numcosmo_path, file, NULL);

  F = g_fopen (full_filename, mode);
  if (F == NULL)
  {
    g_error ("ncm_cfg_fopen: cannot open file %s [%s].", full_filename, g_strerror (errno));
  }

  g_free (file);
  g_free (full_filename);
  return F;
}

/**
 * ncm_cfg_exists:
 * @filename: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_exists (const gchar *filename, ...)
{
  gboolean exists;
  gchar *file;
  gchar *full_filename;
  va_list ap;

  g_assert (numcosmo_init);

  va_start (ap, filename);
  file = g_strdup_vprintf (filename, ap);
  va_end (ap);

  full_filename = g_build_filename (numcosmo_path, file, NULL);
  exists = g_file_test (full_filename, G_FILE_TEST_EXISTS);
  g_free (file);
  g_free (full_filename);
  return exists;
}

/**
 * ncm_cfg_load_spline: (skip)
 * @filename: FIXME
 * @stype: FIXME
 * @s: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_load_spline (const gchar *filename, const gsl_interp_type *stype, NcmSpline **s, ...)
{
  guint64 size;
  NcmVector *xv, *yv;
  va_list ap;

  va_start (ap, s);
  FILE *F = ncm_cfg_vfopen (filename, "r", ap);
  va_end (ap);
  if (F == NULL)
    return FALSE;

  if (fread (&size, sizeof (guint64), 1, F) != 1)
    g_error ("ncm_cfg_save_spline: fwrite error");

  if (*s == NULL)
  {
    xv = ncm_vector_new (size);
    yv = ncm_vector_new (size);
  }
  else
  {
    xv = ncm_spline_get_xv (*s);
    yv = ncm_spline_get_yv (*s);
    g_assert (size == ncm_vector_len(xv));
  }

  if (fread (ncm_vector_ptr (xv, 0), sizeof (gdouble), size, F) != size)
    g_error ("ncm_cfg_save_spline: fwrite error");
  if (fread (ncm_vector_ptr (yv, 0), sizeof (gdouble), size, F) != size)
    g_error ("ncm_cfg_save_spline: fwrite error");

  if (*s == NULL)
  {
    *s = ncm_spline_gsl_new (stype);
    ncm_spline_set (*s, xv, yv, TRUE);
  }
  else
    ncm_spline_prepare (*s);

  ncm_vector_free (xv);
  ncm_vector_free (yv);

  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_save_spline:
 * @filename: FIXME
 * @s: FIXME
 * @...: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gboolean
ncm_cfg_save_spline (const gchar *filename, NcmSpline *s, ...)
{
  guint64 size;
  va_list ap;

  va_start (ap, s);
  FILE *F = ncm_cfg_vfopen (filename, "w", ap);
  va_end (ap);

  if (F == NULL)
    return FALSE;

  size = s->len;

  if (fwrite (&size, sizeof (guint64), 1, F) != 1)
    g_error ("ncm_cfg_save_spline: fwrite error");

  if (fwrite (ncm_vector_ptr (s->xv, 0), sizeof (gdouble), size, F) != size)
    g_error ("ncm_cfg_save_spline: fwrite error");
  if (fwrite (ncm_vector_ptr (s->yv, 0), sizeof (gdouble), size, F) != size)
    g_error ("ncm_cfg_save_spline: fwrite error");

  fclose (F);
  return TRUE;
}

/**
 * ncm_cfg_get_data_filename:
 * @filename: filename to search in the data path.
 * @must_exist: raises an error if @filename is not found.
 *
 * Looks for @filename in the data path and returns
 * the full path if found.
 *
 * Returns: (transfer full): Full path for @filename.
 */
gchar *
ncm_cfg_get_data_filename (const gchar *filename, gboolean must_exist)
{
  const gchar *data_dir = g_getenv (NCM_CFG_DATA_DIR_ENV);
  gchar *full_filename = NULL;

  if (data_dir != NULL)
  {
    full_filename = g_build_filename (data_dir, "data", filename, NULL);

    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
    {
      g_clear_pointer (&full_filename, g_free);
    }
  }

  if (full_filename == NULL)
  {
    full_filename = g_build_filename (PACKAGE_DATA_DIR, "data", filename, NULL);
    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
    {
      g_clear_pointer (&full_filename, g_free);
    }
  }

  if (full_filename == NULL)
    full_filename = g_build_filename (PACKAGE_SOURCE_DIR, "data", filename, NULL);

  if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    if (must_exist)
      g_error ("ncm_cfg_get_data_filename: cannot find `%s'.", filename);
    else
    {
      g_clear_pointer (&full_filename, g_free);
    }
  }

  return full_filename;
}

/**
 * ncm_cfg_command_line:
 * @argv: FIXME
 * @argc: FIXME
 *
 * FIXME
 *
 * Returns: FIXME
 */
gchar *
ncm_cfg_command_line (gchar *argv[], gint argc)
{
  gchar *full_cmd_line;
  gchar *full_cmd_line_ptr;
  guint tsize = (argc - 1) + 1;
  gint argv_size, i;

  for (i = 0; i < argc; i++)
  {
    tsize += strlen (argv[i]);
    if (g_strrstr (argv[i], " ") != NULL)
      tsize += 2;
  }

  full_cmd_line_ptr = full_cmd_line = g_new (gchar, tsize);

  argv_size = strlen(argv[0]);
  memcpy (full_cmd_line_ptr, argv[0], argv_size);
  full_cmd_line_ptr = &full_cmd_line_ptr[argv_size];

  for (i = 1; i < argc; i++)
  {
    gboolean has_space = FALSE;
    full_cmd_line_ptr[0] = ' ';
    full_cmd_line_ptr++;
    argv_size = strlen(argv[i]);
    has_space = (g_strrstr (argv[i], " ") != NULL);
    if (has_space)
      (full_cmd_line_ptr++)[0] = '\'';
    memcpy (full_cmd_line_ptr, argv[i], argv_size);
    full_cmd_line_ptr = &full_cmd_line_ptr[argv_size];
    if (has_space)
      (full_cmd_line_ptr++)[0] = '\'';
  }

  full_cmd_line_ptr[0] = '\0';

  return full_cmd_line;
}

/**
 * ncm_cfg_variant_to_array: (skip)
 * @var: a variant of array type.
 * @esize: element size.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
GArray *
ncm_cfg_variant_to_array (GVariant *var, gsize esize)
{
  g_assert (g_variant_is_of_type (var, G_VARIANT_TYPE_ARRAY));
  {
    GArray *a = g_array_new (FALSE, FALSE, esize);
    ncm_cfg_array_set_variant (a, var);
    return a;
  }
}

/**
 * ncm_cfg_array_set_variant: (skip)
 * @a: a #GArray.
 * @var: a variant of array type.
 *
 * FIXME
 *
 */
void
ncm_cfg_array_set_variant (GArray *a, GVariant *var)
{
  gsize esize = g_array_get_element_size (a);
  gsize n_elements = 0;
  gconstpointer data = g_variant_get_fixed_array (var, &n_elements, esize);
  g_array_set_size (a, n_elements);
  memcpy (a->data, data, n_elements * esize);
}

/**
 * ncm_cfg_array_to_variant: (skip)
 * @a: a #GArray.
 * @etype: element type.
 *
 * FIXME
 *
 * Returns: (transfer full): FIXME
 */
GVariant *
ncm_cfg_array_to_variant (GArray *a, const GVariantType *etype)
{
  gconstpointer data = a->data;
  gsize esize = g_array_get_element_size (a);
  GVariantType *atype = g_variant_type_new_array (etype);
  GVariant *vvar = g_variant_new_from_data (atype,
                                            data,
                                            esize * a->len,
                                            TRUE,
                                            (GDestroyNotify) &g_array_unref,
                                            g_array_ref (a));
  g_variant_type_free (atype);
  return g_variant_ref_sink (vvar);
}

gdouble fftw_default_timeout = 60.0;

#ifdef NUMCOSMO_HAVE_FFTW3
guint fftw_default_flags = FFTW_MEASURE; /* FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE */

/**
 * ncm_cfg_set_fftw_default_flag:
 * @flag: a FFTW library flag
 * @timeout: planner time out in seconds
 *
 * Sets the default FFTW flag (FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE)
 * to be used when building plans. The variable @timeout sets the maximum time spended on
 * planners. 
 *
 */
void
ncm_cfg_set_fftw_default_flag (guint flag, const gdouble timeout)
{
  fftw_default_flags = flag;
  fftw_set_timelimit (10.0);
#ifdef HAVE_FFTW3F
  fftwf_set_timelimit (10.0);
#endif /* HAVE_FFTW3F */
}
#else
guint fftw_default_flags = 0;
#endif /* NUMCOSMO_HAVE_FFTW3 */
