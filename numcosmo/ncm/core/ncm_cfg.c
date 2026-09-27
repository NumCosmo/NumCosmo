/***************************************************************************
 *            ncm_cfg.c
 *
 *  Wed Aug 13 20:59:22 2008
 *  Copyright  2008  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <vitenti@uel.br>
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
 * NcmCfg:
 *
 * Library initialization and configuration.
 *
 * Initialization, log output, thread counts of the linear algebra libraries, FFTW
 * planning and wisdom, data-file lookup, and helpers for #GOptionEntry, #GKeyFile and
 * enumeration types.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif /* HAVE_CONFIG_H */
#include "build_cfg.h"

#include "ncm/core/ncm_cfg.h"
#include "ncm/core/ncm_rng.h"
#include "ncm/core/ncm_memory_pool.h"
#include "ncm/algebra/ncm_complex.h"
#include "ncm/mpi/ncm_mpi_job.h"
#include "ncm/mpi/ncm_mpi_job_test.h"
#include "ncm/mpi/ncm_mpi_job_fit.h"
#include "ncm/mpi/ncm_mpi_job_mcmc.h"
#include "ncm/mpi/ncm_mpi_job_feval.h"
#include "ncm/algebra/ncm_vector.h"
#include "ncm/spline/ncm_spline_bspline.h"
#include "ncm/spline/ncm_spline_gsl.h"
#include "ncm/spline/ncm_spline_cubic.h"
#include "ncm/spline/ncm_spline_cubic_notaknot.h"
#include "ncm/spline/ncm_spline_cubic_d2.h"
#include "ncm/spline/ncm_spline_vec.h"
#include "ncm/stats/ncm_function_sample_set.h"
#include "ncm/spline/ncm_spline2d_bicubic.h"
#include "ncm/spline/ncm_spline2d_gsl.h"
#include "ncm/spline/ncm_spline2d_spline.h"
#include "ncm/integration/ncm_integral1d.h"
#include "ncm/integration/ncm_integral_nd.h"
#include "ncm/core/ncm_pln1d.h"
#include "ncm/powspec/ncm_powspec_corr3d.h"
#include "ncm/powspec/ncm_powspec_filter.h"
#include "ncm/powspec/ncm_powspec_spline2d.h"
#include "ncm/powspec/tests/ncm_powspec_analytic.h"
#include "ncm/powspec/ncm_powspec.h"
#include "ncm/model/ncm_model.h"
#include "ncm/model/ncm_model_ctrl.h"
#include "ncm/model/ncm_model_builder.h"
#include "ncm/model/ncm_model_mvnd.h"
#include "ncm/model/ncm_model_rosenbrock.h"
#include "ncm/model/ncm_model_funnel.h"
#include "ncm/model/ncm_reparam_linear.h"
#include "ncm/data/ncm_data.h"
#include "ncm/data/ncm_data_gauss_cov_mvnd.h"
#include "ncm/data/ncm_data_rosenbrock.h"
#include "ncm/data/ncm_data_funnel.h"
#include "ncm/data/ncm_data_gaussmix2d.h"
#include "ncm/stats/ncm_stats_acorr.h"
#include "ncm/stats/ncm_stats_vec.h"
#include "ncm/fit/ncm_fit_esmcmc_walker_stretch.h"
#include "ncm/data/ncm_data.h"
#include "ncm/stats/ncm_stats_dist1d_epdf.h"
#include "ncm/stats/ncm_stats_dist1d_spline.h"
#include "ncm/data/ncm_dataset.h"
#include "ncm/fit/ncm_fit.h"
#include "ncm/fit/ncm_fit_gsl_ls.h"
#include "ncm/fit/ncm_fit_gsl_mm.h"
#include "ncm/fit/ncm_fit_gsl_mms.h"
#include "ncm/fit/ncm_fit_nlopt.h"
#include "ncm/fit/ncm_prior_gauss_param.h"
#include "ncm/fit/ncm_prior_gauss_func.h"
#include "ncm/fit/ncm_prior_flat_param.h"
#include "ncm/fit/ncm_prior_flat_func.h"
#include "ncm/specfunc/ncm_sbessel_integrator.h"
#include "ncm/specfunc/ncm_sbessel_integrator_gl.h"
#include "ncm/specfunc/ncm_sbessel_integrator_levin.h"
#include "ncm/specfunc/ncm_sbessel_ode_solver.h"
#include "ncm/fftlog/ncm_fftlog_sbessel_j.h"
#include "ncm/algebra/ncm_spectral.h"
#include "nc/background/nc_hicosmo.h"
#include "nc/bbn/nc_bbn.h"
#include "nc/bbn/nc_bbn_parametrized.h"
#include "nc/bbn/nc_bbn_parthenope.h"
#include "nc/cmb/nc_cbe_precision.h"
#include "nc/background/nc_hicosmo_qconst.h"
#include "nc/background/nc_hicosmo_qlinear.h"
#include "nc/background/nc_hicosmo_qspline.h"
#include "nc/background/nc_hicosmo_qrbf.h"
#include "nc/background/nc_hicosmo_lcdm.h"
#include "nc/background/nc_hicosmo_de_xcdm.h"
#include "nc/background/nc_hicosmo_de_wspline.h"
#include "nc/background/nc_hicosmo_de_cpl.h"
#include "nc/background/nc_hicosmo_de_jbp.h"
#include "nc/background/nc_hicosmo_qgrw.h"
#include "nc/background/nc_hicosmo_qgw.h"
#include "nc/background/nc_hicosmo_Vexp.h"
#include "nc/background/nc_hicosmo_de_reparam_ok.h"
#include "nc/background/nc_hicosmo_de_reparam_cmb.h"
#include "nc/primordial/nc_hiprim_atan.h"
#include "nc/primordial/nc_hiprim_bpl.h"
#include "nc/primordial/nc_hiprim_expc.h"
#include "nc/primordial/nc_hiprim_power_law.h"
#include "nc/primordial/nc_hiprim_sbpl.h"
#include "nc/primordial/nc_hiprim_two_fluids.h"
#include "nc/powspec/nc_window_tophat.h"
#include "nc/powspec/nc_window_gaussian.h"
#include "nc/powspec/nc_growth_func.h"
#include "nc/powspec/nc_transfer_func.h"
#include "nc/powspec/nc_transfer_func_bbks.h"
#include "nc/powspec/nc_transfer_func_eh.h"
#include "nc/powspec/nc_transfer_func_eh_no_baryon.h"
#include "nc/powspec/nc_transfer_func_camb.h"
#include "nc/lss/halo/nc_halo_position.h"
#include "nc/lss/halo/nc_halo_density_profile.h"
#include "nc/lss/halo/nc_halo_density_profile_nfw.h"
#include "nc/lss/halo/nc_halo_density_profile_einasto.h"
#include "nc/lss/halo/nc_halo_density_profile_dk14.h"
#include "nc/lss/halo/nc_halo_density_profile_hernquist.h"
#include "nc/lss/halo/nc_halo_mass_summary.h"
#include "nc/lss/halo/nc_halo_cm_param.h"
#include "nc/lss/halo/nc_halo_cm_duffy08.h"
#include "nc/lss/halo/nc_halo_cm_klypin11.h"
#include "nc/lss/halo/nc_halo_cm_prada12.h"
#include "nc/lss/halo/nc_halo_cm_bhattacharya13.h"
#include "nc/lss/halo/nc_halo_cm_dutton14.h"
#include "nc/lss/halo/nc_halo_cm_diemer15.h"
#include "nc/lss/halo/nc_multiplicity_func.h"
#include "nc/lss/halo/nc_multiplicity_func_st.h"
#include "nc/lss/halo/nc_multiplicity_func_ps.h"
#include "nc/lss/halo/nc_multiplicity_func_jenkins.h"
#include "nc/lss/halo/nc_multiplicity_func_warren.h"
#include "nc/lss/halo/nc_multiplicity_func_tinker.h"
#include "nc/lss/halo/nc_multiplicity_func_tinker_mean_normalized.h"
#include "nc/lss/halo/nc_multiplicity_func_crocce.h"
#include "nc/lss/halo/nc_multiplicity_func_bocquet.h"
#include "nc/lss/halo/nc_multiplicity_func_castro.h"
#include "nc/lss/halo/nc_multiplicity_func_despali.h"
#include "nc/lss/halo/nc_multiplicity_func_watson.h"
#include "nc/lss/halo/nc_multiplicity_func_bhattacharya.h"
#include "nc/lss/halo/nc_halo_mass_function.h"
#include "nc/lss/galaxy/nc_galaxy_acf.h"
#include "nc/lss/cluster/nc_cluster_mass.h"
#include "nc/lss/cluster/nc_cluster_mass_nodist.h"
#include "nc/lss/cluster/nc_cluster_mass_lnnormal.h"
#include "nc/lss/cluster/nc_cluster_mass_vanderlinde.h"
#include "nc/lss/cluster/nc_cluster_mass_benson.h"
#include "nc/lss/cluster/nc_cluster_mass_benson_xray.h"
#include "nc/lss/cluster/nc_cluster_mass_plcl.h"
#include "nc/lss/cluster/nc_cluster_mass_ascaso.h"
#include "nc/lss/cluster/nc_cluster_mass_selection.h"
#include "nc/lss/cluster/nc_cluster_redshift.h"
#include "nc/lss/cluster/nc_cluster_redshift_nodist.h"
#include "nc/lss/cluster/nc_cluster_photoz_gauss_global.h"
#include "nc/lss/cluster/nc_cluster_photoz_gauss.h"
#include "nc/lss/halo/nc_halo_bias_castro.h"
#include "nc/lss/halo/nc_halo_bias_despali.h"
#include "nc/lss/halo/nc_halo_bias_ps.h"
#include "nc/lss/halo/nc_halo_bias_st_ellip.h"
#include "nc/lss/halo/nc_halo_bias_st_spher.h"
#include "nc/lss/halo/nc_halo_bias_tinker.h"
#include "nc/lss/cluster/nc_cluster_abundance.h"
#include "nc/lss/cluster/nc_cluster_pseudo_counts.h"
#include "nc/lss/cluster/nc_cor_cluster_cmb_lens_limber.h"
#include "nc/lss/wl/nc_wl_surface_mass_density.h"
#include "nc/lss/wl/nc_reduced_shear_cluster_mass.h"
#include "nc/lss/wl/nc_reduced_shear_calib.h"
#include "nc/lss/wl/nc_reduced_shear_calib_wtg.h"
#include "nc/lss/wl/nc_wl_ellipticity_series.h"
#include "nc/lss/galaxy/nc_galaxy_wl_obs.h"
#include "nc/lss/galaxy/nc_galaxy_position_factor.h"
#include "nc/lss/galaxy/nc_galaxy_position_factor_flat.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_factor.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_factor_composed.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_factor_spline.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_obs.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_obs_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_obs_sel.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_obs_sel_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_pop.h"
#include "nc/lss/galaxy/nc_galaxy_redshift_pop_lsst_srd.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_var_add.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_laplace.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_quad.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_fixed_quad.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_series_lensed.h"
#include "nc/lss/galaxy/nc_galaxy_shape_factor_cgf.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_gauss_local.h"
#include "nc/lss/galaxy/nc_galaxy_shape_pop_beta.h"
#include "nc/background/nc_distance.h"
#include "nc/recomb/nc_recomb.h"
#include "nc/recomb/nc_recomb_cbe.h"
#include "nc/recomb/nc_recomb_seager.h"
#include "nc/reion/nc_hireion.h"
#include "nc/reion/nc_hireion_camb.h"
#include "nc/powspec/nc_powspec_ml_cbe.h"
#include "nc/powspec/nc_powspec_ml_spline.h"
#include "nc/powspec/nc_powspec_ml_transfer.h"
#include "nc/powspec/nc_powspec_ml.h"
#include "nc/powspec/nc_powspec_mnl_halofit.h"
#include "nc/powspec/nc_powspec_mnl.h"
#include "nc/supernova/nc_snia_dist_cov.h"
#include "nc/cmb/nc_planck_fi.h"
#include "nc/cmb/nc_planck_fi_cor_tt.h"
#include "nc/cmb/nc_planck_fi_cor_ttteee.h"
#include "nc/perturbations/nc_hipert_boltzmann_cbe.h"
#include "nc/data/nc_data_bao_a.h"
#include "nc/data/nc_data_bao_dv.h"
#include "nc/data/nc_data_bao_dvdv.h"
#include "nc/data/nc_data_bao_rdv.h"
#include "nc/data/nc_data_bao_empirical_fit.h"
#include "nc/data/nc_data_bao_empirical_fit_2d.h"
#include "nc/data/nc_data_bao_dhr_dar.h"
#include "nc/data/nc_data_bao_dtr_dhr.h"
#include "nc/data/nc_data_bao_dmr_hr.h"
#include "nc/data/nc_data_bao_dvr_dtdh.h"
#include "nc/data/nc_data_dist_mu.h"
#include "nc/data/nc_data_cluster_pseudo_counts.h"
#include "nc/data/nc_data_cluster_ncount.h"
#include "nc/data/nc_data_cluster_ncounts_gauss.h"
#include "nc/data/nc_data_cluster_wl_factor.h"
#include "nc/data/nc_data_cluster_mass_rich.h"
#include "nc/data/nc_data_cluster_mass_rich_count.h"
#include "nc/data/nc_data_cmb_shift_param.h"
#include "nc/data/nc_data_cmb_dist_priors.h"
#include "nc/data/nc_data_hubble.h"
#include "nc/data/nc_data_snia_cov.h"
#include "nc/data/nc_data_xcor.h"
#include "nc/data/nc_data_planck_lkl.h"
#include "nc/data/nc_data_planck_commander.h"
#include "nc/data/nc_data_planck_lensing.h"
#include "nc/data/nc_data_planck_plik_lite.h"
#include "nc/data/nc_data_planck_simall.h"
#include "nc/data/nc_data_planck_smica.h"
#include "nc/xcor/nc_xcor.h"
#include "nc/xcor/nc_xcor_AB.h"
#include "nc/xcor/nc_xcor_solver.h"
#include "nc/xcor/nc_xcor_ssc_sij.h"
#include "nc/xcor/nc_xcor_kernel.h"
#include "nc/xcor/nc_xcor_kernel_component.h"
#include "nc/xcor/nc_xcor_kernel_gal.h"
#include "nc/xcor/nc_xcor_kernel_cluster.h"
#include "nc/xcor/nc_xcor_kernel_cluster_tophat.h"
#include "nc/xcor/nc_xcor_kernel_cmb_isw.h"
#include "nc/xcor/nc_xcor_kernel_CMB_lensing.h"
#include "nc/xcor/nc_xcor_kernel_weak_lensing.h"
#include "nc/xcor/nc_xcor_kernel_tSZ.h"
#include "nc/xcor/nc_xcor_kernel_radial_kdep.h"
#include "nc/xcor/nc_xcor_kernel_radial.h"
#include "nc/xcor/nc_xcor_kernel_table.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_gauss.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_tophat.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_multi.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_student_t.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_power_exp.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_tophat_smooth.h"
#include "nc/xcor/tests/nc_xcor_kernel_analytic_lensing.h"

#ifndef NUMCOSMO_GIR_SCAN
#include <stdlib.h>
#include <gio/gio.h>
#include <fftw3.h>
#include <cuba.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif /* HAVE_MPI */

#ifdef HAVE_BLIS
#include <blis/blis.h>
#endif /* HAVE_BLIS */

#ifndef G_VALUE_INIT
#define G_VALUE_INIT {0}
#endif

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif /* HAVE_EXECINFO_H */
#endif /* NUMCOSMO_GIR_SCAN */

/* *INDENT-OFF* */
G_DEFINE_QUARK (ncm-cfg-error, ncm_cfg_error)
/* *INDENT-ON* */

static gchar *numcosmo_path         = NULL;
static gboolean numcosmo_init       = FALSE;
static FILE *_log_stream            = NULL;
static FILE *_log_stream_err        = NULL;
static guint _log_msg_id            = 0;
static guint _log_err_id            = 0;
static gboolean _enable_msg         = TRUE;
static gboolean _enable_msg_flush   = TRUE;
static gsl_error_handler_t *gsl_err = NULL;

# if (defined (__GNUC__)                                            \
  && ((__GNUC__ == 11 && __GNUC_MINOR__ >= 1) || (__GNUC__ >= 12))) \
  || (defined (__clang__) && (__clang_major__ >= 12))
extern void __gcov_dump (void);
extern void __gcov_reset (void);

#  define __gcov_flush()                   \
        do {                               \
          __gcov_dump (); __gcov_reset (); \
        } while (0)
# else
extern void __gcov_flush (void);

# endif

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
    gsize size    = backtrace (tarray, 30);
    gchar **trace = backtrace_symbols (tarray, size);
    gsize i;

    /* print out all the frames to stderr */
    for (i = 0; i < size; i++)
    {
      fprintf (_log_stream_err, "# (%s): %s-BACKTRACE:[%02zd] %s\n", pname, log_domain, i, trace[i]);
    }

    g_free (trace);
  }
#endif
  fflush (_log_stream_err);

#ifdef USE_GCOV
  __gcov_flush ();
#endif

  abort ();
}

#ifdef HAVE_OPENBLAS_SET_NUM_THREADS
void goto_set_num_threads (gint);
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
void _nc_hicosmo_de_wspline_register_functions (void);
void _nc_hicosmo_qspline_register_functions (void);
void _nc_galaxy_shape_pop_beta_register_functions (void);

#ifdef HAVE_MPI
static void _ncm_cfg_mpi_main_loop (void);
static gboolean _ncm_cfg_mpi_launched (void);

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
  fftw_forget_wisdom ();
}

/**
 * ncm_cfg_init:
 *
 * Initializes the library; it must be called before any other NumCosmo function, and
 * later calls return immediately. It
 *
 * - sets the default FFTW planner flag and time limit from `NCM_FFTW_PLANNER` and
 *   `NCM_FFTW_PLANNER_TIMELIMIT`, see ncm_cfg_set_fftw_default_from_env_str(); the time
 *   limit is 10 s when `NCM_FFTW_PLANNER_TIMELIMIT` is not set;
 * - creates the directory `~/.numcosmo`, see ncm_cfg_get_fullpath();
 * - sets the Cuba library core counts to zero;
 * - turns the GSL error handler off, see ncm_cfg_enable_gsl_err_handler();
 * - installs the NumCosmo log handlers;
 * - registers the library objects and functions;
 * - under an MPI launcher, initializes MPI; every rank except the master then runs the
 *   worker loop and never returns.
 *
 * MPI is detected from the launcher's environment (Open MPI, MPICH/Hydra, Slurm);
 * `NUMCOSMO_MPI_INIT` set to 1 or 0 forces or skips it.
 *
 * Same as ncm_cfg_init_full() without command-line arguments.
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
  {
    return NULL;
  }
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
 * @argc: number of arguments
 * @argv: (array length=argc): the arguments
 *
 * Same as ncm_cfg_init_full_ptr(), for bindings: the arguments are copied and the
 * possibly modified copy is returned.
 *
 * Returns: (transfer full) (array zero-terminated=1): the arguments after MPI
 * initialization.
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

#ifdef HAVE_MPI

/*
 * Whether a parallel launcher started this process. MPI_Init costs about a second of
 * device probing, so it is called only under a launcher. The launcher's environment is
 * used, instead of initializing and backing out, because every rank except the master
 * enters the worker loop during initialization and never returns.
 */
static gboolean
_ncm_cfg_mpi_launched (void)
{
  const gchar *vars[] = {
    "OMPI_COMM_WORLD_SIZE", /* Open MPI      */
    "PMIX_RANK",            /* Open MPI 4+   */
    "PMI_SIZE",             /* MPICH, Hydra  */
    "PMI_RANK",             /* MPICH, Hydra  */
    "SLURM_PROCID",         /* srun          */
  };
  const gchar *force = g_getenv ("NUMCOSMO_MPI_INIT");
  guint i;

  if (force != NULL)
    return g_strcmp0 (force, "0") != 0;

  for (i = 0; i < G_N_ELEMENTS (vars); i++)
  {
    if (g_getenv (vars[i]) != NULL)
      return TRUE;
  }

  return FALSE;
}

#endif /* HAVE_MPI */

/**
 * ncm_cfg_init_full_ptr:
 * @argc: a pointer to the number of arguments
 * @argv: (array length=argc): a pointer to the arguments
 *
 * Initializes the library; it must be called before any other NumCosmo function, and
 * later calls return immediately. It
 *
 * - sets the default FFTW planner flag and time limit from `NCM_FFTW_PLANNER` and
 *   `NCM_FFTW_PLANNER_TIMELIMIT`, see ncm_cfg_set_fftw_default_from_env_str(); the time
 *   limit is 10 s when `NCM_FFTW_PLANNER_TIMELIMIT` is not set;
 * - creates the directory `~/.numcosmo`, see ncm_cfg_get_fullpath();
 * - sets the Cuba library core counts to zero;
 * - turns the GSL error handler off, see ncm_cfg_enable_gsl_err_handler();
 * - installs the NumCosmo log handlers;
 * - registers the library objects and functions;
 * - under an MPI launcher, initializes MPI; every rank except the master then runs the
 *   worker loop and never returns.
 *
 * MPI is detected from the launcher's environment (Open MPI, MPICH/Hydra, Slurm);
 * `NUMCOSMO_MPI_INIT` set to 1 or 0 forces or skips it.
 *
 * @argc and @argv, as received by main(), are passed to MPI_Init().
 */
void
ncm_cfg_init_full_ptr (gint *argc, gchar ***argv)
{
  const gchar *home;

  if (numcosmo_init)
    return;

  ncm_cfg_set_fftw_default_from_env_str (NUMCOSMO_FFTW_PLAN, 10.0, NULL);

  if (sizeof (NcmComplex) != sizeof (fftw_complex))
    g_warning ("NcmComplex is not binary compatible with complex double, expect problems with it!");

  home          = g_get_home_dir ();
  numcosmo_path = g_build_filename (home, ".numcosmo", NULL);

  if (!g_file_test (numcosmo_path, G_FILE_TEST_EXISTS))
    g_mkdir_with_parents (numcosmo_path, 0755);

  g_setenv ("CUBACORES", "0", TRUE);
  g_setenv ("CUBACORESMAX", "0", TRUE);
  g_setenv ("CUBAACCEL", "0", TRUE);
  g_setenv ("CUBAACCELMAX", "0", TRUE);
#ifdef HAVE_LIBCUBA_4_0
  cubaaccel (0, 0);
  cubacores (0, 0);
#endif

  gsl_err = gsl_set_error_handler_off ();

  _log_stream     = stdout;
  _log_stream_err = stderr;

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message, NULL);
  _log_err_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, _ncm_cfg_log_error, NULL);

  ncm_cfg_register_objects ();
  ncm_cfg_register_functions ();

  numcosmo_init = TRUE;

  _mpi_ctrl.initialized    = 0;
  _mpi_ctrl.size           = 1;
  _mpi_ctrl.rank           = 0;
  _mpi_ctrl.nslaves        = 0;
  _mpi_ctrl.working_slaves = 0;

  atexit (_ncm_cfg_exit);

#ifdef HAVE_MPI

  /* Only under a parallel launcher; see _ncm_cfg_mpi_launched(). */
  if (_ncm_cfg_mpi_launched ())
  {
    MPI_Initialized (&_mpi_ctrl.initialized);

    if (!_mpi_ctrl.initialized)
    {
      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] MPI not initialized, calling MPI_Init.\n", _mpi_ctrl.size, _mpi_ctrl.rank);
      MPI_Init (argc, argv);
      MPI_Initialized (&_mpi_ctrl.initialized);
    }
    else
    {
      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] MPI was already initialized!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
    }

    {
      gchar mpi_hostname[MPI_MAX_PROCESSOR_NAME];
      gint len = 0;

      MPI_Comm_size (MPI_COMM_WORLD, &_mpi_ctrl.size);
      MPI_Comm_rank (MPI_COMM_WORLD, &_mpi_ctrl.rank);
      MPI_Get_processor_name (mpi_hostname, &len);

      NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] We have %d MPI process!! My rank is %d and I'm running on `%s'.\n", _mpi_ctrl.size, _mpi_ctrl.rank, _mpi_ctrl.size, _mpi_ctrl.rank, mpi_hostname);

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
  }
  else
  {
    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] No MPI launcher detected, skipping MPI_Init.\n",
                             _mpi_ctrl.size, _mpi_ctrl.rank);
  }

#endif /* HAVE_MPI */

  return;
}

/**
 * ncm_cfg_register_objects:
 *
 * Registers the library types with ncm_cfg_register_obj(). Called by ncm_cfg_init().
 */
void
ncm_cfg_register_objects (void)
{
  ncm_cfg_register_obj (NCM_TYPE_RNG);

  ncm_cfg_register_obj (NCM_TYPE_VECTOR);
  ncm_cfg_register_obj (NCM_TYPE_MATRIX);

  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_TEST);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_FIT);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_MCMC);
  ncm_cfg_register_obj (NCM_TYPE_MPI_JOB_FEVAL);

  ncm_cfg_register_obj (NCM_TYPE_SPLINE);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC_NOTAKNOT);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_CUBIC_D2);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_BSPLINE);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_GSL);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE_VEC);
  ncm_cfg_register_obj (NCM_TYPE_FUNCTION_SAMPLE_SET);

  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_BICUBIC);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_GSL);
  ncm_cfg_register_obj (NCM_TYPE_SPLINE2D_SPLINE);

  ncm_cfg_register_obj (NCM_TYPE_INTEGRAL1D);
  ncm_cfg_register_obj (NCM_TYPE_INTEGRAL1D_PTR);
  ncm_cfg_register_obj (NCM_TYPE_INTEGRAL_ND);

  ncm_cfg_register_obj (NCM_TYPE_PLN1D);
  ncm_cfg_register_obj (NCM_TYPE_SPECTRAL);

  ncm_cfg_register_obj (NCM_TYPE_SBESSEL_INTEGRATOR_GL);
  ncm_cfg_register_obj (NCM_TYPE_SBESSEL_INTEGRATOR_LEVIN);
  ncm_cfg_register_obj (NCM_TYPE_SBESSEL_ODE_SOLVER);

  ncm_cfg_register_obj (NCM_TYPE_POWSPEC);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_SPLINE2D);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_FILTER);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_CORR3D);
  ncm_cfg_register_obj (NCM_TYPE_POWSPEC_ANALYTIC);

  ncm_cfg_register_obj (NCM_TYPE_MODEL);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_CTRL);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_BUILDER);

  ncm_cfg_register_obj (NCM_TYPE_MODEL_MVND);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_ROSENBROCK);
  ncm_cfg_register_obj (NCM_TYPE_MODEL_FUNNEL);

  ncm_cfg_register_obj (NCM_TYPE_REPARAM);
  ncm_cfg_register_obj (NCM_TYPE_REPARAM_LINEAR);

  ncm_cfg_register_obj (NCM_TYPE_BOOTSTRAP);
  ncm_cfg_register_obj (NCM_TYPE_STATS_ACORR);
  ncm_cfg_register_obj (NCM_TYPE_STATS_VEC);

  ncm_cfg_register_obj (NCM_TYPE_FIT_ESMCMC_WALKER_STRETCH);

  ncm_cfg_register_obj (NCM_TYPE_DATA);
  ncm_cfg_register_obj (NCM_TYPE_DATASET);

  ncm_cfg_register_obj (NCM_TYPE_DATA_GAUSS_COV_MVND);
  ncm_cfg_register_obj (NCM_TYPE_DATA_ROSENBROCK);
  ncm_cfg_register_obj (NCM_TYPE_DATA_FUNNEL);
  ncm_cfg_register_obj (NCM_TYPE_DATA_GAUSSMIX2D);

  ncm_cfg_register_obj (NCM_TYPE_FIT);

  ncm_cfg_register_obj (NCM_TYPE_FIT_GSL_LS);
  ncm_cfg_register_obj (NCM_TYPE_FIT_GSL_MM);
  ncm_cfg_register_obj (NCM_TYPE_FIT_GSL_MMS);

  ncm_cfg_register_obj (NCM_TYPE_FIT_NLOPT);

  ncm_cfg_register_obj (NCM_TYPE_PRIOR_GAUSS_PARAM);
  ncm_cfg_register_obj (NCM_TYPE_PRIOR_GAUSS_FUNC);

  ncm_cfg_register_obj (NCM_TYPE_PRIOR_FLAT_PARAM);
  ncm_cfg_register_obj (NCM_TYPE_PRIOR_FLAT_FUNC);


  ncm_cfg_register_obj (NCM_TYPE_FFTLOG_SBESSEL_J);

  ncm_cfg_register_obj (NCM_TYPE_DATA);

  ncm_cfg_register_obj (NCM_TYPE_STATS_DIST1D_EPDF);
  ncm_cfg_register_obj (NCM_TYPE_STATS_DIST1D_SPLINE);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QCONST);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QLINEAR);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QSPLINE);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QSPLINE_CONT_PRIOR);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QRBF);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_LCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_XCDM);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_WSPLINE);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_CPL);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_JBP);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QGRW);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_QGW);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_VEXP);

  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_REPARAM_OK);
  ncm_cfg_register_obj (NC_TYPE_HICOSMO_DE_REPARAM_CMB);



  ncm_cfg_register_obj (NC_TYPE_HIPRIM_ATAN);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_BPL);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_EXPC);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_POWER_LAW);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_SBPL);
  ncm_cfg_register_obj (NC_TYPE_HIPRIM_TWO_FLUIDS);

  ncm_cfg_register_obj (NC_TYPE_CBE_PRECISION);

  ncm_cfg_register_obj (NC_TYPE_WINDOW);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_TOPHAT);
  ncm_cfg_register_obj (NC_TYPE_WINDOW_GAUSSIAN);

  ncm_cfg_register_obj (NC_TYPE_GROWTH_FUNC);

  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_BBKS);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_EH);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_EH_NO_BARYON);
  ncm_cfg_register_obj (NC_TYPE_TRANSFER_FUNC_CAMB);

  ncm_cfg_register_obj (NC_TYPE_HALO_POSITION);

  ncm_cfg_register_obj (NC_TYPE_HALO_DENSITY_PROFILE);
  ncm_cfg_register_obj (NC_TYPE_HALO_DENSITY_PROFILE_NFW);
  ncm_cfg_register_obj (NC_TYPE_HALO_DENSITY_PROFILE_EINASTO);
  ncm_cfg_register_obj (NC_TYPE_HALO_DENSITY_PROFILE_DK14);
  ncm_cfg_register_obj (NC_TYPE_HALO_DENSITY_PROFILE_HERNQUIST);

  ncm_cfg_register_obj (NC_TYPE_HALO_MASS_SUMMARY);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_PARAM);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_DUFFY08);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_KLYPIN11);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_PRADA12);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_BHATTACHARYA13);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_DUTTON14);
  ncm_cfg_register_obj (NC_TYPE_HALO_CM_DIEMER15);


  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_PS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_ST);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_JENKINS);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_WARREN);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_TINKER_MEAN_NORMALIZED);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_CROCCE);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_BOCQUET);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_CASTRO);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_DESPALI);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_WATSON);
  ncm_cfg_register_obj (NC_TYPE_MULTIPLICITY_FUNC_BHATTACHARYA);

  ncm_cfg_register_obj (NC_TYPE_HALO_MASS_FUNCTION);

  ncm_cfg_register_obj (NC_TYPE_GALAXY_ACF);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_LNNORMAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_VANDERLINDE);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_BENSON_XRAY);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_PLCL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_ASCASO);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_MASS_SELECTION);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_REDSHIFT_NODIST);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS_GLOBAL);
  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PHOTOZ_GAUSS);

  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_CASTRO);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_DESPALI);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_PS);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_ST_ELLIP);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_ST_SPHER);
  ncm_cfg_register_obj (NC_TYPE_HALO_BIAS_TINKER);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_ABUNDANCE);

  ncm_cfg_register_obj (NC_TYPE_CLUSTER_PSEUDO_COUNTS);

  ncm_cfg_register_obj (NC_TYPE_COR_CLUSTER_CMB_LENS_LIMBER);

  ncm_cfg_register_obj (NC_TYPE_WL_SURFACE_MASS_DENSITY);

  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CLUSTER_MASS);

  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CALIB);
  ncm_cfg_register_obj (NC_TYPE_REDUCED_SHEAR_CALIB_WTG);
  ncm_cfg_register_obj (NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE);
  ncm_cfg_register_obj (NC_TYPE_WL_ELLIPTICITY_SERIES_TRACE_DET);

  ncm_cfg_register_obj (NC_TYPE_GALAXY_WL_OBS);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_POSITION_FACTOR);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_POSITION_FACTOR_FLAT);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_FACTOR);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_FACTOR_COMPOSED);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_FACTOR_SPLINE);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_OBS);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_OBS_GAUSS);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_OBS_SEL);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_OBS_SEL_GAUSS);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_POP);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_REDSHIFT_POP_LSST_SRD);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_VAR_ADD);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_LAPLACE);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_QUAD);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_FIXED_QUAD);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_SERIES_LENSED);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_FACTOR_CGF);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_POP);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_POP_GAUSS);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_POP_GAUSS_LOCAL);
  ncm_cfg_register_obj (NC_TYPE_GALAXY_SHAPE_POP_BETA);

  ncm_cfg_register_obj (NC_TYPE_DISTANCE);

  ncm_cfg_register_obj (NC_TYPE_RECOMB);
  ncm_cfg_register_obj (NC_TYPE_RECOMB_CBE);
  ncm_cfg_register_obj (NC_TYPE_RECOMB_SEAGER);

  ncm_cfg_register_obj (NC_TYPE_BBN);
  ncm_cfg_register_obj (NC_TYPE_BBN_PARAMETRIZED);
  ncm_cfg_register_obj (NC_TYPE_BBN_PARTHENOPE);

  ncm_cfg_register_obj (NC_TYPE_HIREION);
  ncm_cfg_register_obj (NC_TYPE_HIREION_CAMB);

  ncm_cfg_register_obj (NC_TYPE_POWSPEC_ML);
  ncm_cfg_register_obj (NC_TYPE_POWSPEC_ML_SPLINE);
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
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DTR_DHR);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DMR_HR);
  ncm_cfg_register_obj (NC_TYPE_DATA_BAO_DVR_DTDH);

  ncm_cfg_register_obj (NC_TYPE_DATA_DIST_MU);

  ncm_cfg_register_obj (NC_TYPE_DATA_HUBBLE);

  ncm_cfg_register_obj (NC_TYPE_DATA_SNIA_COV);

  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_NCOUNT);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_NCOUNTS_GAUSS);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_PSEUDO_COUNTS);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_WL_FACTOR);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_MASS_RICH);
  ncm_cfg_register_obj (NC_TYPE_DATA_CLUSTER_MASS_RICH_COUNT);

  ncm_cfg_register_obj (NC_TYPE_DATA_CMB_SHIFT_PARAM);
  ncm_cfg_register_obj (NC_TYPE_DATA_CMB_DIST_PRIORS);

  ncm_cfg_register_obj (NC_TYPE_XCOR);
  ncm_cfg_register_obj (NC_TYPE_XCOR_SOLVER);
  ncm_cfg_register_obj (NC_TYPE_XCOR_SSC_SIJ);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_COMPONENT);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_RADIAL_KDEP);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_RADIAL_KDEP_GROWTH);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_RADIAL);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_TABLE);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_GAUSS);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_MULTI);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_STUDENT_T);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_POWER_EXP);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_TOPHAT_SMOOTH);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_ANALYTIC_LENSING);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_GAL);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_CLUSTER);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_CLUSTER_TOPHAT);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_TSZ);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_CMB_LENSING);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_CMB_ISW);
  ncm_cfg_register_obj (NC_TYPE_XCOR_KERNEL_WEAK_LENSING);
  ncm_cfg_register_obj (NC_TYPE_DATA_XCOR);
  ncm_cfg_register_obj (NC_TYPE_XCOR_AB);

  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_LKL);
  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_COMMANDER);
  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_LENSING);
  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_PLIK_LITE);
  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_SIMALL);
  ncm_cfg_register_obj (NC_TYPE_DATA_PLANCK_SMICA);

  return;
}

static gsize _functions_initialized = 0;

/**
 * ncm_cfg_register_functions:
 *
 * Registers, once per process, the #NcmMSetFuncList functions of the library models.
 * Called by ncm_cfg_init().
 */
void
ncm_cfg_register_functions (void)
{
  if (g_once_init_enter (&_functions_initialized))
  {
    _nc_hicosmo_register_functions ();
    _nc_hicosmo_de_register_functions ();
    _nc_hiprim_register_functions ();
    _nc_hireion_register_functions ();
    _nc_distance_register_functions ();
    _nc_planck_fi_cor_tt_register_functions ();
    _nc_hicosmo_de_wspline_register_functions ();
    _nc_hicosmo_qspline_register_functions ();
    _nc_galaxy_shape_pop_beta_register_functions ();

    g_once_init_leave (&_functions_initialized, TRUE);
  }

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

  NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Dying slave!\n", _mpi_ctrl.size, _mpi_ctrl.rank);
  exit (0);
}

static gboolean
_ncm_cfg_mpi_cmd_handler (gpointer user_data)
{
  enum buf_type
  {
    input_type,
    ret_type,
    msg_type,
  };
  struct buf_desc
  {
    gpointer obj;
    gpointer buf;
    enum buf_type t;
  };
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
    gboolean end = FALSE;
    gint cmd     = 0;
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

          NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave received a work request, command %d.\n",
                                   _mpi_ctrl.size, _mpi_ctrl.rank, cmd);

          MPI_Recv (input_buf, input_len, input_dtype, NCM_MPI_CTRL_MASTER_ID, NCM_MPI_CTRL_TAG_WORK_INPUT, MPI_COMM_WORLD, &status);
          MPI_Get_count (&status, input_dtype, &input_recv);

          NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Slave received work data: %d-bytes, working...\n", _mpi_ctrl.size, _mpi_ctrl.rank, input_len);

          g_assert_cmpint (input_recv, ==, input_len);

          ncm_mpi_job_unpack_input (mpi_job, input_buf, input);

          bd.obj = ncm_mpi_job_create_return (mpi_job);

          ncm_mpi_job_run (mpi_job, input, bd.obj);
          bd.buf = ncm_mpi_job_pack_return (mpi_job, bd.obj);

          NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Job done, sending result, length %d!\n",
                                   _mpi_ctrl.size, _mpi_ctrl.rank, return_len);

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
        NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Send %d is %s!\n",
                                 _mpi_ctrl.size, _mpi_ctrl.rank, i, done ? "done" : "not done");

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

    NCM_MPI_JOB_DEBUG_PRINT ("#[%3d %3d] Finished command %d, %d requests left%s\n",
                             _mpi_ctrl.size, _mpi_ctrl.rank, cmd, work_ret_request->len, end ? ", exiting!" : ".");

    if (end)
      break;
  }

  if (work_ret_request->len > 0)
  {
    guint i;

    MPI_Waitall (work_ret_request->len, (MPI_Request *) work_ret_request->data, MPI_STATUSES_IGNORE);

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

  ncm_serialize_free (ser);

  return normal_exit;
}

#endif /* HAVE_MPI */

/**
 * ncm_cfg_enable_gsl_err_handler:
 *
 * Restores the GSL error handler that ncm_cfg_init() turned off, so that GSL errors
 * abort.
 */
void
ncm_cfg_enable_gsl_err_handler (void)
{
  g_assert (numcosmo_init);
  gsl_set_error_handler (gsl_err);
}

static guint nreg_model = 0;

/**
 * ncm_cfg_register_obj:
 * @obj: a #GType
 *
 * Registers @obj and initializes its class, so that it can be found by name, for
 * example when deserializing.
 */
void
ncm_cfg_register_obj (GType obj)
{
  g_type_ensure (obj);
  {
    gpointer obj_class = g_type_class_ref (obj);

    g_type_class_unref (obj_class);
    nreg_model++;
  }
}

/**
 * ncm_cfg_mpi_nslaves:
 *
 * Returns: the number of MPI worker ranks, zero outside an MPI launcher.
 */
guint
ncm_cfg_mpi_nslaves (void)
{
  return _mpi_ctrl.nslaves;
}

/**
 * ncm_cfg_set_logfile:
 * @filename: the file name
 *
 * Sends the log messages to @filename, which is truncated. Aborts if it cannot be
 * opened.
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
 * Sends the log messages to @stream.
 */
void
ncm_cfg_set_logstream (FILE *stream)
{
  g_assert (stream != NULL);
  _log_stream = stream;
}

typedef struct _NcmCfgLoggerFuncContainer
{
  NcmCfgLoggerFunc logger;
} NcmCfgLoggerFuncContainer;

static void
_ncm_cfg_log_message_logger (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  NcmCfgLoggerFuncContainer *container = (NcmCfgLoggerFuncContainer *) user_data;

  if (_enable_msg && _log_stream)
    container->logger (message);
}

/* Errors and criticals reach the logger whatever ncm_cfg_logfile() says */
static void
_ncm_cfg_log_error_logger (const gchar *log_domain, GLogLevelFlags log_level, const gchar *message, gpointer user_data)
{
  NcmCfgLoggerFuncContainer *container = (NcmCfgLoggerFuncContainer *) user_data;

  container->logger (message);
}

/**
 * ncm_cfg_set_log_handler:
 * @logger: (scope notified): a logger function
 *
 * Sends the message, info and debug log messages to @logger instead of the log stream.
 */
void
ncm_cfg_set_log_handler (NcmCfgLoggerFunc logger)
{
  static NcmCfgLoggerFuncContainer container = {NULL};

  container.logger = logger;

  _log_msg_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE | G_LOG_LEVEL_INFO | G_LOG_LEVEL_DEBUG, _ncm_cfg_log_message_logger, &container);
}

/**
 * ncm_cfg_set_error_log_handler:
 * @logger: (scope notified): a logger function
 *
 * Sends the error and critical log messages to @logger, also while ncm_cfg_logfile() is
 * off.
 */
void
ncm_cfg_set_error_log_handler (NcmCfgLoggerFunc logger)
{
  static NcmCfgLoggerFuncContainer container = {NULL};

  container.logger = logger;

  _log_err_id = g_log_set_handler (G_LOG_DOMAIN, G_LOG_LEVEL_ERROR | G_LOG_LEVEL_CRITICAL | G_LOG_FLAG_FATAL | G_LOG_FLAG_RECURSION, _ncm_cfg_log_error_logger, &container);
}

/**
 * ncm_cfg_set_openmp_nthreads:
 * @n: number of threads
 *
 * Sets the number of OpenMP threads, if NumCosmo was built with OpenMP.
 */
void
ncm_cfg_set_openmp_nthreads (gint n)
{
#ifdef _OPENMP
  omp_set_num_threads (n);
#endif /* _OPENMP */
}

/**
 * ncm_cfg_set_openblas_nthreads:
 * @n: number of threads
 *
 * Sets the number of OpenBLAS threads, if NumCosmo was built with OpenBLAS.
 */
void
ncm_cfg_set_openblas_nthreads (gint n)
{
#ifdef HAVE_OPENBLAS_SET_NUM_THREADS
  openblas_set_num_threads (n);
  goto_set_num_threads (n);
#endif /* HAVE_OPENBLAS_SET_NUM_THREADS */
}

/**
 * ncm_cfg_set_blis_nthreads:
 * @n: number of threads
 *
 * Sets the number of BLIS threads, if NumCosmo was built with BLIS.
 */
void
ncm_cfg_set_blis_nthreads (gint n)
{
#ifdef HAVE_BLIS
  bli_thread_set_num_threads (n);
#endif /* HAVE_BLIS */
}

/**
 * ncm_cfg_set_mkl_nthreads:
 * @n: number of threads
 *
 * Sets the number of MKL threads, if NumCosmo was built with MKL.
 */
void
ncm_cfg_set_mkl_nthreads (gint n)
{
#ifdef HAVE_MKL_SET_NUM_THREADS
  MKL_Set_Num_Threads (n);
#endif /* HAVE_MKL_SET_NUM_THREADS */
}

/**
 * ncm_cfg_logfile:
 * @on: whether messages are logged
 *
 * Turns the log messages on or off.
 */
void
ncm_cfg_logfile (gboolean on)
{
  _enable_msg = on;
}

/**
 * ncm_cfg_logfile_flush:
 * @on: whether to flush
 *
 * Turns on or off flushing the log stream after each message.
 */
void
ncm_cfg_logfile_flush (gboolean on)
{
  _enable_msg_flush = on;
}

/**
 * ncm_cfg_logfile_flush_now:
 *
 * Flushes the log stream.
 */
void
ncm_cfg_logfile_flush_now (void)
{
  fflush (_log_stream);
}

/**
 * ncm_message_str:
 * @msg: a string
 *
 * Logs @msg as a message of the NumCosmo log domain.
 */
void
ncm_message_str (const gchar *msg)
{
  g_log (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, "%s", msg);
}

/**
 * ncm_message:
 * @msg: a printf format string
 * @...: arguments for @msg
 *
 * Logs the formatted message as a message of the NumCosmo log domain.
 */
void
ncm_message (const gchar *msg, ...)
{
  va_list ap;

  va_start (ap, msg);
  g_logv (G_LOG_DOMAIN, G_LOG_LEVEL_MESSAGE, msg, ap);
  va_end (ap);
}

/**
 * ncm_string_ww:
 * @msg: a string
 * @first: prefix of the first line
 * @rest: prefix of the other lines
 * @ncols: number of columns
 *
 * Wraps @msg at spaces into lines of at most @ncols columns, each ending with a newline.
 *
 * Returns: (transfer full): the wrapped string.
 */
gchar *
ncm_string_ww (const gchar *msg, const gchar *first, const gchar *rest, guint ncols)
{
  gchar **msg_split   = g_strsplit (msg, " ", 0);
  guint size_print    = strlen (msg_split[0]);
  guint first_size    = strlen (first);
  guint rest_size     = strlen (rest);
  guint msg_len       = strlen (msg);
  guint msg_final_len = msg_len + first_size + msg_len / ncols * rest_size;
  GString *msg_ww     = g_string_sized_new (msg_final_len);
  guint i             = 1;

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
 * @msg: a string
 * @first: prefix of the first line
 * @rest: prefix of the other lines
 * @ncols: number of columns
 *
 * Logs @msg wrapped by ncm_string_ww().
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
 * Logs a separator line.
 */
void
ncm_cfg_msg_sepa (void)
{
  g_message ("#----------------------------------------------------------------------------------\n");
}

/**
 * ncm_cfg_get_fullpath:
 * @filename: a printf format string
 * @...: arguments for @filename
 *
 * Returns: (transfer full): the path of the formatted file name inside `~/.numcosmo`.
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
 * ncm_cfg_get_fullpath_base:
 *
 * Returns: (transfer none): the path of `~/.numcosmo`.
 */
const gchar *
ncm_cfg_get_fullpath_base (void)
{
  g_assert (numcosmo_init);

  return numcosmo_path;
}

/**
 * ncm_cfg_keyfile_to_arg:
 * @kfile: a #GKeyFile
 * @group_name: group name
 * @entries: a %NULL-terminated array of #GOptionEntry
 * @argv: an array of strings
 * @argc: the number of strings in @argv
 *
 * Appends to @argv, starting at position *@argc, the command-line arguments equivalent
 * to the keys of @group_name in @kfile that match @entries, and updates *@argc. Boolean
 * keys become a flag when true, empty values are skipped and list keys repeat the
 * option. @argv must have room for them. Aborts if a key cannot be parsed.
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
        if ((entries[i].arg == G_OPTION_ARG_STRING_ARRAY) || (entries[i].arg == G_OPTION_ARG_FILENAME_ARRAY))
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
              argv[argc[0]++] = g_strdup_printf ("--%s", entries[i].long_name);

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
 * @str: a string
 *
 * Returns: (transfer full): @str wrapped at 80 columns below a separator line, as a
 * #GKeyFile comment.
 */
gchar *
ncm_cfg_string_to_comment (const gchar *str)
{
  g_assert (str != NULL);
  {
    gchar *desc_ww = ncm_string_ww (str, "  ", "  ", 80);
    gchar *desc    = g_strdup_printf ("###############################################################################\n\n%s\n", desc_ww);

    g_free (desc_ww);

    return desc;
  }
}

/**
 * ncm_cfg_entries_to_keyfile:
 * @kfile: a #GKeyFile
 * @group_name: group name
 * @entries: a %NULL-terminated array of #GOptionEntry
 *
 * Writes the current values of @entries as keys of @group_name in @kfile, each commented
 * with the entry description. Callback entries are skipped.
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
        gboolean arg_b = ((gboolean *) entries[i].arg_data)[0];

        g_key_file_set_boolean (kfile, group_name, entries[i].long_name, arg_b);
        break;
      }
      case G_OPTION_ARG_STRING:
      case G_OPTION_ARG_FILENAME:
      {
        gchar **arg_s = (gchar **) entries[i].arg_data;

        g_key_file_set_string (kfile, group_name, entries[i].long_name, *arg_s != NULL ? *arg_s : "");
        break;
      }
      case G_OPTION_ARG_STRING_ARRAY:
      case G_OPTION_ARG_FILENAME_ARRAY:
      {
        const gchar ***arg_as = (const gchar ***) entries[i].arg_data;
        gchar ***arg_cas      = (gchar ***) entries[i].arg_data;

        if (*arg_cas != NULL)
          g_key_file_set_string_list (kfile, group_name, entries[i].long_name, *arg_as, g_strv_length (*arg_cas));
        else
          g_key_file_set_string_list (kfile, group_name, entries[i].long_name, NULL, 0);

        break;
      }
      case G_OPTION_ARG_INT:
      {
        gint arg_i = ((gint *) entries[i].arg_data)[0];

        g_key_file_set_integer (kfile, group_name, entries[i].long_name, arg_i);
        break;
      }
      case G_OPTION_ARG_INT64:
      {
        gint64 arg_l = ((gint64 *) entries[i].arg_data)[0];

        g_key_file_set_int64 (kfile, group_name, entries[i].long_name, arg_l);
        break;
      }
      case G_OPTION_ARG_DOUBLE:
      {
        gdouble arg_d = ((double *) entries[i].arg_data)[0];

        g_key_file_set_double (kfile, group_name, entries[i].long_name, arg_d);
        break;
      }
      case G_OPTION_ARG_CALLBACK:
      default:
        skip_comment = TRUE;
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
 * @enum_type: an enumeration #GType
 * @id_name_nick: a value, name or nick
 *
 * Looks up @id_name_nick in @enum_type: a decimal integer is taken as the value,
 * otherwise as the name, then as the nick.
 *
 * Returns: (transfer none) (nullable): the #GEnumValue, or %NULL if not found.
 */
const GEnumValue *
ncm_cfg_get_enum_by_id_name_nick (GType enum_type, const gchar *id_name_nick)
{
  g_assert (id_name_nick != NULL);
  {
    GEnumValue *res        = NULL;
    gchar *endptr          = NULL;
    gint64 id              = g_ascii_strtoll (id_name_nick, &endptr, 10);
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
    {
      res = g_enum_get_value (enum_class, id);
    }

    g_type_class_unref (enum_class);

    return res;
  }
}

/**
 * ncm_cfg_enum_get_value:
 * @enum_type: an enumeration #GType
 * @n: the value
 *
 * Returns: (transfer none) (nullable): the #GEnumValue of value @n, or %NULL if not found.
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
 * @enum_type: an enumeration #GType
 * @header: header string
 *
 * Prints to standard output a table of the values, names and nicks of every member of
 * @enum_type, in declaration order.
 */
void
ncm_cfg_enum_print_all (GType enum_type, const gchar *header)
{
  GEnumClass *enum_class;
  GEnumValue *snia;
  gint i            = 0;
  gint name_max_len = 4;
  gint nick_max_len = 4;
  gint pad;

  g_assert (G_TYPE_IS_ENUM (enum_type));

  enum_class = g_type_class_ref (enum_type);

  for (i = 0; i < (gint) enum_class->n_values; i++)
  {
    snia         = &enum_class->values[i];
    name_max_len = GSL_MAX (name_max_len, (gint) strlen (snia->value_name));
    nick_max_len = GSL_MAX (nick_max_len, (gint) strlen (snia->value_nick));
  }

  printf ("# %s:\n", header);
  pad = 10 + name_max_len + nick_max_len;
  printf ("#");

  while (pad-- != 0)
    printf ("-");

  printf ("#\n");
  printf ("# Id | %-*s | %-*s |\n", name_max_len, "Name", nick_max_len, "Nick");

  for (i = 0; i < (gint) enum_class->n_values; i++)
  {
    snia = &enum_class->values[i];
    printf ("# %02d | %-*s | %-*s |\n", snia->value, name_max_len, snia->value_name, nick_max_len, snia->value_nick);
  }

  pad = 10 + name_max_len + nick_max_len;
  printf ("#");

  while (pad-- != 0)
    printf ("-");

  printf ("#\n");

  g_type_class_unref (enum_class);
}

G_LOCK_DEFINE_STATIC (fftw_saveload_lock);

G_LOCK_DEFINE_STATIC (fftw_plan_lock);

/*
 * One wisdom file per MPI rank. FFTW's wisdom is global and cumulative, so it is loaded
 * once per process and saved only when it changed. Protected by fftw_saveload_lock.
 */
static gboolean _wisdom_loaded_d = FALSE;
static gboolean _wisdom_loaded_f = FALSE;
static gchar *_wisdom_saved_d    = NULL;
static gchar *_wisdom_saved_f    = NULL;

/**
 * ncm_cfg_lock_plan_fftw:
 *
 * Locks the global lock that serializes FFTW planning, which is not thread-safe.
 */
void
ncm_cfg_lock_plan_fftw (void)
{
  G_LOCK (fftw_plan_lock);
}

/**
 * ncm_cfg_unlock_plan_fftw:
 *
 * Unlocks the lock taken by ncm_cfg_lock_plan_fftw().
 */
void
ncm_cfg_unlock_plan_fftw (void)
{
  G_UNLOCK (fftw_plan_lock);
}

static void _ncm_cfg_load_fftw_wisdom (void);
static void _ncm_cfg_save_fftw_wisdom (void);

G_LOCK_DEFINE_STATIC (fftw_keys_lock);

static GHashTable *_fftw_planned_keys = NULL;

/**
 * ncm_cfg_fftw_plan_begin: (skip)
 * @key: a printf format string naming what is planned
 * @...: arguments for @key
 *
 * Starts creating FFTW plans: loads the FFTW wisdom of this MPI rank, once per process, from
 * `~/.numcosmo/ncm_cfg_wisdom_rank<rank>.fftw3` (and `.fftw3f`), and takes the
 * planning lock, see ncm_cfg_lock_plan_fftw(). @key identifies the plans: the caller and
 * everything that makes a plan different, such as the transform sizes and kinds and the
 * number of transforms; the current default planner flag is added to it. It only tells
 * whether the process planned the same before, and the wisdom file is the same for every
 * key. Pass the return value to ncm_cfg_fftw_plan_end().
 *
 * Returns: whether @key is planned for the first time in this process.
 */
gboolean
ncm_cfg_fftw_plan_begin (const gchar *key, ...)
{
  gchar *key_str;
  gboolean first;
  va_list ap;

  {
    gchar *site_key;

    va_start (ap, key);
    site_key = g_strdup_vprintf (key, ap);
    va_end (ap);

    /* A new planner flag plans anew, so it is part of the key */
    key_str = g_strdup_printf ("%s|%s", site_key, ncm_cfg_get_fftw_default_flag_str ());
    g_free (site_key);
  }

  _ncm_cfg_load_fftw_wisdom ();

  G_LOCK (fftw_keys_lock);

  if (_fftw_planned_keys == NULL)
    _fftw_planned_keys = g_hash_table_new_full (g_str_hash, g_str_equal, g_free, NULL);

  first = !g_hash_table_contains (_fftw_planned_keys, key_str);

  if (first)
    g_hash_table_add (_fftw_planned_keys, key_str);
  else
    g_free (key_str);

  G_UNLOCK (fftw_keys_lock);

  ncm_cfg_lock_plan_fftw ();

  return first;
}

/**
 * ncm_cfg_fftw_plan_end: (skip)
 * @first: the value returned by ncm_cfg_fftw_plan_begin()
 *
 * Releases the planning lock and, when @first is %TRUE, rewrites the wisdom file if the
 * wisdom changed. A key planned before adds no wisdom, so its save, which exports the whole
 * wisdom to compare it, is skipped. Wisdom is neither loaded nor saved under FFTW_ESTIMATE.
 */
void
ncm_cfg_fftw_plan_end (gboolean first)
{
  ncm_cfg_unlock_plan_fftw ();

  if (first)
    _ncm_cfg_save_fftw_wisdom ();
}

/**
 * ncm_cfg_fftw_plan_destroy: (skip)
 * @plan: (nullable): a double-precision FFTW plan
 *
 * Destroys @plan holding the planning lock, see ncm_cfg_lock_plan_fftw(): FFTW's plan
 * destruction is not thread-safe either. Does nothing for %NULL, so it can be the free
 * function of a #GPtrArray of plans.
 */
void
ncm_cfg_fftw_plan_destroy (gpointer plan)
{
  if (plan == NULL)
    return;

  ncm_cfg_lock_plan_fftw ();
  fftw_destroy_plan (plan);
  ncm_cfg_unlock_plan_fftw ();
}

/**
 * ncm_cfg_fftwf_plan_destroy: (skip)
 * @plan: (nullable): a single-precision FFTW plan
 *
 * Same as ncm_cfg_fftw_plan_destroy() for single precision. Aborts if NumCosmo was built
 * without single-precision FFTW.
 */
void
ncm_cfg_fftwf_plan_destroy (gpointer plan)
{
  if (plan == NULL)
    return;

#ifdef HAVE_FFTW3F
  ncm_cfg_lock_plan_fftw ();
  fftwf_destroy_plan (plan);
  ncm_cfg_unlock_plan_fftw ();
#else /* HAVE_FFTW3F */
  g_error ("ncm_cfg_fftwf_plan_destroy: NumCosmo was built without single-precision FFTW.");
#endif /* HAVE_FFTW3F */
}

/*
 * Imports the FFTW wisdom of this MPI rank, once per process, from
 * ~/.numcosmo/ncm_cfg_wisdom_rank<rank>.fftw3 (and .fftw3f). Does nothing under
 * FFTW_ESTIMATE, which uses no wisdom. Thread-safe.
 */
static void
_ncm_cfg_load_fftw_wisdom (void)
{
  gchar *file, *file_ext;
  gchar *full_filename;

  g_assert (numcosmo_init);

  /* FFTW_ESTIMATE neither consumes nor produces useful wisdom; skip the file I/O. */
  if (ncm_cfg_get_fftw_default_flag () == FFTW_ESTIMATE)
    return;

  G_LOCK (fftw_saveload_lock);

  if (_wisdom_loaded_d && _wisdom_loaded_f)
  {
    /* Already loaded once this process -- FFTW's wisdom registry is
     * global and cumulative, so re-parsing the same file again would
     * teach it nothing new. */
    G_UNLOCK (fftw_saveload_lock);

    return;
  }

  file = g_strdup_printf ("ncm_cfg_wisdom_rank%d", _mpi_ctrl.rank);

  file_ext      = g_strdup_printf ("%s.fftw3", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  if (!_wisdom_loaded_d)
  {
    if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
      fftw_import_wisdom_from_filename (full_filename);

    _wisdom_loaded_d = TRUE;
  }

#ifdef HAVE_FFTW3F
  g_free (file_ext);
  g_free (full_filename);

  file_ext      = g_strdup_printf ("%s.fftw3f", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  if (!_wisdom_loaded_f)
  {
    if (g_file_test (full_filename, G_FILE_TEST_EXISTS))
      fftwf_import_wisdom_from_filename (full_filename);

    _wisdom_loaded_f = TRUE;
  }

#else
  _wisdom_loaded_f = TRUE; /* no single-precision FFTW3 build, nothing to load */
#endif

  g_free (file);
  g_free (file_ext);
  g_free (full_filename);

  G_UNLOCK (fftw_saveload_lock);
}

/*
 * Writes the FFTW wisdom of this MPI rank to the file _ncm_cfg_load_fftw_wisdom() reads,
 * only if it changed since the last save. Does nothing under FFTW_ESTIMATE. Thread-safe.
 */
static void
_ncm_cfg_save_fftw_wisdom (void)
{
  gchar *file, *file_ext;
  gchar *full_filename;

  g_assert (numcosmo_init);

  /* FFTW_ESTIMATE neither consumes nor produces useful wisdom; skip the file I/O. */
  if (ncm_cfg_get_fftw_default_flag () == FFTW_ESTIMATE)
    return;

  G_LOCK (fftw_saveload_lock);

  file = g_strdup_printf ("ncm_cfg_wisdom_rank%d", _mpi_ctrl.rank);

  file_ext      = g_strdup_printf ("%s.fftw3", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  {
    char *wisdom_str = fftw_export_wisdom_to_string ();

    if (wisdom_str != NULL)
    {
      if ((_wisdom_saved_d != NULL) && g_str_equal (_wisdom_saved_d, wisdom_str))
      {
        /* Nothing learned since the last save -- skip the rewrite. */
        g_free (wisdom_str);
      }
      else
      {
        gssize len  = strlen (wisdom_str);
        gboolean OK = FALSE;

#if GLIB_CHECK_VERSION (2, 66, 0)
        OK = g_file_set_contents_full (full_filename, wisdom_str, len,
                                       G_FILE_SET_CONTENTS_CONSISTENT,
                                       0666, NULL);
#else /* GLIB_CHECK_VERSION (2, 66, 0) */
        OK = g_file_set_contents (full_filename, wisdom_str, len, NULL);
#endif /* GLIB_CHECK_VERSION (2, 66, 0) */
        g_assert (OK);
        g_free (_wisdom_saved_d);
        _wisdom_saved_d = wisdom_str; /* keep as the new comparison baseline */
      }
    }
  }

#ifdef HAVE_FFTW3F
  g_free (file_ext);
  g_free (full_filename);

  file_ext      = g_strdup_printf ("%s.fftw3f", file);
  full_filename = g_build_filename (numcosmo_path, file_ext, NULL);

  {
    char *wisdom_str = fftwf_export_wisdom_to_string ();

    if (wisdom_str != NULL)
    {
      if ((_wisdom_saved_f != NULL) && g_str_equal (_wisdom_saved_f, wisdom_str))
      {
        /* Nothing learned since the last save -- skip the rewrite. */
        g_free (wisdom_str);
      }
      else
      {
        gssize len  = strlen (wisdom_str);
        gboolean OK = FALSE;

#if GLIB_CHECK_VERSION (2, 66, 0)
        OK = g_file_set_contents_full (full_filename, wisdom_str, len,
                                       G_FILE_SET_CONTENTS_CONSISTENT,
                                       0666, NULL);
#else /* GLIB_CHECK_VERSION (2, 66, 0) */
        OK = g_file_set_contents (full_filename, wisdom_str, len, NULL);
#endif /* GLIB_CHECK_VERSION (2, 66, 0) */

        g_assert (OK);
        g_free (_wisdom_saved_f);
        _wisdom_saved_f = wisdom_str; /* keep as the new comparison baseline */
      }
    }
  }
#endif

  g_free (file);
  g_free (file_ext);
  g_free (full_filename);

  G_UNLOCK (fftw_saveload_lock);
}

/**
 * ncm_cfg_exists:
 * @filename: a printf format string
 * @...: arguments for @filename
 *
 * Returns: whether the formatted file name exists inside `~/.numcosmo`.
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
  exists        = g_file_test (full_filename, G_FILE_TEST_EXISTS);
  g_free (file);
  g_free (full_filename);

  return exists;
}

/**
 * ncm_cfg_get_data_filename:
 * @filename: a path relative to the data directory
 * @must_exist: whether to abort if @filename is not found
 *
 * Looks for @filename in the `data` directory under, in order, the path in
 * #NCM_CFG_DATA_DIR_ENV, the installed package data directory and the source directory.
 *
 * Returns: (transfer full) (nullable): the full path, or %NULL if not found and
 * @must_exist is %FALSE.
 */
gchar *
ncm_cfg_get_data_filename (const gchar *filename, gboolean must_exist)
{
  const gchar *data_dir = g_getenv (NCM_CFG_DATA_DIR_ENV);
  gchar *full_filename  = NULL;

  if (data_dir != NULL)
  {
    full_filename = g_build_filename (data_dir, "data", filename, NULL);

    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
      g_clear_pointer (&full_filename, g_free);
  }

  if (full_filename == NULL)
  {
    full_filename = g_build_filename (PACKAGE_DATA_DIR, "data", filename, NULL);

    if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
      g_clear_pointer (&full_filename, g_free);
  }

  if (full_filename == NULL)
    full_filename = g_build_filename (PACKAGE_SOURCE_DIR, "data", filename, NULL);

  if (!g_file_test (full_filename, G_FILE_TEST_EXISTS))
  {
    if (must_exist)
      g_error ("ncm_cfg_get_data_filename: cannot find `%s'.", filename);
    else
      g_clear_pointer (&full_filename, g_free);
  }

  return full_filename;
}

/**
 * ncm_cfg_get_data_directory:
 *
 * Looks for the `data` directory in the same places as ncm_cfg_get_data_filename().
 * Aborts if none exists.
 *
 * Returns: (transfer full): the path of the data directory.
 */
gchar *
ncm_cfg_get_data_directory (void)
{
  const gchar *data_dir = g_getenv (NCM_CFG_DATA_DIR_ENV);
  gchar *full_directory = NULL;

  if (data_dir != NULL)
  {
    full_directory = g_build_filename (data_dir, "data", NULL);

    if (!g_file_test (full_directory, G_FILE_TEST_IS_DIR))
      g_clear_pointer (&full_directory, g_free);
  }

  if (full_directory == NULL)
  {
    full_directory = g_build_filename (PACKAGE_DATA_DIR, "data", NULL);

    if (!g_file_test (full_directory, G_FILE_TEST_IS_DIR))
      g_clear_pointer (&full_directory, g_free);
  }

  if (full_directory == NULL)
    full_directory = g_build_filename (PACKAGE_SOURCE_DIR, "data", NULL);

  if (!g_file_test (full_directory, G_FILE_TEST_IS_DIR))
  {
    g_clear_pointer (&full_directory, g_free);
    g_error ("ncm_cfg_get_data_directory: cannot determine data directory.");
  }


  return full_directory;
}

/**
 * ncm_cfg_command_line:
 * @argv: array of strings
 * @argc: number of strings in @argv
 *
 * Joins @argv with spaces, quoting with single quotes the arguments that contain a space.
 *
 * Returns: (transfer full): the command line.
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

  argv_size = strlen (argv[0]);
  memcpy (full_cmd_line_ptr, argv[0], argv_size);
  full_cmd_line_ptr = &full_cmd_line_ptr[argv_size];

  for (i = 1; i < argc; i++)
  {
    gboolean has_space = FALSE;

    full_cmd_line_ptr[0] = ' ';
    full_cmd_line_ptr++;
    argv_size = strlen (argv[i]);
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
 * ncm_cfg_array_set_variant: (skip)
 * @a: a #GArray
 * @var: a #GVariant of array type
 *
 * Resizes @a to the length of @var and copies its elements, which must have the element
 * size of @a.
 */
void
ncm_cfg_array_set_variant (GArray *a, GVariant *var)
{
  gsize esize        = g_array_get_element_size (a);
  gsize n_elements   = 0;
  gconstpointer data = g_variant_get_fixed_array (var, &n_elements, esize);

  g_array_set_size (a, n_elements);
  memcpy (a->data, data, n_elements * esize);
}

/**
 * ncm_cfg_array_to_variant: (skip)
 * @a: a #GArray
 * @etype: the element type
 *
 * Creates a #GVariant array of @etype sharing the data of @a, which it holds a reference
 * to.
 *
 * Returns: (transfer full): the #GVariant.
 */
GVariant *
ncm_cfg_array_to_variant (GArray *a, const GVariantType *etype)
{
  gconstpointer data  = a->data;
  gsize esize         = g_array_get_element_size (a);
  GVariantType *atype = g_variant_type_new_array (etype);
  GVariant *vvar      = g_variant_new_from_data (atype,
                                                 data,
                                                 esize * a->len,
                                                 TRUE,
                                                 (GDestroyNotify) & g_array_unref,
                                                 g_array_ref (a));

  g_variant_type_free (atype);

  return g_variant_ref_sink (vvar);
}

static guint __fftw_default_flags = FFTW_MEASURE;
static gdouble __fftw_timelimit   = 60.0;

/**
 * ncm_cfg_set_fftw_default_flag:
 * @flag: FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE
 * @timeout: planner time limit in seconds
 * @error: a #GError location
 *
 * Sets the planner flag used for new FFTW plans and FFTW's planner time limit; a negative
 * @timeout means no limit. Sets @error, or aborts if @error is %NULL, for any other
 * @flag.
 */
void
ncm_cfg_set_fftw_default_flag (guint flag, const gdouble timeout, GError **error)
{
  switch (flag)
  {
    case FFTW_ESTIMATE:
    case FFTW_MEASURE:
    case FFTW_PATIENT:
    case FFTW_EXHAUSTIVE:
      break;

    default:
      ncm_util_set_or_call_error (error,
                                  NCM_CFG_ERROR,
                                  NCM_CFG_ERROR_INVALID_FFTW_FLAG,
                                  "Invalid FFTW flag '%d'", flag);

      return;
  }

  __fftw_default_flags = flag;
  __fftw_timelimit     = timeout;

  fftw_set_timelimit (timeout);
#ifdef HAVE_FFTW3F
  fftwf_set_timelimit (timeout);
#endif /* HAVE_FFTW3F */
}

/**
 * ncm_cfg_set_fftw_default_flag_str:
 * @flag_str: "estimate", "measure", "patient" or "exhaustive", in any case
 * @timeout: planner time limit in seconds
 * @error: a #GError location
 *
 * Same as ncm_cfg_set_fftw_default_flag() with the flag named by @flag_str. Sets @error,
 * or aborts if @error is %NULL, for any other string.
 */
void
ncm_cfg_set_fftw_default_flag_str (const gchar *flag_str, const gdouble timeout, GError **error)
{
  guint flag = 0;

  if (g_ascii_strcasecmp (flag_str, "estimate") == 0)
  {
    flag = FFTW_ESTIMATE;
  }
  else if (g_ascii_strcasecmp (flag_str, "measure") == 0)
  {
    flag = FFTW_MEASURE;
  }
  else if (g_ascii_strcasecmp (flag_str, "patient") == 0)
  {
    flag = FFTW_PATIENT;
  }
  else if (g_ascii_strcasecmp (flag_str, "exhaustive") == 0)
  {
    flag = FFTW_EXHAUSTIVE;
  }
  else
  {
    ncm_util_set_or_call_error (error,
                                NCM_CFG_ERROR,
                                NCM_CFG_ERROR_INVALID_FFTW_FLAG_STRING,
                                "Invalid FFTW flag string '%s'", flag_str);

    return;
  }

  ncm_cfg_set_fftw_default_flag (flag, timeout, error);
}

/**
 * ncm_cfg_set_fftw_default_from_env:
 * @fallback_flag: the flag used if `NCM_FFTW_PLANNER` is not set
 * @fallback_timeout: the time limit used if `NCM_FFTW_PLANNER_TIMELIMIT` is not set
 * @error: a #GError location
 *
 * Same as ncm_cfg_set_fftw_default_flag() with the flag named in `NCM_FFTW_PLANNER`
 * and the time limit in `NCM_FFTW_PLANNER_TIMELIMIT`. Sets @error, or aborts if @error
 * is %NULL, if either is invalid.
 */
void
ncm_cfg_set_fftw_default_from_env (guint fallback_flag, const gdouble fallback_timeout, GError **error)
{
  const gchar *fftw_planner_env   = g_getenv ("NCM_FFTW_PLANNER");
  const gchar *fftw_timelimit_env = g_getenv ("NCM_FFTW_PLANNER_TIMELIMIT");
  gdouble timeout                 = fallback_timeout;

  if (fftw_timelimit_env != NULL)
  {
    gchar *endptr;
    gdouble timelimit = g_ascii_strtod (fftw_timelimit_env, &endptr);

    if (endptr == fftw_timelimit_env)
    {
      ncm_util_set_or_call_error (error,
                                  NCM_CFG_ERROR,
                                  NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT,
                                  "Invalid FFTW planner timelimit '%s'", fftw_timelimit_env);

      return;
    }

    timeout = timelimit;
  }

  if (fftw_planner_env != NULL)
    ncm_cfg_set_fftw_default_flag_str (fftw_planner_env, timeout, error);
  else
    ncm_cfg_set_fftw_default_flag (fallback_flag, timeout, error);
}

/**
 * ncm_cfg_set_fftw_default_from_env_str:
 * @fallback_flag_str: the flag name used if `NCM_FFTW_PLANNER` is not set
 * @fallback_timeout: the time limit used if `NCM_FFTW_PLANNER_TIMELIMIT` is not set
 * @error: a #GError location
 *
 * Same as ncm_cfg_set_fftw_default_from_env() with the fallback flag given by name.
 */
void
ncm_cfg_set_fftw_default_from_env_str (const gchar *fallback_flag_str, const gdouble fallback_timeout, GError **error)
{
  const gchar *fftw_planner_env   = g_getenv ("NCM_FFTW_PLANNER");
  const gchar *fftw_timelimit_env = g_getenv ("NCM_FFTW_PLANNER_TIMELIMIT");
  gdouble timeout                 = fallback_timeout;

  if (fftw_timelimit_env != NULL)
  {
    gchar *endptr;
    gdouble timelimit = g_ascii_strtod (fftw_timelimit_env, &endptr);

    if (endptr == fftw_timelimit_env)
    {
      ncm_util_set_or_call_error (error,
                                  NCM_CFG_ERROR,
                                  NCM_CFG_ERROR_INVALID_FFTW_TIMELIMIT,
                                  "Invalid FFTW planner timelimit '%s'", fftw_timelimit_env);

      return;
    }

    timeout = timelimit;
  }

  if (fftw_planner_env != NULL)
    ncm_cfg_set_fftw_default_flag_str (fftw_planner_env, timeout, error);
  else
    ncm_cfg_set_fftw_default_flag_str (fallback_flag_str, timeout, error);
}

/**
 * ncm_cfg_get_fftw_default_flag:
 *
 * Returns: the planner flag used for new FFTW plans.
 */
guint
ncm_cfg_get_fftw_default_flag (void)
{
  return __fftw_default_flags;
}

/**
 * ncm_cfg_get_fftw_default_flag_str:
 *
 * Returns: (transfer none): the name of the planner flag used for new FFTW plans.
 */
const gchar *
ncm_cfg_get_fftw_default_flag_str (void)
{
  const gchar *flag_str = NULL;

  switch (__fftw_default_flags)
  {
    case FFTW_ESTIMATE:
      flag_str = "estimate";
      break;
    case FFTW_MEASURE:
      flag_str = "measure";
      break;
    case FFTW_PATIENT:
      flag_str = "patient";
      break;
    case FFTW_EXHAUSTIVE:
      flag_str = "exhaustive";
      break;
    default:                   /* LCOV_EXCL_LINE */
      g_assert_not_reached (); /* LCOV_EXCL_LINE */
  }

  return flag_str;
}

/**
 * ncm_cfg_get_fftw_timelimit:
 *
 * Returns: the planner time limit last set through ncm_cfg_set_fftw_default_flag(), in
 * seconds; negative means no limit.
 */
gdouble
ncm_cfg_get_fftw_timelimit (void)
{
  return __fftw_timelimit;
}

/**
 * ncm_cfg_get_version:
 * @major: (out) (optional): the major version
 * @minor: (out) (optional): the minor version
 * @micro: (out) (optional): the micro version
 *
 * Returns: the version as $10^4\,\mathrm{major} + 10^2\,\mathrm{minor} + \mathrm{micro}$.
 */
guint
ncm_cfg_get_version (guint *major, guint *minor, guint *micro)
{
  if (major != NULL)
    *major = NUMCOSMO_MAJOR_VERSION;

  if (minor != NULL)
    *minor = NUMCOSMO_MINOR_VERSION;

  if (micro != NULL)
    *micro = NUMCOSMO_MICRO_VERSION;

  return NUMCOSMO_VERSION_NUMBER;
}

/**
 * ncm_cfg_get_version_string:
 *
 * Returns: (transfer full): the version as "major.minor.micro".
 */
gchar *
ncm_cfg_get_version_string (void)
{
  return g_strdup_printf ("%d.%d.%d", NUMCOSMO_MAJOR_VERSION, NUMCOSMO_MINOR_VERSION, NUMCOSMO_MICRO_VERSION);
}

/**
 * ncm_cfg_version_check:
 * @major: major version
 * @minor: minor version
 * @micro: micro version
 *
 * Returns: whether the library version is at least @major.@minor.@micro.
 */
gboolean
ncm_cfg_version_check (guint major, guint minor, guint micro)
{
  /* We can remove the casing once MAJOR is larger than 0 */
  return ((gint) NUMCOSMO_MAJOR_VERSION > (gint) major) ||
         ((NUMCOSMO_MAJOR_VERSION == major) && (NUMCOSMO_MINOR_VERSION > minor)) ||
         ((NUMCOSMO_MAJOR_VERSION == major) && (NUMCOSMO_MINOR_VERSION == minor) && (NUMCOSMO_MICRO_VERSION >= micro));
}

/**
 * ncm_cfg_get_commit_hash:
 *
 * Returns: (transfer none): the git commit the library was built from.
 */
const gchar *
ncm_cfg_get_commit_hash (void)
{
  return NUMCOSMO_GIT_COMMIT;
}

