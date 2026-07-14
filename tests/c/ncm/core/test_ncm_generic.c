/***************************************************************************
 *            test_ncm_generic.c
 *
 *  Wed July 26 12:04:44 2017
 *  Copyright  2017  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_generic.c
 *
 * Copyright (C) 2017 - Sandro Dias Pinto Vitenti
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

void test_ncm_data_basic (void);
void test_ncm_data_funnel_basic (void);
void test_ncm_data_gaussmix2d_basic (void);
void test_ncm_data_rosenbrock_basic (void);
void test_ncm_dataset_basic (void);
void test_ncm_fftlog_basic (void);
void test_ncm_sphere_nn (void);
void test_ncm_timer_basic (void);
void test_ncm_sbessel_ode_solver_basic (void);
void test_ncm_sbessel_integrator_gl_basic (void);
void test_ncm_sbessel_integrator_fftl_basic (void);
void test_ncm_sbessel_integrator_levin_basic (void);
void test_ncm_fftlog_sbessel_j_basic (void);
void test_ncm_fftlog_sbessel_jljm_basic (void);
void test_ncm_bootstrap_basic (void);
void test_ncm_stats_vec_basic (void);
void test_ncm_mpi_job_basic (void);
void test_ncm_mpi_job_test_basic (void);
void test_ncm_mpi_job_fit_basic (void);
void test_ncm_mpi_job_mcmc_basic (void);
void test_ncm_mpi_job_feval_basic (void);
void test_ncm_powspec_spline2d_basic (void);
void test_ncm_pln1d_basic (void);

void test_nc_data_cluster_mass_rich_basic (void);
void test_nc_data_cluster_mass_rich_count_basic (void);
void test_nc_data_cluster_wl_basic (void);
void test_nc_de_cont_basic (void);
void test_nc_distance_basic (void);
void test_nc_growth_func_basic (void);
void test_nc_transfer_func_basic (void);
void test_nc_galaxy_sd_obs_redshift_basic (void);
void test_nc_galaxy_sd_obs_redshift_gauss_basic (void);
void test_nc_galaxy_sd_obs_redshift_spec_basic (void);
void test_nc_galaxy_sd_obs_redshift_pz_basic (void);
void test_nc_galaxy_sd_position_basic (void);
void test_nc_galaxy_sd_position_flat_basic (void);
void test_nc_galaxy_sd_true_redshift_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_y1_source_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_y1_lens_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_y10_source_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_y10_lens_basic (void);
void test_nc_galaxy_sd_true_redshift_lsst_srd_from_type_basic (void);
void test_nc_galaxy_sd_obs_redshift_gauss_lsst_srd_bins_basic (void);
void test_nc_galaxy_sd_shape_basic (void);
void test_nc_galaxy_sd_shape_hsm_gauss_global_basic (void);
void test_nc_galaxy_sd_shape_hsm_gauss_basic (void);
void test_nc_galaxy_sd_z_proxy_basic (void);
void test_nc_galaxy_sd_z_proxy_dirac_basic (void);
void test_nc_galaxy_sd_z_proxy_gauss_basic (void);
void test_nc_galaxy_wl_obs_basic (void);
void test_nc_halo_position_basic (void);
void test_nc_hicosmo_qgw_basic (void);
void test_nc_hipert_adiab_basic (void);
void test_nc_hipert_em_basic (void);
void test_nc_hipert_gw_basic (void);
void test_nc_hipert_two_fluids_basic (void);
void test_nc_hiprim_two_fluids_basic (void);
void test_nc_halo_mass_summary_basic (void);
void test_nc_halo_cm_param_basic (void);
void test_nc_halo_cm_klypin11_basic (void);
void test_nc_halo_cm_duffy08_basic (void);
void test_nc_halo_cm_bhattacharya13_basic (void);
void test_nc_halo_cm_diemer15_basic (void);
void test_nc_halo_cm_prada12_basic (void);
void test_nc_halo_cm_dutton14_basic (void);
void test_nc_halo_bias_despali_basic (void);
void test_nc_multiplicity_func_bhattacharya_basic (void);
void test_nc_multiplicity_func_bhattacharya_convention (void);
void test_nc_multiplicity_func_bhattacharya_mean_mdef (void);
void test_nc_multiplicity_func_bhattacharya_critical_mdef (void);
void test_nc_multiplicity_func_bhattacharya_virial_mdef (void);
void test_nc_cluster_mass_ascaso_basic (void);
void test_nc_cluster_mass_selection_basic (void);
void test_nc_cluster_photoz_gauss_basic (void);
void test_nc_xcor_basic (void);

void test_nc_galaxy_position_factor_flat_basic (void);
void test_nc_galaxy_redshift_factor_composed_basic (void);
void test_nc_galaxy_redshift_factor_spline_basic (void);
void test_nc_galaxy_redshift_binning_basic (void);
void test_nc_galaxy_redshift_pop_basic (void);
void test_nc_galaxy_shape_factor_basic (void);
void test_nc_galaxy_shape_factor_quad_basic (void);
void test_nc_galaxy_shape_factor_series_lensed_basic (void);
void test_nc_galaxy_shape_factor_fixed_quad_basic (void);
void test_nc_data_cluster_wl_factor_basic (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/data/basic", test_ncm_data_basic);
  g_test_add_func ("/ncm/data_funnel/basic", test_ncm_data_funnel_basic);
  g_test_add_func ("/ncm/data_gaussmix2d/basic", test_ncm_data_gaussmix2d_basic);
  g_test_add_func ("/ncm/data_rosenbrock/basic", test_ncm_data_rosenbrock_basic);
  g_test_add_func ("/ncm/dataset/basic", test_ncm_dataset_basic);
  g_test_add_func ("/ncm/fftlog/basic", test_ncm_fftlog_basic);
  g_test_add_func ("/ncm/sphere/nn", test_ncm_sphere_nn);
  g_test_add_func ("/ncm/timer/basic", test_ncm_timer_basic);
  g_test_add_func ("/ncm/sbessel_ode_solver/basic", test_ncm_sbessel_ode_solver_basic);
  g_test_add_func ("/ncm/sbessel_integrator_gl/basic", test_ncm_sbessel_integrator_gl_basic);
  g_test_add_func ("/ncm/sbessel_integrator_fftl/basic", test_ncm_sbessel_integrator_fftl_basic);
  g_test_add_func ("/ncm/sbessel_integrator_levin/basic", test_ncm_sbessel_integrator_levin_basic);
  g_test_add_func ("/ncm/fftlog_sbessel_j/basic", test_ncm_fftlog_sbessel_j_basic);
  g_test_add_func ("/ncm/fftlog_sbessel_jljm/basic", test_ncm_fftlog_sbessel_jljm_basic);
  g_test_add_func ("/ncm/bootstrap/basic", test_ncm_bootstrap_basic);
  g_test_add_func ("/ncm/stats_vec/basic", test_ncm_stats_vec_basic);
  g_test_add_func ("/ncm/mpi_job/basic", test_ncm_mpi_job_basic);
  g_test_add_func ("/ncm/mpi_job_test/basic", test_ncm_mpi_job_test_basic);
  g_test_add_func ("/ncm/mpi_job_fit/basic", test_ncm_mpi_job_fit_basic);
  g_test_add_func ("/ncm/mpi_job_mcmc/basic", test_ncm_mpi_job_mcmc_basic);
  g_test_add_func ("/ncm/mpi_job_feval/basic", test_ncm_mpi_job_feval_basic);
  g_test_add_func ("/ncm/powspec_spline2d/basic", test_ncm_powspec_spline2d_basic);
  g_test_add_func ("/ncm/pln1d/basic", test_ncm_pln1d_basic);

  g_test_add_func ("/nc/data/cluster_mass_rich/basic", test_nc_data_cluster_mass_rich_basic);
  g_test_add_func ("/nc/data/cluster_mass_rich_count/basic", test_nc_data_cluster_mass_rich_count_basic);
  g_test_add_func ("/nc/data/cluster_wl/basic", test_nc_data_cluster_wl_basic);
  g_test_add_func ("/nc/de_cont/basic", test_nc_de_cont_basic);
  g_test_add_func ("/nc/distance/basic", test_nc_distance_basic);
  g_test_add_func ("/nc/growth_func/basic", test_nc_growth_func_basic);
  g_test_add_func ("/nc/transfer_func/basic", test_nc_transfer_func_basic);
  g_test_add_func ("/nc/galaxy/sd_obs_redshift/basic", test_nc_galaxy_sd_obs_redshift_basic);
  g_test_add_func ("/nc/galaxy/sd_obs_redshift_gauss/basic", test_nc_galaxy_sd_obs_redshift_gauss_basic);
  g_test_add_func ("/nc/galaxy/sd_obs_redshift_spec/basic", test_nc_galaxy_sd_obs_redshift_spec_basic);
  g_test_add_func ("/nc/galaxy/sd_obs_redshift_pz/basic", test_nc_galaxy_sd_obs_redshift_pz_basic);
  g_test_add_func ("/nc/galaxy/sd_position/basic", test_nc_galaxy_sd_position_basic);
  g_test_add_func ("/nc/galaxy/sd_position_flat/basic", test_nc_galaxy_sd_position_flat_basic);
  g_test_add_func ("/nc/galaxy/sd_shape/basic", test_nc_galaxy_sd_shape_basic);
  g_test_add_func ("/nc/galaxy/sd_shape_hsm_gauss_global/basic", test_nc_galaxy_sd_shape_hsm_gauss_global_basic);
  g_test_add_func ("/nc/galaxy/sd_shape_hsm_gauss/basic", test_nc_galaxy_sd_shape_hsm_gauss_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift/basic", test_nc_galaxy_sd_true_redshift_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd_y1_source/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_y1_source_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd_y1_lens/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_y1_lens_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd_y10_source/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_y10_source_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd_y10_lens/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_y10_lens_basic);
  g_test_add_func ("/nc/galaxy/sd_true_redshift_lsst_srd_from_type/basic", test_nc_galaxy_sd_true_redshift_lsst_srd_from_type_basic);
  g_test_add_func ("/nc/galaxy/sd_obs_redshift_gauss_lsst_srd_bins/basic", test_nc_galaxy_sd_obs_redshift_gauss_lsst_srd_bins_basic);
  g_test_add_func ("/nc/galaxy/wl_obs/basic", test_nc_galaxy_wl_obs_basic);
  g_test_add_func ("/nc/halo_position/basic", test_nc_halo_position_basic);
  g_test_add_func ("/nc/hicosmo/qgw/basic", test_nc_hicosmo_qgw_basic);

  g_test_add_func ("/nc/hipert/adiab/basic", test_nc_hipert_adiab_basic);
  g_test_add_func ("/nc/hipert/em/basic", test_nc_hipert_em_basic);
  g_test_add_func ("/nc/hipert/gw/basic", test_nc_hipert_gw_basic);
  g_test_add_func ("/nc/hipert/two_fluids/basic", test_nc_hipert_two_fluids_basic);
  g_test_add_func ("/nc/hiprim/two_fluids/basic", test_nc_hiprim_two_fluids_basic);

  g_test_add_func ("/nc/halo_mass_summary/basic", test_nc_halo_mass_summary_basic);
  g_test_add_func ("/nc/halo_cm_param/basic", test_nc_halo_cm_param_basic);
  g_test_add_func ("/nc/halo_cm_klypin11/basic", test_nc_halo_cm_klypin11_basic);
  g_test_add_func ("/nc/halo_cm_duffy08/basic", test_nc_halo_cm_duffy08_basic);
  g_test_add_func ("/nc/halo_cm_bhattacharya13/basic", test_nc_halo_cm_bhattacharya13_basic);
  g_test_add_func ("/nc/halo_cm_diemer15/basic", test_nc_halo_cm_diemer15_basic);
  g_test_add_func ("/nc/halo_cm_prada12/basic", test_nc_halo_cm_prada12_basic);
  g_test_add_func ("/nc/halo_cm_dutton14/basic", test_nc_halo_cm_dutton14_basic);

  g_test_add_func ("/nc/halo_bias_despali/basic", test_nc_halo_bias_despali_basic);

  g_test_add_func ("/nc/multiplicity_func_bhattacharya/basic", test_nc_multiplicity_func_bhattacharya_basic);
  g_test_add_func ("/nc/multiplicity_func_bhattacharya/convention", test_nc_multiplicity_func_bhattacharya_convention);
  g_test_add_func ("/nc/multiplicity_func_bhattacharya/mean_mdef", test_nc_multiplicity_func_bhattacharya_mean_mdef);
  g_test_add_func ("/nc/multiplicity_func_bhattacharya/critical_mdef", test_nc_multiplicity_func_bhattacharya_critical_mdef);
  g_test_add_func ("/nc/multiplicity_func_bhattacharya/virial_mdef", test_nc_multiplicity_func_bhattacharya_virial_mdef);

  g_test_add_func ("/nc/cluster_mass_ascaso/basic", test_nc_cluster_mass_ascaso_basic);
  g_test_add_func ("/nc/cluster_mass_selection/basic", test_nc_cluster_mass_selection_basic);
  g_test_add_func ("/nc/cluster_photoz_gauss/basic", test_nc_cluster_photoz_gauss_basic);

  g_test_add_func ("/nc/xcor/basic", test_nc_xcor_basic);

  g_test_add_func ("/nc/galaxy/position_factor_flat/basic", test_nc_galaxy_position_factor_flat_basic);
  g_test_add_func ("/nc/galaxy/redshift_factor_composed/basic", test_nc_galaxy_redshift_factor_composed_basic);
  g_test_add_func ("/nc/galaxy/redshift_factor_spline/basic", test_nc_galaxy_redshift_factor_spline_basic);
  g_test_add_func ("/nc/galaxy/redshift_binning/basic", test_nc_galaxy_redshift_binning_basic);
  g_test_add_func ("/nc/galaxy/redshift_pop/basic", test_nc_galaxy_redshift_pop_basic);
  g_test_add_func ("/nc/galaxy/shape_factor/basic", test_nc_galaxy_shape_factor_basic);
  g_test_add_func ("/nc/galaxy/shape_factor_quad/basic", test_nc_galaxy_shape_factor_quad_basic);
  g_test_add_func ("/nc/galaxy/shape_factor_series_lensed/basic", test_nc_galaxy_shape_factor_series_lensed_basic);
  g_test_add_func ("/nc/galaxy/shape_factor_fixed_quad/basic", test_nc_galaxy_shape_factor_fixed_quad_basic);
  g_test_add_func ("/nc/data/cluster_wl_factor/basic", test_nc_data_cluster_wl_factor_basic);

  g_test_run ();
}

void
test_ncm_timer_basic (void)
{
  NcmTimer *nt = ncm_timer_new ();
  NcmTimer *nt2;

  g_assert_true (nt != NULL);
  g_assert_true (NCM_IS_TIMER (nt));

  nt2 = ncm_timer_ref (nt);
  ncm_timer_clear (&nt2);
  g_assert_true (nt2 == NULL);

  g_assert_true (NCM_IS_TIMER (nt));

  NCM_TEST_FREE (ncm_timer_free, nt);
}

void
test_ncm_sbessel_ode_solver_basic (void)
{
  NcmSBesselOdeSolver *solver = ncm_sbessel_ode_solver_new ();
  NcmSBesselOdeSolver *solver2;

  g_assert_true (solver != NULL);
  g_assert_true (NCM_IS_SBESSEL_ODE_SOLVER (solver));

  solver2 = ncm_sbessel_ode_solver_ref (solver);
  ncm_sbessel_ode_solver_clear (&solver2);
  g_assert_true (solver2 == NULL);

  g_assert_true (NCM_IS_SBESSEL_ODE_SOLVER (solver));

  NCM_TEST_FREE (ncm_sbessel_ode_solver_free, solver);
}

void
test_ncm_sbessel_integrator_gl_basic (void)
{
  NcmSBesselIntegratorGL *sbigl = ncm_sbessel_integrator_gl_new (0, 10);
  NcmSBesselIntegratorGL *sbigl2;

  g_assert_true (sbigl != NULL);
  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl));

  sbigl2 = ncm_sbessel_integrator_gl_ref (sbigl);
  ncm_sbessel_integrator_gl_clear (&sbigl2);
  g_assert_true (sbigl2 == NULL);

  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_GL (sbigl));

  NCM_TEST_FREE (ncm_sbessel_integrator_gl_free, sbigl);
}

void
test_ncm_sbessel_integrator_fftl_basic (void)
{
  NcmSBesselIntegratorFFTL *sbilf = ncm_sbessel_integrator_fftl_new (0, 10);
  NcmSBesselIntegratorFFTL *sbilf2;

  g_assert_true (sbilf != NULL);
  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_FFTL (sbilf));

  sbilf2 = ncm_sbessel_integrator_fftl_ref (sbilf);
  ncm_sbessel_integrator_fftl_clear (&sbilf2);
  g_assert_true (sbilf2 == NULL);

  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_FFTL (sbilf));

  NCM_TEST_FREE (ncm_sbessel_integrator_fftl_free, sbilf);
}

void
test_ncm_sbessel_integrator_levin_basic (void)
{
  NcmSBesselIntegratorLevin *sbilv = ncm_sbessel_integrator_levin_new (0, 10);
  NcmSBesselIntegratorLevin *sbilv2;

  g_assert_true (sbilv != NULL);
  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbilv));

  sbilv2 = ncm_sbessel_integrator_levin_ref (sbilv);
  ncm_sbessel_integrator_levin_clear (&sbilv2);
  g_assert_true (sbilv2 == NULL);

  g_assert_true (NCM_IS_SBESSEL_INTEGRATOR_LEVIN (sbilv));

  NCM_TEST_FREE (ncm_sbessel_integrator_levin_free, sbilv);
}

void
test_ncm_fftlog_sbessel_j_basic (void)
{
  NcmFftlogSBesselJ *fftlog_jl = ncm_fftlog_sbessel_j_new (2, 0.0, 0.0, 1.0, 128);
  NcmFftlog *fftlog_jl2;

  g_assert_true (fftlog_jl != NULL);
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_J (fftlog_jl));

  fftlog_jl2 = ncm_fftlog_ref (NCM_FFTLOG (fftlog_jl));
  ncm_fftlog_clear (&fftlog_jl2);
  g_assert_true (fftlog_jl2 == NULL);

  g_assert_true (NCM_IS_FFTLOG_SBESSEL_J (fftlog_jl));

  NCM_TEST_FREE (ncm_fftlog_free, NCM_FFTLOG (fftlog_jl));
}

void
test_ncm_fftlog_sbessel_jljm_basic (void)
{
  NcmFftlogSBesselJLJM *fftlog_jljm = ncm_fftlog_sbessel_jljm_new (2, 0, 0.0, 0.0, 0.0, 1.0, 128);
  NcmFftlog *fftlog_jljm2;

  g_assert_true (fftlog_jljm != NULL);
  g_assert_true (NCM_IS_FFTLOG_SBESSEL_JLJM (fftlog_jljm));

  fftlog_jljm2 = ncm_fftlog_ref (NCM_FFTLOG (fftlog_jljm));
  ncm_fftlog_clear (&fftlog_jljm2);
  g_assert_true (fftlog_jljm2 == NULL);

  g_assert_true (NCM_IS_FFTLOG_SBESSEL_JLJM (fftlog_jljm));

  NCM_TEST_FREE (ncm_fftlog_free, NCM_FFTLOG (fftlog_jljm));
}

void
test_ncm_bootstrap_basic (void)
{
  NcmBootstrap *bstrap = ncm_bootstrap_new ();
  NcmBootstrap *bstrap2;

  g_assert_true (bstrap != NULL);
  g_assert_true (NCM_IS_BOOTSTRAP (bstrap));

  bstrap2 = ncm_bootstrap_ref (bstrap);
  ncm_bootstrap_clear (&bstrap2);
  g_assert_true (bstrap2 == NULL);

  g_assert_true (NCM_IS_BOOTSTRAP (bstrap));

  NCM_TEST_FREE (ncm_bootstrap_free, bstrap);
}

void
test_ncm_stats_vec_basic (void)
{
  NcmStatsVec *svec = ncm_stats_vec_new (3, NCM_STATS_VEC_MEAN, FALSE);
  NcmStatsVec *svec2;

  g_assert_true (svec != NULL);
  g_assert_true (NCM_IS_STATS_VEC (svec));

  svec2 = ncm_stats_vec_ref (svec);
  ncm_stats_vec_clear (&svec2);
  g_assert_true (svec2 == NULL);

  g_assert_true (NCM_IS_STATS_VEC (svec));

  NCM_TEST_FREE (ncm_stats_vec_free, svec);
}

void
test_ncm_data_basic (void)
{
  NcmData *data = NCM_DATA (ncm_data_funnel_new ());
  NcmData *data2;

  g_assert_true (data != NULL);
  g_assert_true (NCM_IS_DATA (data));

  data2 = ncm_data_ref (data);
  ncm_data_clear (&data2);
  g_assert_true (data2 == NULL);

  g_assert_true (NCM_IS_DATA (data));

  NCM_TEST_FREE (ncm_data_free, data);
}

void
test_ncm_data_funnel_basic (void)
{
  NcmDataFunnel *funnel = ncm_data_funnel_new ();
  NcmDataFunnel *funnel2;

  g_assert_true (funnel != NULL);
  g_assert_true (NCM_IS_DATA_FUNNEL (funnel));

  funnel2 = ncm_data_funnel_ref (funnel);
  ncm_data_funnel_clear (&funnel2);
  g_assert_true (funnel2 == NULL);

  g_assert_true (NCM_IS_DATA_FUNNEL (funnel));

  NCM_TEST_FREE (ncm_data_funnel_free, funnel);
}

void
test_ncm_data_gaussmix2d_basic (void)
{
  NcmDataGaussMix2D *gm2d = ncm_data_gaussmix2d_new ();
  NcmDataGaussMix2D *gm2d2;

  g_assert_true (gm2d != NULL);
  g_assert_true (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  gm2d2 = ncm_data_gaussmix2d_ref (gm2d);
  ncm_data_gaussmix2d_clear (&gm2d2);
  g_assert_true (gm2d2 == NULL);

  g_assert_true (NCM_IS_DATA_GAUSSMIX2D (gm2d));

  NCM_TEST_FREE (ncm_data_gaussmix2d_free, gm2d);
}

void
test_ncm_data_rosenbrock_basic (void)
{
  NcmDataRosenbrock *rosenbrock = ncm_data_rosenbrock_new ();
  NcmDataRosenbrock *rosenbrock2;

  g_assert_true (rosenbrock != NULL);
  g_assert_true (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  rosenbrock2 = ncm_data_rosenbrock_ref (rosenbrock);
  ncm_data_rosenbrock_clear (&rosenbrock2);
  g_assert_true (rosenbrock2 == NULL);

  g_assert_true (NCM_IS_DATA_ROSENBROCK (rosenbrock));

  NCM_TEST_FREE (ncm_data_rosenbrock_free, rosenbrock);
}

void
test_ncm_dataset_basic (void)
{
  NcmDataset *dset = ncm_dataset_new ();
  NcmDataset *dset2;

  g_assert_true (dset != NULL);
  g_assert_true (NCM_IS_DATASET (dset));

  dset2 = ncm_dataset_ref (dset);
  ncm_dataset_clear (&dset2);
  g_assert_true (dset2 == NULL);

  g_assert_true (NCM_IS_DATASET (dset));

  NCM_TEST_FREE (ncm_dataset_free, dset);
}

void
test_ncm_fftlog_basic (void)
{
  NcmFftlog *fftlog = NCM_FFTLOG (ncm_fftlog_gausswin2_new (0.0, 0.0, 20.0, 1000));
  NcmFftlog *fftlog2;

  g_assert_true (fftlog != NULL);
  g_assert_true (NCM_IS_FFTLOG (fftlog));

  fftlog2 = ncm_fftlog_ref (fftlog);
  ncm_fftlog_clear (&fftlog2);
  g_assert_true (fftlog2 == NULL);

  g_assert_true (NCM_IS_FFTLOG (fftlog));

  NCM_TEST_FREE (ncm_fftlog_free, fftlog);
}

void
test_ncm_sphere_nn (void)
{
  NcmSphereNN *snn = ncm_sphere_nn_new ();
  NcmSphereNN *snn2;

  g_assert_true (snn != NULL);
  g_assert_true (NCM_IS_SPHERE_NN (snn));

  snn2 = ncm_sphere_nn_ref (snn);
  ncm_sphere_nn_clear (&snn2);
  g_assert_true (snn2 == NULL);

  g_assert_true (NCM_IS_SPHERE_NN (snn));

  NCM_TEST_FREE (ncm_sphere_nn_free, snn);
}

void
test_ncm_mpi_job_basic (void)
{
  NcmMPIJob *mpi_job = NCM_MPI_JOB (ncm_mpi_job_test_new ());
  NcmMPIJob *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB (mpi_job));

  mpi_job2 = ncm_mpi_job_ref (mpi_job);
  ncm_mpi_job_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_free, mpi_job);
}

void
test_ncm_mpi_job_test_basic (void)
{
  NcmMPIJobTest *mpi_job = ncm_mpi_job_test_new ();
  NcmMPIJobTest *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_TEST (mpi_job));

  mpi_job2 = ncm_mpi_job_test_ref (mpi_job);
  ncm_mpi_job_test_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_TEST (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_test_free, mpi_job);
}

void
test_ncm_mpi_job_fit_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobFit *mpi_job     = ncm_mpi_job_fit_new (fit, NULL);
  NcmMPIJobFit *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_FIT (mpi_job));

  mpi_job2 = ncm_mpi_job_fit_ref (mpi_job);
  ncm_mpi_job_fit_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_FIT (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_fit_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_mpi_job_mcmc_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobMCMC *mpi_job    = ncm_mpi_job_mcmc_new (fit, NULL);
  NcmMPIJobMCMC *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_MCMC (mpi_job));

  mpi_job2 = ncm_mpi_job_mcmc_ref (mpi_job);
  ncm_mpi_job_mcmc_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_MCMC (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_mcmc_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_mpi_job_feval_basic (void)
{
  NcmMSet *mset             = ncm_mset_empty_new ();
  NcmData *data             = NCM_DATA (ncm_data_funnel_new ());
  NcmDataset *dset          = ncm_dataset_new_list (data, NULL);
  NcmLikelihood *likelihood = ncm_likelihood_new (dset);
  NcmFit *fit               = ncm_fit_factory (NCM_FIT_TYPE_GSL_MMS, NULL, likelihood, mset, NCM_FIT_GRAD_NUMDIFF_FORWARD);
  NcmMPIJobFEval *mpi_job   = ncm_mpi_job_feval_new (fit, NULL);
  NcmMPIJobFEval *mpi_job2;

  g_assert_true (mpi_job != NULL);
  g_assert_true (NCM_IS_MPI_JOB_FEVAL (mpi_job));

  mpi_job2 = ncm_mpi_job_feval_ref (mpi_job);
  ncm_mpi_job_feval_clear (&mpi_job2);
  g_assert_true (mpi_job2 == NULL);

  g_assert_true (NCM_IS_MPI_JOB_FEVAL (mpi_job));

  NCM_TEST_FREE (ncm_mpi_job_feval_free, mpi_job);

  NCM_TEST_FREE (ncm_fit_free, fit);
  NCM_TEST_FREE (ncm_likelihood_free, likelihood);
  NCM_TEST_FREE (ncm_dataset_free, dset);
  NCM_TEST_FREE (ncm_data_free, data);
  NCM_TEST_FREE (ncm_mset_free, mset);
}

void
test_ncm_powspec_spline2d_basic (void)
{
  gdouble x[6]      = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  gdouble y[7]      = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  gdouble z[42]     = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, };
  NcmVector *xv     = ncm_vector_new_data_static (x, 6, 1);
  NcmVector *yv     = ncm_vector_new_data_static (y, 7, 1);
  NcmMatrix *zm     = ncm_matrix_new_data_static (z, 7, 6);
  NcmSpline *sc     = NCM_SPLINE (ncm_spline_cubic_notaknot_new ());
  NcmSpline2d *sb2d = ncm_spline2d_bicubic_new (sc);

  ncm_spline2d_set (NCM_SPLINE2D (sb2d), xv, yv, zm, FALSE);

  {
    NcmPowspecSpline2d *ps_s2d = ncm_powspec_spline2d_new (NCM_SPLINE2D (sb2d));
    NcmPowspecSpline2d *ps_s2d2;

    g_assert_true (ps_s2d != NULL);
    g_assert_true (NCM_IS_POWSPEC_SPLINE2D (ps_s2d));

    ps_s2d2 = ncm_powspec_spline2d_ref (ps_s2d);
    ncm_powspec_spline2d_clear (&ps_s2d2);
    g_assert_true (ps_s2d2 == NULL);

    g_assert_true (NCM_IS_POWSPEC_SPLINE2D (ps_s2d));

    ncm_vector_free (xv);
    ncm_vector_free (yv);
    ncm_matrix_free (zm);
    ncm_spline_free (NCM_SPLINE (sc));

    NCM_TEST_FREE (ncm_powspec_spline2d_free, ps_s2d);
    NCM_TEST_FREE (ncm_spline2d_free, NCM_SPLINE2D (sb2d));
  }
}

void
test_ncm_pln1d_basic (void)
{
  NcmPLN1D *pln1d = ncm_pln1d_new (60);
  NcmPLN1D *pln1d2;

  g_assert_true (pln1d != NULL);
  g_assert_true (NCM_IS_PLN1D (pln1d));

  pln1d2 = ncm_pln1d_ref (pln1d);
  ncm_pln1d_clear (&pln1d2);
  g_assert_true (pln1d2 == NULL);

  g_assert_true (NCM_IS_PLN1D (pln1d));

  NCM_TEST_FREE (ncm_pln1d_free, pln1d);
}

void
test_nc_galaxy_sd_obs_redshift_basic (void)
{
  NcGalaxySDTrueRedshift *gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshift *gsdor  = NC_GALAXY_SD_OBS_REDSHIFT (nc_galaxy_sd_obs_redshift_gauss_new (gsdtr, 0.1, 1.0));
  NcGalaxySDObsRedshift *gsdor2;

  g_assert_true (gsdor != NULL);
  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  gsdor2 = nc_galaxy_sd_obs_redshift_ref (gsdor);
  nc_galaxy_sd_obs_redshift_clear (&gsdor2);
  g_assert_true (gsdor2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT (gsdor));

  nc_galaxy_sd_true_redshift_clear (&gsdtr);
  NCM_TEST_FREE (nc_galaxy_sd_obs_redshift_free, gsdor);
}

void
test_nc_galaxy_sd_obs_redshift_gauss_basic (void)
{
  NcGalaxySDTrueRedshift *gsdtr      = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshiftGauss *gsdorg = nc_galaxy_sd_obs_redshift_gauss_new (gsdtr, 0.1, 1.0);
  NcGalaxySDObsRedshiftGauss *gsdorg2;

  g_assert_true (gsdorg != NULL);
  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorg));

  gsdorg2 = nc_galaxy_sd_obs_redshift_gauss_ref (gsdorg);
  nc_galaxy_sd_obs_redshift_gauss_clear (&gsdorg2);
  g_assert_true (gsdorg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (gsdorg));

  nc_galaxy_sd_true_redshift_clear (&gsdtr);
  NCM_TEST_FREE (nc_galaxy_sd_obs_redshift_gauss_free, gsdorg);
}

void
test_nc_galaxy_sd_obs_redshift_spec_basic (void)
{
  NcGalaxySDTrueRedshift *gsdtr     = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDObsRedshiftSpec *gsdors = nc_galaxy_sd_obs_redshift_spec_new (gsdtr, 0.1, 1.0);
  NcGalaxySDObsRedshiftSpec *gsdors2;

  g_assert_true (gsdors != NULL);
  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdors));

  gsdors2 = nc_galaxy_sd_obs_redshift_spec_ref (gsdors);
  nc_galaxy_sd_obs_redshift_spec_clear (&gsdors2);
  g_assert_true (gsdors2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_SPEC (gsdors));

  nc_galaxy_sd_true_redshift_clear (&gsdtr);
  NCM_TEST_FREE (nc_galaxy_sd_obs_redshift_spec_free, gsdors);
}

void
test_nc_galaxy_sd_obs_redshift_pz_basic (void)
{
  NcGalaxySDObsRedshiftPz *gsdorpz = nc_galaxy_sd_obs_redshift_pz_new ();
  NcGalaxySDObsRedshiftPz *gsdorpz2;

  g_assert_true (gsdorpz != NULL);
  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (gsdorpz));

  gsdorpz2 = nc_galaxy_sd_obs_redshift_pz_ref (gsdorpz);
  nc_galaxy_sd_obs_redshift_pz_clear (&gsdorpz2);
  g_assert_true (gsdorpz2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_PZ (gsdorpz));

  NCM_TEST_FREE (nc_galaxy_sd_obs_redshift_pz_free, gsdorpz);
}

void
test_nc_halo_position_basic (void)
{
  NcHaloPosition *hp;
  NcHaloPosition *hp2;
  NcDistance *dist = nc_distance_new (1100.0);
  NcHICosmo *cosmo = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());

  nc_distance_prepare (dist, cosmo);

  hp = nc_halo_position_new (dist);

  g_assert_true (hp != NULL);
  g_assert_true (NC_IS_HALO_POSITION (hp));

  hp2 = nc_halo_position_ref (hp);
  nc_halo_position_clear (&hp2);
  g_assert_true (hp2 == NULL);

  g_assert_true (NC_IS_HALO_POSITION (hp));

  nc_distance_clear (&dist);
  nc_hicosmo_clear (&cosmo);
  NCM_TEST_FREE (nc_halo_position_free, hp);
}

void
test_nc_hicosmo_qgw_basic (void)
{
  NcHICosmoQGW *qgw = nc_hicosmo_qgw_new ();

  g_assert_true (qgw != NULL);
  g_assert_true (NC_IS_HICOSMO_QGW (qgw));

  NCM_TEST_FREE (nc_hicosmo_free, NC_HICOSMO (qgw));
}

void
test_nc_data_cluster_mass_rich_basic (void)
{
  NcDataClusterMassRich *dcmr = nc_data_cluster_mass_rich_new ();
  NcDataClusterMassRich *dcmr2;

  g_assert_true (dcmr != NULL);
  g_assert_true (NC_IS_DATA_CLUSTER_MASS_RICH (dcmr));

  dcmr2 = nc_data_cluster_mass_rich_ref (dcmr);
  nc_data_cluster_mass_rich_clear (&dcmr2);
  g_assert_true (dcmr2 == NULL);

  g_assert_true (NC_IS_DATA_CLUSTER_MASS_RICH (dcmr));

  NCM_TEST_FREE (nc_data_cluster_mass_rich_free, dcmr);
}

void
test_nc_data_cluster_mass_rich_count_basic (void)
{
  NcDataClusterMassRichCount *dmrc = nc_data_cluster_mass_rich_count_new ();
  NcDataClusterMassRichCount *dmrc2;

  g_assert_true (dmrc != NULL);
  g_assert_true (NC_IS_DATA_CLUSTER_MASS_RICH_COUNT (dmrc));

  dmrc2 = nc_data_cluster_mass_rich_count_ref (dmrc);
  nc_data_cluster_mass_rich_count_clear (&dmrc2);
  g_assert_true (dmrc2 == NULL);

  g_assert_true (NC_IS_DATA_CLUSTER_MASS_RICH_COUNT (dmrc));

  NCM_TEST_FREE (nc_data_cluster_mass_rich_count_free, dmrc);
}

void
test_nc_data_cluster_wl_basic (void)
{
  NcDataClusterWL *dcwl = nc_data_cluster_wl_new ();
  NcDataClusterWL *dcwl2;

  g_assert_true (dcwl != NULL);
  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));

  dcwl2 = nc_data_cluster_wl_ref (dcwl);
  nc_data_cluster_wl_clear (&dcwl2);
  g_assert_true (dcwl2 == NULL);

  g_assert_true (NC_IS_DATA_CLUSTER_WL (dcwl));

  NCM_TEST_FREE (nc_data_cluster_wl_free, dcwl);
}

void
test_nc_distance_basic (void)
{
  NcDistance *dist = nc_distance_new (6.0);
  NcDistance *dist2;

  g_assert_true (dist != NULL);
  g_assert_true (NC_IS_DISTANCE (dist));

  dist2 = nc_distance_ref (dist);
  nc_distance_clear (&dist2);
  g_assert_true (dist2 == NULL);

  g_assert_true (NC_IS_DISTANCE (dist));

  NCM_TEST_FREE (nc_distance_free, dist);
}

void
test_nc_growth_func_basic (void)
{
  NcGrowthFunc *gf = nc_growth_func_new ();
  NcGrowthFunc *gf2;

  g_assert_true (gf != NULL);
  g_assert_true (NC_IS_GROWTH_FUNC (gf));

  gf2 = nc_growth_func_ref (gf);
  nc_growth_func_clear (&gf2);
  g_assert_true (gf2 == NULL);

  g_assert_true (NC_IS_GROWTH_FUNC (gf));

  NCM_TEST_FREE (nc_growth_func_free, gf);
}

void
test_nc_transfer_func_basic (void)
{
  NcTransferFunc *tf = NC_TRANSFER_FUNC (nc_transfer_func_eh_new ());
  NcTransferFunc *tf2;

  g_assert_true (tf != NULL);
  g_assert_true (NC_IS_TRANSFER_FUNC (tf));

  tf2 = nc_transfer_func_ref (tf);
  nc_transfer_func_clear (&tf2);
  g_assert_true (tf2 == NULL);

  g_assert_true (NC_IS_TRANSFER_FUNC (tf));

  NCM_TEST_FREE (nc_transfer_func_free, tf);
}

void
test_nc_de_cont_basic (void)
{
  NcDECont *dec = nc_de_cont_new (0.3, 0.7, 0.01, 0.01);
  NcDECont *dec2;

  g_assert_true (dec != NULL);
  g_assert_true (NC_IS_DE_CONT (dec));

  dec2 = nc_de_cont_ref (dec);
  nc_de_cont_clear (&dec2);
  g_assert_true (dec2 == NULL);

  g_assert_true (NC_IS_DE_CONT (dec));

  nc_de_cont_set_k (dec, 0.1);
  g_assert_cmpfloat (nc_de_cont_get_k (dec), ==, 0.1);

  NCM_TEST_FREE (nc_de_cont_free, dec);
}

void
test_nc_galaxy_sd_position_basic (void)
{
  NcGalaxySDPosition *sdpos = NC_GALAXY_SD_POSITION (nc_galaxy_sd_position_flat_new (0.0, 2.0, -0.5, 0.5));
  NcGalaxySDPosition *sdpos2;

  g_assert_true (sdpos != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION (sdpos));

  sdpos2 = nc_galaxy_sd_position_ref (sdpos);
  nc_galaxy_sd_position_clear (&sdpos2);
  g_assert_true (sdpos2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION (sdpos));

  NCM_TEST_FREE (nc_galaxy_sd_position_free, sdpos);
}

void
test_nc_galaxy_sd_position_flat_basic (void)
{
  NcGalaxySDPositionFlat *gsdpf = nc_galaxy_sd_position_flat_new (0.0, 2.0, -0.5, 0.5);
  NcGalaxySDPositionFlat *gsdpf2;

  g_assert_true (gsdpf != NULL);
  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpf));

  gsdpf2 = nc_galaxy_sd_position_flat_ref (gsdpf);
  nc_galaxy_sd_position_flat_clear (&gsdpf2);
  g_assert_true (gsdpf2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_POSITION_FLAT (gsdpf));

  NCM_TEST_FREE (nc_galaxy_sd_position_flat_free, gsdpf);
}

void
test_nc_galaxy_sd_shape_basic (void)
{
  NcGalaxySDShape *gsds = NC_GALAXY_SD_SHAPE (nc_galaxy_sd_shape_hsm_gauss_global_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));
  NcGalaxySDShape *gsds2;

  g_assert_true (gsds != NULL);
  g_assert_true (NC_IS_GALAXY_SD_SHAPE (gsds));

  gsds2 = nc_galaxy_sd_shape_ref (gsds);
  nc_galaxy_sd_shape_clear (&gsds2);
  g_assert_true (gsds2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE (gsds));

  NCM_TEST_FREE (nc_galaxy_sd_shape_free, gsds);
}

void
test_nc_galaxy_sd_shape_hsm_gauss_global_basic (void)
{
  NcGalaxySDShapeHSMGaussGlobal *gsdsg = nc_galaxy_sd_shape_hsm_gauss_global_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxySDShapeHSMGaussGlobal *gsdsg2;

  g_assert_true (gsdsg != NULL);
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (gsdsg));

  gsdsg2 = nc_galaxy_sd_shape_hsm_gauss_global_ref (gsdsg);
  nc_galaxy_sd_shape_hsm_gauss_global_clear (&gsdsg2);
  g_assert_true (gsdsg2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS_GLOBAL (gsdsg));

  NCM_TEST_FREE (nc_galaxy_sd_shape_hsm_gauss_global_free, gsdsg);
}

void
test_nc_galaxy_sd_shape_hsm_gauss_basic (void)
{
  NcGalaxySDShapeHSMGauss *gsdsgh = nc_galaxy_sd_shape_hsm_gauss_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxySDShapeHSMGauss *gsdsgh2;

  g_assert_true (gsdsgh != NULL);
  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS (gsdsgh));

  gsdsgh2 = nc_galaxy_sd_shape_hsm_gauss_ref (gsdsgh);
  nc_galaxy_sd_shape_hsm_gauss_clear (&gsdsgh2);
  g_assert_true (gsdsgh2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_SHAPE_HSM_GAUSS (gsdsgh));

  NCM_TEST_FREE (nc_galaxy_sd_shape_hsm_gauss_free, gsdsgh);
}

void
test_nc_galaxy_sd_true_redshift_basic (void)
{
  NcGalaxySDTrueRedshift *gsdtr = NC_GALAXY_SD_TRUE_REDSHIFT (nc_galaxy_sd_true_redshift_lsst_srd_new ());
  NcGalaxySDTrueRedshift *gsdtr2;

  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT (gsdtr));

  gsdtr2 = nc_galaxy_sd_true_redshift_ref (gsdtr);
  nc_galaxy_sd_true_redshift_clear (&gsdtr2);
  g_assert_true (gsdtr2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT (gsdtr));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_free, gsdtr);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst = nc_galaxy_sd_true_redshift_lsst_srd_new ();
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtrlsst2;

  g_assert_true (gsdtrlsst != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtrlsst));

  gsdtrlsst2 = nc_galaxy_sd_true_redshift_lsst_srd_ref (gsdtrlsst);
  nc_galaxy_sd_true_redshift_lsst_srd_clear (&gsdtrlsst2);
  g_assert_true (gsdtrlsst2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtrlsst));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtrlsst);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_y1_source_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = nc_galaxy_sd_true_redshift_lsst_srd_new_y1_source ();
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr2;

  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  gsdtr2 = nc_galaxy_sd_true_redshift_lsst_srd_ref (gsdtr);
  nc_galaxy_sd_true_redshift_lsst_srd_clear (&gsdtr2);
  g_assert_true (gsdtr2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_y1_lens_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = nc_galaxy_sd_true_redshift_lsst_srd_new_y1_lens ();
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr2;

  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  gsdtr2 = nc_galaxy_sd_true_redshift_lsst_srd_ref (gsdtr);
  nc_galaxy_sd_true_redshift_lsst_srd_clear (&gsdtr2);
  g_assert_true (gsdtr2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_y10_source_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = nc_galaxy_sd_true_redshift_lsst_srd_new_y10_source ();
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr2;

  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  gsdtr2 = nc_galaxy_sd_true_redshift_lsst_srd_ref (gsdtr);
  nc_galaxy_sd_true_redshift_lsst_srd_clear (&gsdtr2);
  g_assert_true (gsdtr2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_y10_lens_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = nc_galaxy_sd_true_redshift_lsst_srd_new_y10_lens ();
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr2;

  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  gsdtr2 = nc_galaxy_sd_true_redshift_lsst_srd_ref (gsdtr);
  nc_galaxy_sd_true_redshift_lsst_srd_clear (&gsdtr2);
  g_assert_true (gsdtr2 == NULL);

  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr);
}

void
test_nc_galaxy_sd_true_redshift_lsst_srd_from_type_basic (void)
{
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr_y1_source  = nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE);
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr_y1_lens    = nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS);
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr_y10_source = nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE);
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr_y10_lens   = nc_galaxy_sd_true_redshift_lsst_srd_new_from_type (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS);

  g_assert_true (gsdtr_y1_source != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr_y1_source));

  g_assert_true (gsdtr_y1_lens != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr_y1_lens));

  g_assert_true (gsdtr_y10_source != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr_y10_source));

  g_assert_true (gsdtr_y10_lens != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr_y10_lens));

  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr_y1_source);
  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr_y1_lens);
  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr_y10_source);
  NCM_TEST_FREE (nc_galaxy_sd_true_redshift_lsst_srd_free, gsdtr_y10_lens);
}

void
test_nc_galaxy_sd_obs_redshift_gauss_lsst_srd_bins_basic (void)
{
  /* Test Y1 source bins (5 bins) */
  GPtrArray *y1_source_bins = nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE, NULL);

  g_assert_true (y1_source_bins != NULL);
  g_assert_cmpuint (y1_source_bins->len, ==, 5);

  for (guint i = 0; i < y1_source_bins->len; i++)
  {
    NcGalaxySDObsRedshiftGauss *bin = g_ptr_array_index (y1_source_bins, i);

    g_assert_true (bin != NULL);
    g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (bin));
  }

  g_ptr_array_unref (y1_source_bins);

  /* Test Y1 lens bins (5 bins) */
  GPtrArray *y1_lens_bins = nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_LENS, NULL);

  g_assert_true (y1_lens_bins != NULL);
  g_assert_cmpuint (y1_lens_bins->len, ==, 5);

  for (guint i = 0; i < y1_lens_bins->len; i++)
  {
    NcGalaxySDObsRedshiftGauss *bin = g_ptr_array_index (y1_lens_bins, i);

    g_assert_true (bin != NULL);
    g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (bin));
  }

  g_ptr_array_unref (y1_lens_bins);

  /* Test Y10 source bins (5 bins) */
  GPtrArray *y10_source_bins = nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_SOURCE, NULL);

  g_assert_true (y10_source_bins != NULL);
  g_assert_cmpuint (y10_source_bins->len, ==, 5);

  for (guint i = 0; i < y10_source_bins->len; i++)
  {
    NcGalaxySDObsRedshiftGauss *bin = g_ptr_array_index (y10_source_bins, i);

    g_assert_true (bin != NULL);
    g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (bin));
  }

  g_ptr_array_unref (y10_source_bins);

  /* Test Y10 lens bins (10 bins) */
  GPtrArray *y10_lens_bins = nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y10_LENS, NULL);

  g_assert_true (y10_lens_bins != NULL);
  g_assert_cmpuint (y10_lens_bins->len, ==, 10);

  for (guint i = 0; i < y10_lens_bins->len; i++)
  {
    NcGalaxySDObsRedshiftGauss *bin = g_ptr_array_index (y10_lens_bins, i);

    g_assert_true (bin != NULL);
    g_assert_true (NC_IS_GALAXY_SD_OBS_REDSHIFT_GAUSS (bin));
  }

  g_ptr_array_unref (y10_lens_bins);

  /* Test retrieving the true redshift distribution */
  NcGalaxySDTrueRedshiftLSSTSRD *gsdtr = NULL;
  GPtrArray *bins_with_gsdtr           = nc_galaxy_sd_obs_redshift_gauss_new_lsst_srd_bins (NC_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD_Y1_SOURCE, &gsdtr);

  g_assert_true (bins_with_gsdtr != NULL);
  g_assert_true (gsdtr != NULL);
  g_assert_true (NC_IS_GALAXY_SD_TRUE_REDSHIFT_LSST_SRD (gsdtr));

  g_ptr_array_unref (bins_with_gsdtr);
  nc_galaxy_sd_true_redshift_lsst_srd_free (gsdtr);
}

void
test_nc_galaxy_wl_obs_basic (void)
{
  GStrv names        = g_strsplit ("e1 e2 ra dec z sz", " ", -1);
  NcGalaxyWLObs *gwl = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, NC_GALAXY_WL_OBS_COORD_CELESTIAL, 100, names);
  NcGalaxyWLObs *gwl2;

  g_assert_true (gwl != NULL);
  g_assert_true (NC_IS_GALAXY_WL_OBS (gwl));

  gwl2 = nc_galaxy_wl_obs_ref (gwl);
  nc_galaxy_wl_obs_clear (&gwl2);
  g_assert_true (gwl2 == NULL);

  g_assert_true (NC_IS_GALAXY_WL_OBS (gwl));

  NCM_TEST_FREE (nc_galaxy_wl_obs_free, gwl);
  g_strfreev (names);
}

void
test_nc_hipert_adiab_basic (void)
{
  NcHIPertAdiab *adiab = nc_hipert_adiab_new ();
  NcHIPertAdiab *adiab2;

  g_assert_true (adiab != NULL);
  g_assert_true (NC_IS_HIPERT_ADIAB (adiab));

  adiab2 = nc_hipert_adiab_ref (adiab);
  nc_hipert_adiab_clear (&adiab2);
  g_assert_true (adiab2 == NULL);

  g_assert_true (NC_IS_HIPERT_ADIAB (adiab));

  nc_hipert_adiab_set_k (adiab, 0.1);
  g_assert_cmpfloat (nc_hipert_adiab_get_k (adiab), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_adiab_free, adiab);
}

void
test_nc_hipert_em_basic (void)
{
  NcHIPertEM *em = nc_hipert_em_new ();
  NcHIPertEM *em2;

  g_assert_true (em != NULL);
  g_assert_true (NC_IS_HIPERT_EM (em));

  em2 = nc_hipert_em_ref (em);
  nc_hipert_em_clear (&em2);
  g_assert_true (em2 == NULL);

  g_assert_true (NC_IS_HIPERT_EM (em));

  nc_hipert_em_set_k (em, 0.1);
  g_assert_cmpfloat (nc_hipert_em_get_k (em), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_em_free, em);
}

void
test_nc_hipert_gw_basic (void)
{
  NcHIPertGW *gw = nc_hipert_gw_new ();
  NcHIPertGW *gw2;

  g_assert_true (gw != NULL);
  g_assert_true (NC_IS_HIPERT_GW (gw));

  gw2 = nc_hipert_gw_ref (gw);
  nc_hipert_gw_clear (&gw2);
  g_assert_true (gw2 == NULL);

  g_assert_true (NC_IS_HIPERT_GW (gw));

  nc_hipert_gw_set_k (gw, 0.1);
  g_assert_cmpfloat (nc_hipert_gw_get_k (gw), ==, 0.1);

  NCM_TEST_FREE (nc_hipert_gw_free, gw);
}

void
test_nc_hipert_two_fluids_basic (void)
{
  NcHIPertTwoFluids *tf = nc_hipert_two_fluids_new ();
  NcHIPertTwoFluids *tf2;

  g_assert_true (tf != NULL);
  g_assert_true (NC_IS_HIPERT_TWO_FLUIDS (tf));

  tf2 = nc_hipert_two_fluids_ref (tf);
  nc_hipert_two_fluids_clear (&tf2);
  g_assert_true (tf2 == NULL);

  g_assert_true (NC_IS_HIPERT_TWO_FLUIDS (tf));

  NCM_TEST_FREE (nc_hipert_two_fluids_free, tf);
}

void
test_nc_hiprim_two_fluids_basic (void)
{
  NcHIPrimTwoFluids *tf = nc_hiprim_two_fluids_new ();
  NcHIPrimTwoFluids *tf2;

  g_assert_true (tf != NULL);
  g_assert_true (NC_IS_HIPRIM_TWO_FLUIDS (tf));

  tf2 = nc_hiprim_two_fluids_ref (tf);
  nc_hiprim_two_fluids_clear (&tf2);
  g_assert_true (tf2 == NULL);

  g_assert_true (NC_IS_HIPRIM_TWO_FLUIDS (tf));

  NCM_TEST_FREE (nc_hiprim_two_fluids_free, tf);
}

void
test_nc_halo_mass_summary_basic (void)
{
  NcHaloMassSummary *hms = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0));
  NcHaloMassSummary *hms2;

  g_assert_true (hms != NULL);
  g_assert_true (NC_IS_HALO_MASS_SUMMARY (hms));

  hms2 = nc_halo_mass_summary_ref (hms);
  nc_halo_mass_summary_clear (&hms2);
  g_assert_true (hms2 == NULL);

  g_assert_true (NC_IS_HALO_MASS_SUMMARY (hms));

  NCM_TEST_FREE (nc_halo_mass_summary_free, hms);
}

void
test_nc_halo_cm_param_basic (void)
{
  NcHaloCMParam *hmp = nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0);
  NcHaloCMParam *hmp2;

  g_assert_true (hmp != NULL);
  g_assert_true (NC_IS_HALO_CM_PARAM (hmp));

  hmp2 = nc_halo_cm_param_ref (hmp);
  nc_halo_cm_param_clear (&hmp2);
  g_assert_true (hmp2 == NULL);

  g_assert_true (NC_IS_HALO_CM_PARAM (hmp));

  NCM_TEST_FREE (nc_halo_cm_param_free, hmp);
}

void
test_nc_halo_cm_klypin11_basic (void)
{
  NcHaloCMKlypin11 *hcmk = nc_halo_cm_klypin11_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0);
  NcHaloCMKlypin11 *hcmk2;

  g_assert_true (hcmk != NULL);
  g_assert_true (NC_IS_HALO_CM_KLYPIN11 (hcmk));

  hcmk2 = nc_halo_cm_klypin11_ref (hcmk);
  nc_halo_cm_klypin11_clear (&hcmk2);
  g_assert_true (hcmk2 == NULL);

  g_assert_true (NC_IS_HALO_CM_KLYPIN11 (hcmk));

  NCM_TEST_FREE (nc_halo_cm_klypin11_free, hcmk);
}

void
test_nc_halo_cm_duffy08_basic (void)
{
  NcHaloCMDuffy08 *hcmd = nc_halo_cm_duffy08_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0);
  NcHaloCMDuffy08 *hcmd2;

  g_assert_true (hcmd != NULL);
  g_assert_true (NC_IS_HALO_CM_DUFFY08 (hcmd));

  hcmd2 = nc_halo_cm_duffy08_ref (hcmd);
  nc_halo_cm_duffy08_clear (&hcmd2);
  g_assert_true (hcmd2 == NULL);

  g_assert_true (NC_IS_HALO_CM_DUFFY08 (hcmd));

  NCM_TEST_FREE (nc_halo_cm_duffy08_free, hcmd);
}

void
test_nc_halo_cm_bhattacharya13_basic (void)
{
  NcDistance *dist        = nc_distance_new (1100.0);
  NcTransferFunc *tf      = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml      = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf   = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mf  = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new ());
  NcHaloMassFunction *hmf = nc_halo_mass_function_new (dist, psf, mf);
  NcGrowthFunc *gf        = nc_growth_func_new ();
  NcHaloCMBhattacharya13 *hcmb;
  NcHaloCMBhattacharya13 *hcmb2;

  nc_multiplicity_func_set_mdef (mf, NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL);
  nc_multiplicity_func_set_Delta (mf, 200.0);

  hcmb = nc_halo_cm_bhattacharya13_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0, hmf, gf);

  g_assert_true (hcmb != NULL);
  g_assert_true (NC_IS_HALO_CM_BHATTACHARYA13 (hcmb));

  hcmb2 = nc_halo_cm_bhattacharya13_ref (hcmb);
  nc_halo_cm_bhattacharya13_clear (&hcmb2);
  g_assert_true (hcmb2 == NULL);

  g_assert_true (NC_IS_HALO_CM_BHATTACHARYA13 (hcmb));

  nc_distance_clear (&dist);
  nc_transfer_func_clear (&tf);
  nc_powspec_ml_clear (&ps_ml);
  ncm_powspec_filter_clear (&psf);
  nc_multiplicity_func_clear (&mf);
  nc_halo_mass_function_clear (&hmf);
  nc_growth_func_clear (&gf);

  NCM_TEST_FREE (nc_halo_cm_bhattacharya13_free, hcmb);
}

void
test_nc_halo_cm_diemer15_basic (void)
{
  NcDistance *dist        = nc_distance_new (1100.0);
  NcTransferFunc *tf      = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml      = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf   = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mf  = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new ());
  NcHaloMassFunction *hmf = nc_halo_mass_function_new (dist, psf, mf);
  NcHaloCMDiemer15 *hcmd;
  NcHaloCMDiemer15 *hcmd2;

  hcmd = nc_halo_cm_diemer15_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0, hmf);

  g_assert_true (hcmd != NULL);
  g_assert_true (NC_IS_HALO_CM_DIEMER15 (hcmd));

  hcmd2 = nc_halo_cm_diemer15_ref (hcmd);
  nc_halo_cm_diemer15_clear (&hcmd2);
  g_assert_true (hcmd2 == NULL);

  g_assert_true (NC_IS_HALO_CM_DIEMER15 (hcmd));

  nc_distance_clear (&dist);
  nc_transfer_func_clear (&tf);
  nc_powspec_ml_clear (&ps_ml);
  ncm_powspec_filter_clear (&psf);
  nc_multiplicity_func_clear (&mf);
  nc_halo_mass_function_clear (&hmf);

  NCM_TEST_FREE (nc_halo_cm_diemer15_free, hcmd);
}

void
test_nc_halo_cm_prada12_basic (void)
{
  NcDistance *dist        = nc_distance_new (1100.0);
  NcTransferFunc *tf      = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml      = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf   = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mf  = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_tinker_new ());
  NcHaloMassFunction *hmf = nc_halo_mass_function_new (dist, psf, mf);
  NcHaloCMPrada12 *hcmp;
  NcHaloCMPrada12 *hcmp2;

  hcmp = nc_halo_cm_prada12_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0, hmf);

  g_assert_true (hcmp != NULL);
  g_assert_true (NC_IS_HALO_CM_PRADA12 (hcmp));

  hcmp2 = nc_halo_cm_prada12_ref (hcmp);
  nc_halo_cm_prada12_clear (&hcmp2);
  g_assert_true (hcmp2 == NULL);

  g_assert_true (NC_IS_HALO_CM_PRADA12 (hcmp));

  nc_distance_clear (&dist);
  nc_transfer_func_clear (&tf);
  nc_powspec_ml_clear (&ps_ml);
  ncm_powspec_filter_clear (&psf);
  nc_multiplicity_func_clear (&mf);
  nc_halo_mass_function_clear (&hmf);

  NCM_TEST_FREE (nc_halo_cm_prada12_free, hcmp);
}

void
test_nc_halo_cm_dutton14_basic (void)
{
  NcHaloCMDutton14 *hcmd = nc_halo_cm_dutton14_new (NC_HALO_MASS_SUMMARY_MASS_DEF_CRITICAL, 200.0);
  NcHaloCMDutton14 *hcmd2;

  g_assert_true (hcmd != NULL);
  g_assert_true (NC_IS_HALO_CM_DUTTON14 (hcmd));

  hcmd2 = nc_halo_cm_dutton14_ref (hcmd);
  nc_halo_cm_dutton14_clear (&hcmd2);
  g_assert_true (hcmd2 == NULL);

  g_assert_true (NC_IS_HALO_CM_DUTTON14 (hcmd));

  NCM_TEST_FREE (nc_halo_cm_dutton14_free, hcmd);
}

void
test_nc_xcor_basic (void)
{
  NcmPowspec *ps   = NCM_POWSPEC (nc_powspec_ml_cbe_new ());
  NcDistance *dist = nc_distance_new (1100.0);
  NcXcor *xc       = nc_xcor_new (dist, ps, NC_XCOR_METHOD_LIMBER_Z_GSL);
  NcXcor *xc2;

  g_assert_true (xc != NULL);
  g_assert_true (NC_IS_XCOR (xc));

  xc2 = nc_xcor_ref (xc);
  nc_xcor_clear (&xc2);
  g_assert_true (xc2 == NULL);

  g_assert_true (NC_IS_XCOR (xc));

  ncm_powspec_clear (&ps);
  nc_distance_clear (&dist);

  NCM_TEST_FREE (nc_xcor_free, xc);
}

void
test_nc_halo_bias_despali_basic (void)
{
  NcDistance *dist        = nc_distance_new (1100.0);
  NcTransferFunc *tf      = nc_transfer_func_eh_new ();
  NcPowspecML *ps_ml      = NC_POWSPEC_ML (nc_powspec_ml_transfer_new (tf));
  NcmPowspecFilter *psf   = ncm_powspec_filter_new (NCM_POWSPEC (ps_ml), NCM_POWSPEC_FILTER_TYPE_TOPHAT);
  NcMultiplicityFunc *mf  = NC_MULTIPLICITY_FUNC (nc_multiplicity_func_despali_new ());
  NcHaloMassFunction *hmf = nc_halo_mass_function_new (dist, psf, mf);
  NcHaloBiasDespali *hbd  = nc_halo_bias_despali_new (hmf);
  NcHaloBiasDespali *hbd2;

  g_assert_true (hbd != NULL);
  g_assert_true (NC_IS_HALO_BIAS_DESPALI (hbd));

  hbd2 = nc_halo_bias_despali_ref (hbd);
  nc_halo_bias_despali_clear (&hbd2);
  g_assert_true (hbd2 == NULL);

  g_assert_true (NC_IS_HALO_BIAS_DESPALI (hbd));

  nc_distance_clear (&dist);
  nc_transfer_func_clear (&tf);
  nc_powspec_ml_clear (&ps_ml);
  ncm_powspec_filter_clear (&psf);
  nc_multiplicity_func_clear (&mf);
  nc_halo_mass_function_clear (&hmf);

  NCM_TEST_FREE (nc_halo_bias_despali_free, hbd);
}

void
test_nc_multiplicity_func_bhattacharya_basic (void)
{
  NcMultiplicityFuncBhattacharya *mbt = nc_multiplicity_func_bhattacharya_new ();
  NcMultiplicityFunc *mulf            = NC_MULTIPLICITY_FUNC (mbt);
  NcHICosmo *cosmo                    = NC_HICOSMO (nc_hicosmo_lcdm_new ());
  NcMultiplicityFuncBhattacharya *mbt2;

  g_assert_true (mbt != NULL);
  g_assert_true (NC_IS_MULTIPLICITY_FUNC_BHATTACHARYA (mbt));
  g_assert_cmpint (nc_multiplicity_func_get_mdef (mulf), ==, NC_MULTIPLICITY_FUNC_MASS_DEF_FOF);
  g_assert_cmpint (nc_multiplicity_func_bhattacharya_get_convention (mbt), ==, NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011);

  mbt2 = nc_multiplicity_func_bhattacharya_ref (mbt);
  nc_multiplicity_func_bhattacharya_clear (&mbt2);
  g_assert_true (mbt2 == NULL);

  g_assert_true (NC_IS_MULTIPLICITY_FUNC_BHATTACHARYA (mbt));

  nc_multiplicity_func_bhattacharya_set_A (mbt, 0.35);
  nc_multiplicity_func_bhattacharya_set_a (mbt, 0.8);
  nc_multiplicity_func_bhattacharya_set_p (mbt, 0.9);
  nc_multiplicity_func_bhattacharya_set_q (mbt, 1.7);
  nc_multiplicity_func_bhattacharya_set_delta_c (mbt, 1.686);

  g_assert_cmpfloat (nc_multiplicity_func_bhattacharya_get_A (mbt), ==, 0.35);
  g_assert_cmpfloat (nc_multiplicity_func_bhattacharya_get_a (mbt), ==, 0.8);
  g_assert_cmpfloat (nc_multiplicity_func_bhattacharya_get_p (mbt), ==, 0.9);
  g_assert_cmpfloat (nc_multiplicity_func_bhattacharya_get_q (mbt), ==, 1.7);
  g_assert_cmpfloat (nc_multiplicity_func_bhattacharya_get_delta_c (mbt), ==, 1.686);

  g_assert_cmpfloat (nc_multiplicity_func_eval (mulf, cosmo, 1.0, 0.5), >, 0.0);

  nc_hicosmo_clear (&cosmo);

  NCM_TEST_FREE (nc_multiplicity_func_bhattacharya_free, mbt);
}

void
test_nc_multiplicity_func_bhattacharya_convention (void)
{
  NcMultiplicityFuncBhattacharya *mbt_b = nc_multiplicity_func_bhattacharya_new_full (NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011);
  NcMultiplicityFuncBhattacharya *mbt_h = nc_multiplicity_func_bhattacharya_new_full (NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_HEITMANN2019);
  NcMultiplicityFunc *mulf_b            = NC_MULTIPLICITY_FUNC (mbt_b);
  NcMultiplicityFunc *mulf_h            = NC_MULTIPLICITY_FUNC (mbt_h);
  NcHICosmo *cosmo                      = NC_HICOSMO (nc_hicosmo_lcdm_new ());

  g_assert_cmpint (nc_multiplicity_func_bhattacharya_get_convention (mbt_b), ==, NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011);
  g_assert_cmpint (nc_multiplicity_func_bhattacharya_get_convention (mbt_h), ==, NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_HEITMANN2019);

  /* The two conventions agree at z = 0 and differ for z > 0. */
  ncm_assert_cmpdouble_e (nc_multiplicity_func_eval (mulf_b, cosmo, 1.0, 0.0), ==,
                          nc_multiplicity_func_eval (mulf_h, cosmo, 1.0, 0.0), 1.0e-15, 0.0);
  g_assert_cmpfloat (nc_multiplicity_func_eval (mulf_b, cosmo, 1.0, 1.0), !=,
                     nc_multiplicity_func_eval (mulf_h, cosmo, 1.0, 1.0));

  nc_multiplicity_func_bhattacharya_set_convention (mbt_h, NC_MULTIPLICITY_FUNC_BHATTACHARYA_CONVENTION_BHATTACHARYA2011);
  ncm_assert_cmpdouble_e (nc_multiplicity_func_eval (mulf_b, cosmo, 1.0, 1.0), ==,
                          nc_multiplicity_func_eval (mulf_h, cosmo, 1.0, 1.0), 1.0e-15, 0.0);

  nc_hicosmo_clear (&cosmo);
  nc_multiplicity_func_bhattacharya_free (mbt_b);
  NCM_TEST_FREE (nc_multiplicity_func_bhattacharya_free, mbt_h);
}

void
test_nc_multiplicity_func_bhattacharya_mean_mdef (void)
{
  if (g_test_subprocess ())
  {
    NcMultiplicityFuncBhattacharya *mbt = nc_multiplicity_func_bhattacharya_new ();

    nc_multiplicity_func_set_mdef (NC_MULTIPLICITY_FUNC (mbt), NC_MULTIPLICITY_FUNC_MASS_DEF_MEAN);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not support mean mass def*");
}

void
test_nc_multiplicity_func_bhattacharya_critical_mdef (void)
{
  if (g_test_subprocess ())
  {
    NcMultiplicityFuncBhattacharya *mbt = nc_multiplicity_func_bhattacharya_new ();

    nc_multiplicity_func_set_mdef (NC_MULTIPLICITY_FUNC (mbt), NC_MULTIPLICITY_FUNC_MASS_DEF_CRITICAL);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not support critical mass def*");
}

void
test_nc_multiplicity_func_bhattacharya_virial_mdef (void)
{
  if (g_test_subprocess ())
  {
    NcMultiplicityFuncBhattacharya *mbt = nc_multiplicity_func_bhattacharya_new ();

    nc_multiplicity_func_set_mdef (NC_MULTIPLICITY_FUNC (mbt), NC_MULTIPLICITY_FUNC_MASS_DEF_VIRIAL);

    return; /* LCOV_EXCL_LINE */
  }

  g_test_trap_subprocess (NULL, 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not support virial mass def*");
}

void
test_nc_cluster_mass_ascaso_basic (void)
{
  NcClusterMassAscaso *cluster_m = g_object_new (NC_TYPE_CLUSTER_MASS_ASCASO, NULL);
  NcClusterMassAscaso *cluster_m2;

  g_assert_true (cluster_m != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS_ASCASO (cluster_m));

  cluster_m2 = g_object_ref (cluster_m);
  g_clear_object (&cluster_m2);
  g_assert_true (cluster_m2 == NULL);

  g_assert_true (NC_IS_CLUSTER_MASS_ASCASO (cluster_m));

  NCM_TEST_FREE (g_object_unref, cluster_m);
}

void
test_nc_cluster_mass_selection_basic (void)
{
  NcClusterMassSelection *cluster_m = g_object_new (NC_TYPE_CLUSTER_MASS_SELECTION, NULL);
  NcClusterMassSelection *cluster_m2;

  g_assert_true (cluster_m != NULL);
  g_assert_true (NC_IS_CLUSTER_MASS_SELECTION (cluster_m));

  cluster_m2 = g_object_ref (cluster_m);
  g_clear_object (&cluster_m2);
  g_assert_true (cluster_m2 == NULL);

  g_assert_true (NC_IS_CLUSTER_MASS_SELECTION (cluster_m));

  NCM_TEST_FREE (g_object_unref, cluster_m);
}

void
test_nc_cluster_photoz_gauss_basic (void)
{
  NcClusterPhotozGauss *cluster_z = nc_cluster_photoz_gauss_new ();
  NcClusterPhotozGauss *cluster_z2;

  g_assert_true (cluster_z != NULL);
  g_assert_true (NC_IS_CLUSTER_PHOTOZ_GAUSS (cluster_z));
  cluster_z2 = g_object_ref (cluster_z);
  g_clear_object (&cluster_z2);
  g_assert_true (cluster_z2 == NULL);
  g_assert_true (NC_IS_CLUSTER_PHOTOZ_GAUSS (cluster_z));
  NCM_TEST_FREE (g_object_unref, cluster_z);
}

void
test_nc_galaxy_position_factor_flat_basic (void)
{
  NcGalaxyPositionFactorFlat *gpff = nc_galaxy_position_factor_flat_new (0.0, 2.0, -0.5, 0.5);
  NcGalaxyPositionFactorFlat *gpff2;

  g_assert_true (gpff != NULL);
  g_assert_true (NC_IS_GALAXY_POSITION_FACTOR_FLAT (gpff));

  gpff2 = nc_galaxy_position_factor_flat_ref (gpff);
  nc_galaxy_position_factor_flat_clear (&gpff2);
  g_assert_true (gpff2 == NULL);

  g_assert_true (NC_IS_GALAXY_POSITION_FACTOR_FLAT (gpff));

  NCM_TEST_FREE (nc_galaxy_position_factor_flat_free, gpff);
}

void
test_nc_galaxy_redshift_factor_composed_basic (void)
{
  NcGalaxyRedshiftFactorComposed *gsdrc = nc_galaxy_redshift_factor_composed_new (0.0, 5.0);
  NcGalaxyRedshiftFactorComposed *gsdrc2;

  g_assert_true (gsdrc != NULL);
  g_assert_true (NC_IS_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdrc));

  gsdrc2 = nc_galaxy_redshift_factor_composed_ref (gsdrc);
  nc_galaxy_redshift_factor_composed_clear (&gsdrc2);
  g_assert_true (gsdrc2 == NULL);

  g_assert_true (NC_IS_GALAXY_REDSHIFT_FACTOR_COMPOSED (gsdrc));

  NCM_TEST_FREE (nc_galaxy_redshift_factor_composed_free, gsdrc);
}

void
test_nc_galaxy_redshift_factor_spline_basic (void)
{
  NcGalaxyRedshiftFactorSpline *gsdrs = nc_galaxy_redshift_factor_spline_new ();
  NcGalaxyRedshiftFactorSpline *gsdrs2;

  g_assert_true (gsdrs != NULL);
  g_assert_true (NC_IS_GALAXY_REDSHIFT_FACTOR_SPLINE (gsdrs));

  gsdrs2 = nc_galaxy_redshift_factor_spline_ref (gsdrs);
  nc_galaxy_redshift_factor_spline_clear (&gsdrs2);
  g_assert_true (gsdrs2 == NULL);

  g_assert_true (NC_IS_GALAXY_REDSHIFT_FACTOR_SPLINE (gsdrs));

  NCM_TEST_FREE (nc_galaxy_redshift_factor_spline_free, gsdrs);
}

void
test_nc_galaxy_redshift_binning_basic (void)
{
  NcGalaxyRedshiftBinning *gsdrb = nc_galaxy_redshift_binning_new ();
  NcGalaxyRedshiftBinning *gsdrb2;

  g_assert_true (gsdrb != NULL);
  g_assert_true (NC_IS_GALAXY_REDSHIFT_BINNING (gsdrb));

  gsdrb2 = nc_galaxy_redshift_binning_ref (gsdrb);
  nc_galaxy_redshift_binning_clear (&gsdrb2);
  g_assert_true (gsdrb2 == NULL);

  g_assert_true (NC_IS_GALAXY_REDSHIFT_BINNING (gsdrb));

  NCM_TEST_FREE (nc_galaxy_redshift_binning_free, gsdrb);
}

void
test_nc_galaxy_redshift_pop_basic (void)
{
  /* NcGalaxyRedshiftPop is abstract: use a concrete subclass instance, but
   * exercise the base class's own ref/free/clear directly (every concrete
   * subclass conventionally shadows these with its own same-named wrapper,
   * so a subclass-specific _basic test would never reach the base's). */
  NcGalaxyRedshiftPop *gsdrp = NC_GALAXY_REDSHIFT_POP (nc_galaxy_redshift_pop_lsst_srd_new_y1_source ());
  NcGalaxyRedshiftPop *gsdrp2;

  g_assert_true (gsdrp != NULL);
  g_assert_true (NC_IS_GALAXY_REDSHIFT_POP (gsdrp));

  gsdrp2 = nc_galaxy_redshift_pop_ref (gsdrp);
  nc_galaxy_redshift_pop_clear (&gsdrp2);
  g_assert_true (gsdrp2 == NULL);

  g_assert_true (NC_IS_GALAXY_REDSHIFT_POP (gsdrp));

  NCM_TEST_FREE (nc_galaxy_redshift_pop_free, gsdrp);
}

void
test_nc_galaxy_shape_factor_basic (void)
{
  /* NcGalaxyShapeFactor is abstract: same rationale as
   * test_nc_galaxy_redshift_pop_basic above. */
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_quad_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));
  NcGalaxyShapeFactor *gsf2;

  g_assert_true (gsf != NULL);
  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  gsf2 = nc_galaxy_shape_factor_ref (gsf);
  nc_galaxy_shape_factor_clear (&gsf2);
  g_assert_true (gsf2 == NULL);

  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR (gsf));

  NCM_TEST_FREE (nc_galaxy_shape_factor_free, gsf);
}

void
test_nc_galaxy_shape_factor_quad_basic (void)
{
  NcGalaxyShapeFactorQuad *gsfq = nc_galaxy_shape_factor_quad_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapeFactorQuad *gsfq2;

  g_assert_true (gsfq != NULL);
  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_QUAD (gsfq));

  gsfq2 = nc_galaxy_shape_factor_quad_ref (gsfq);
  nc_galaxy_shape_factor_quad_clear (&gsfq2);
  g_assert_true (gsfq2 == NULL);

  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_QUAD (gsfq));

  NCM_TEST_FREE (nc_galaxy_shape_factor_quad_free, gsfq);
}

void
test_nc_galaxy_shape_factor_series_lensed_basic (void)
{
  NcGalaxyShapeFactorSeriesLensed *gsfsl = nc_galaxy_shape_factor_series_lensed_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, 4);
  NcGalaxyShapeFactorSeriesLensed *gsfsl2;

  g_assert_true (gsfsl != NULL);
  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_SERIES_LENSED (gsfsl));

  gsfsl2 = nc_galaxy_shape_factor_series_lensed_ref (gsfsl);
  nc_galaxy_shape_factor_series_lensed_clear (&gsfsl2);
  g_assert_true (gsfsl2 == NULL);

  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_SERIES_LENSED (gsfsl));

  NCM_TEST_FREE (nc_galaxy_shape_factor_series_lensed_free, gsfsl);
}

void
test_nc_galaxy_shape_factor_fixed_quad_basic (void)
{
  NcGalaxyShapeFactorFixedQuad *gsffq = nc_galaxy_shape_factor_fixed_quad_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapeFactorFixedQuad *gsffq2;

  g_assert_true (gsffq != NULL);
  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_FIXED_QUAD (gsffq));

  gsffq2 = nc_galaxy_shape_factor_fixed_quad_ref (gsffq);
  nc_galaxy_shape_factor_fixed_quad_clear (&gsffq2);
  g_assert_true (gsffq2 == NULL);

  g_assert_true (NC_IS_GALAXY_SHAPE_FACTOR_FIXED_QUAD (gsffq));

  NCM_TEST_FREE (nc_galaxy_shape_factor_fixed_quad_free, gsffq);
}

void
test_nc_data_cluster_wl_factor_basic (void)
{
  NcGalaxyPositionFactor *position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (0.0, 2.0, -0.5, 0.5));
  NcGalaxyRedshiftFactor *redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 5.0));
  NcGalaxyShapeFactor *shape_factor       = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_quad_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET));
  NcDataClusterWLFactor *dcwlf            = nc_data_cluster_wl_factor_new (position_factor, redshift_factor, shape_factor);
  NcDataClusterWLFactor *dcwlf2;

  g_assert_true (dcwlf != NULL);
  g_assert_true (NC_IS_DATA_CLUSTER_WL_FACTOR (dcwlf));

  dcwlf2 = nc_data_cluster_wl_factor_ref (dcwlf);
  nc_data_cluster_wl_factor_clear (&dcwlf2);
  g_assert_true (dcwlf2 == NULL);

  g_assert_true (NC_IS_DATA_CLUSTER_WL_FACTOR (dcwlf));

  /* new() dups its own refs via the construct properties. */
  nc_galaxy_position_factor_clear (&position_factor);
  nc_galaxy_redshift_factor_clear (&redshift_factor);
  nc_galaxy_shape_factor_clear (&shape_factor);

  NCM_TEST_FREE (nc_data_cluster_wl_factor_free, dcwlf);
}

