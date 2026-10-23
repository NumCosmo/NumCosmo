/***************************************************************************
 *            test_nc_galaxy_shape_factor_moments_tilt.c
 *
 *  Fri Oct 23 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2026 <vitenti@uel.br>
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

/*
 * Tests of NcGalaxyShapeFactorMomentsTilt:newton-maxit (a table built with a
 * small iteration cap, where solves end at the cap, agrees with the default
 * table and reports no failed solve) and of NcGalaxyShapeFactorMomentsTilt:
 * restrict-range (the bound is dropped when a parameter it does not cover
 * is free).
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>

#define MOMENTS_TILT_SIGMA_POP (0.3)
#define MOMENTS_TILT_STD_NOISE (0.1)
#define MOMENTS_TILT_SMALL_MAXIT (2)
#define MOMENTS_TILT_DEFAULT_MAXIT (40)
#define MOMENTS_TILT_RANGE_SIGMA_POP (0.4658)
#define MOMENTS_TILT_RANGE_LOG10M_UPPER (14.5)
#define MOMENTS_TILT_RANGE_RA (0.8)
#define MOMENTS_TILT_RANGE_Z_MAX (3.0)

typedef struct _TestNcGalaxyShapeFactorMomentsTilt
{
  NcmMSet *mset;
  NcGalaxyShapePop *pop;
  NcGalaxyPositionFactor *position_factor;
  NcGalaxyRedshiftFactor *redshift_factor;
} TestNcGalaxyShapeFactorMomentsTilt;

void test_nc_galaxy_shape_factor_moments_tilt_new (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_moments_tilt_free (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_default (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_cap (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_moments_tilt_restrict_range_free_params (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/nc/galaxy_shape_factor_moments_tilt/newton_maxit/default", TestNcGalaxyShapeFactorMomentsTilt, NULL,
              &test_nc_galaxy_shape_factor_moments_tilt_new,
              &test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_default,
              &test_nc_galaxy_shape_factor_moments_tilt_free);

  g_test_add ("/nc/galaxy_shape_factor_moments_tilt/newton_maxit/cap", TestNcGalaxyShapeFactorMomentsTilt, NULL,
              &test_nc_galaxy_shape_factor_moments_tilt_new,
              &test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_cap,
              &test_nc_galaxy_shape_factor_moments_tilt_free);

  g_test_add ("/nc/galaxy_shape_factor_moments_tilt/restrict_range/free_params", TestNcGalaxyShapeFactorMomentsTilt, NULL,
              &test_nc_galaxy_shape_factor_moments_tilt_new,
              &test_nc_galaxy_shape_factor_moments_tilt_restrict_range_free_params,
              &test_nc_galaxy_shape_factor_moments_tilt_free);

  g_test_run ();

  return 0;
}

/* Builds a MomentsTilt with @newton_maxit and @restrict_range, and the
 * per-galaxy data of one galaxy at right ascension @ra (degrees from the
 * halo center), prepared. */
static NcGalaxyShapeFactorMomentsTilt *
_moments_tilt_build (TestNcGalaxyShapeFactorMomentsTilt *test, const guint newton_maxit, const gboolean restrict_range,
                     const gdouble ra, NcGalaxyShapeFactorData **data)
{
  NcGalaxyShapeFactorMomentsTilt *gsfmt = g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_MOMENTS_TILT,
                                                        "ellip-conv", NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET,
                                                        "restrict-range", restrict_range,
                                                        "newton-maxit", newton_maxit,
                                                        NULL);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  GPtrArray *data_array              = g_ptr_array_new ();

  p_data->ra  = ra;
  p_data->dec = 0.0;

  *data = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (NC_GALAXY_SHAPE_FACTOR (gsfmt), *data, 0.0, 0.0, MOMENTS_TILT_STD_NOISE, 0.0, 0.0, 0.0,
                                   NC_WL_ELLIPTICITY_FRAME_CELESTIAL);

  g_ptr_array_add (data_array, *data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset, data_array, TRUE, TRUE));

  g_ptr_array_unref (data_array);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);

  return gsfmt;
}

/* Returns the top of the table @gsfmt holds for @data. */
static gdouble
_moments_tilt_top (NcGalaxyShapeFactorMomentsTilt *gsfmt, TestNcGalaxyShapeFactorMomentsTilt *test, NcGalaxyShapeFactorData *data)
{
  guint n_panels  = 0;
  gdouble top     = 0.0;
  GArray *degrees = NULL;

  nc_galaxy_shape_factor_moments_tilt_peek_layout (gsfmt, test->pop, data, &n_panels, &top, &degrees);
  g_array_unref (degrees);

  return top;
}

void
test_nc_galaxy_shape_factor_moments_tilt_new (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata)
{
  NcHICosmo *cosmo                  = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                  = nc_distance_new (100.0);
  NcHaloMassSummary *hms            = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp          = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd       = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z   = nc_galaxy_redshift_obs_gauss_new ();

  test->pop = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, MOMENTS_TILT_SIGMA_POP);
  ncm_model_param_set (NCM_MODEL (hms), NC_HALO_CM_PARAM_LOG10M_DELTA, 14.0);
  ncm_model_param_set (NCM_MODEL (hp), NC_HALO_POSITION_Z, 0.2);
  nc_halo_position_prepare (hp, cosmo);

  test->mset = ncm_mset_empty_new ();
  ncm_mset_set (test->mset, NCM_MODEL (cosmo), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (dp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (hp), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (smd), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (test->pop), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (pop_z), NULL);
  ncm_mset_set (test->mset, NCM_MODEL (obs_z), NULL);

  test->position_factor = NC_GALAXY_POSITION_FACTOR (nc_galaxy_position_factor_flat_new (-1.0, 1.0, -1.0, 1.0));
  test->redshift_factor = NC_GALAXY_REDSHIFT_FACTOR (nc_galaxy_redshift_factor_composed_new (0.0, 20.0));

  nc_hicosmo_free (cosmo);
  nc_distance_free (dist);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
  nc_halo_position_free (hp);
  nc_wl_surface_mass_density_free (smd);
  nc_galaxy_redshift_pop_lsst_srd_free (pop_z);
  nc_galaxy_redshift_obs_gauss_free (obs_z);
}

void
test_nc_galaxy_shape_factor_moments_tilt_free (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  nc_galaxy_shape_pop_clear (&test->pop);
  nc_galaxy_position_factor_clear (&test->position_factor);
  nc_galaxy_redshift_factor_clear (&test->redshift_factor);
}

void
test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_default (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorMomentsTilt *gsfmt = nc_galaxy_shape_factor_moments_tilt_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  guint newton_maxit                    = 0;

  g_object_get (gsfmt, "newton-maxit", &newton_maxit, NULL);
  g_assert_cmpuint (newton_maxit, ==, MOMENTS_TILT_DEFAULT_MAXIT);

  nc_galaxy_shape_factor_moments_tilt_free (gsfmt);
}

/* With newton-maxit = 2 most node solves end at the cap. Each must return
 * the log-normalization of the lambda it returns, so no solve is reported
 * as failed and ln P agrees with the default table to within accuracy-gate
 * (measured: 3.4e-6 against a gate of 1e-3). */
void
test_nc_galaxy_shape_factor_moments_tilt_newton_maxit_cap (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata)
{
  const gdouble ghat[]  = { 0.2, 0.6, 0.9, 0.99 };
  const gdouble e_obs[] = { 0.3, 0.2, -0.4, 0.5, 0.05, -0.1, 0.8, 0.1 };
  NcGalaxyShapeFactorData *data_ref, *data_cap;
  NcGalaxyShapeFactorMomentsTilt *gsf_ref = _moments_tilt_build (test, MOMENTS_TILT_DEFAULT_MAXIT, FALSE, 0.0, &data_ref);
  NcGalaxyShapeFactorMomentsTilt *gsf_cap = _moments_tilt_build (test, MOMENTS_TILT_SMALL_MAXIT, FALSE, 0.0, &data_cap);
  gdouble gate                            = 0.0;
  gdouble max_diff                        = 0.0;
  guint i, j;

  g_object_get (gsf_cap, "accuracy-gate", &gate, NULL);

  for (i = 0; i < G_N_ELEMENTS (ghat); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (e_obs); j += 2)
    {
      const gdouble ln_ref = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsf_ref), test->pop, data_ref, ghat[i], 0.0, e_obs[j], e_obs[j + 1]);
      const gdouble ln_cap = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsf_cap), test->pop, data_cap, ghat[i], 0.0, e_obs[j], e_obs[j + 1]);

      g_assert_true (gsl_finite (ln_cap));
      max_diff = MAX (max_diff, fabs (ln_cap - ln_ref));
    }
  }

  g_assert_cmpuint (nc_galaxy_shape_factor_moments_tilt_get_solve_error_count (gsf_cap), ==, 0);
  g_assert_cmpfloat (max_diff, <, gate);

  nc_galaxy_shape_factor_data_unref (data_ref);
  nc_galaxy_shape_factor_data_unref (data_cap);
  nc_galaxy_shape_factor_moments_tilt_free (gsf_ref);
  nc_galaxy_shape_factor_moments_tilt_free (gsf_cap);
}

/* A galaxy that cannot reach the critical curve for any mass below the
 * log10MDelta upper bound gets a table cut below 1. Freeing a parameter the
 * bound does not cover (Omega_c, then cDelta) restores the full table;
 * fixing it again restores the cut; freeing log10MDelta keeps the cut. */
void
test_nc_galaxy_shape_factor_moments_tilt_restrict_range_free_params (TestNcGalaxyShapeFactorMomentsTilt *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcGalaxyShapeFactorData *data;
  NcGalaxyShapeFactorMomentsTilt *gsfmt;
  guint log10M_idx = 0;
  guint cDelta_idx = 0;
  gdouble top_cut;

  ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, MOMENTS_TILT_RANGE_SIGMA_POP);
  g_assert_true (ncm_model_param_index_from_name (hms, "log10MDelta", &log10M_idx, NULL));
  g_assert_true (ncm_model_param_index_from_name (hms, "cDelta", &cDelta_idx, NULL));
  ncm_model_param_set_upper_bound (hms, log10M_idx, MOMENTS_TILT_RANGE_LOG10M_UPPER);

  gsfmt = _moments_tilt_build (test, MOMENTS_TILT_DEFAULT_MAXIT, TRUE, MOMENTS_TILT_RANGE_RA, &data);
  nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset, data, MOMENTS_TILT_RANGE_Z_MAX);

  top_cut = _moments_tilt_top (gsfmt, test, data);
  g_assert_cmpfloat (top_cut, <, 1.0);

  ncm_mset_param_set_ftype (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset);
  g_assert_cmpfloat (_moments_tilt_top (gsfmt, test, data), ==, 1.0);

  ncm_mset_param_set_ftype (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset);
  g_assert_cmpfloat (_moments_tilt_top (gsfmt, test, data), ==, top_cut);

  ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), log10M_idx, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset);
  g_assert_cmpfloat (_moments_tilt_top (gsfmt, test, data), ==, top_cut);

  ncm_mset_param_set_ftype (test->mset, nc_halo_mass_summary_id (), cDelta_idx, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfmt), test->mset);
  g_assert_cmpfloat (_moments_tilt_top (gsfmt, test, data), ==, 1.0);

  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_moments_tilt_free (gsfmt);
}

