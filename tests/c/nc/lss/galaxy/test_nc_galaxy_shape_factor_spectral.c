/***************************************************************************
 *            test_nc_galaxy_shape_factor_spectral.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 *  Copyright  2026  Caio Lima de Oliveira
 *  <caiolimadeoliveira@pm.me>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_factor_spectral.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
 * Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
 *
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
 * Tests of NcGalaxyShapeFactorSpectral: ln P against the values certified by the Arb
 * reference tests/tools/nc_galaxy_shape_marginal_arb.c, the |g| = 1 closed form, the
 * |g| > 1 reduction, shears next to |g| = 1, the rotation of the frame, the hooks of NcGalaxyShapeFactor, the
 * accepted panels, the interpolant of ln P in g, the properties and the errors for an input
 * that cannot be integrated and for an interpolation range that is not accepted.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "test_nc_galaxy_shape_marginal_arb_cases.h"

#include <math.h>
#include <string.h>
#include <glib.h>
#include <glib-object.h>

/*
 * Largest |ln P - ln P_Arb| accepted at the default reltol (1e-5). The largest difference
 * measured is 6.9e-9, on the 20,062 certified points of the certification grid with
 * (|eps_obs| - 1) / std_noise <= 8.
 */
#define TEST_SPECTRAL_TOL (1.0e-7)

/*
 * A loop that vectorizes and a scalar computation of the same ln P agree within rounding,
 * not bit for bit (nc_galaxy_shape_factor_eval_nodes()): the compiler may contract
 * a * b + c into an fma differently in each, and a difference in the last bit of an
 * intermediate is amplified by what follows. They are compared within 1e-13 of
 * max (1, |ln P|).
 */
#define TEST_SPECTRAL_ASSERT_LOOP(loop, scalar) g_assert_cmpfloat (fabs ((loop) - (scalar)), <=, 1.0e-13 * MAX (1.0, fabs (scalar)))

/*
 * Largest |ln P (interpolant) - ln P (integral)| accepted at the default interp-abstol
 * (1e-5). The largest difference measured is 1.0e-8, on 150 HSC galaxies in both
 * conventions with Beta and truncated Gaussian populations.
 */
#define TEST_SPECTRAL_INTERP_TOL (1.0e-7)

/*
 * interp-abstol of the prepared-catalog tests, which check the record, the loading and the
 * mismatch aborts, not accuracy; /interp and /interp/rebuild keep the default.
 */
#define TEST_SPECTRAL_OBS_INTERP_ABSTOL (1.0e-2)

/*
 * Points with 1 - |g| from 1e-5 to 1e-8, certified by the Arb reference with --coords=uv
 * (radius of ln P below 1e-15): HWL16a-002 rows 475 and 7139 with the Beta population
 * fitted to that catalog. Below 1 - |g| = 1e-5 the class forms the value from integrals at
 * that distance. The largest difference of the class measured on them is 1.9e-9.
 */
static const TestShapeArbCell test_spectral_near_unit_cells[] = {
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591169317498968, 0.28297475387508314, 0.691024972903032, 0.0, -18.15200444413578233},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591264271024746, 0.28297755535316127, 0.691024972903032, 0.0, -18.15200449764241692},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591169317498968, 0.28297475387508314, 0.691024972903032, 0.0, -18.14997346866772457},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591264271024746, 0.28297755535316127, 0.691024972903032, 0.0, -18.15199473594852941},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.999999, 0.0, 0.178728848695755, 0.3698130249977112, -59.93454784752150307},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.99999999, 0.0, 0.178728848695755, 0.3698130249977112, -79.77604456776660732},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.99999999, 0.0, 0.178728848695755, 0.3698130249977112, -47.75380931482856171},
};

/* The changes of the fit after a catalog is prepared, and the message each aborts with */
static const gchar *test_spectral_obs_mismatch_cases[][2] = {
  {"pop",        "*built for other settings:*the population (NcGalaxyShapePopBeta*"},
  {"free",       "*built for other settings:*the population parameter * is free*"},
  {"abstol",     "*built for other settings:*interp-abstol*"},
  {"conv",       "*built for other settings:*the ellipticity convention*"},
  {"use_interp", "*the catalog carries interpolants, but use-interp is not set*"},
  {"eval",       "*built for other settings:*the population (NcGalaxyShapePopBeta*"},
};

/* The galaxies of the prepared catalogs: observed ellipticity and noise */
static const gdouble test_spectral_obs_gal[][3] = {
  { 0.3,   0.2, 0.1},
  { 0.05, -0.4, 0.05},
  {-0.6,   0.1, 0.2},
};

typedef struct _TestNcGalaxyShapeFactorSpectral
{
  NcmMSet *mset;
  NcGalaxyShapePop *pop;
  NcGalaxyPositionFactor *position_factor;
  NcGalaxyRedshiftFactor *redshift_factor;
} TestNcGalaxyShapeFactorSpectral;

void test_nc_galaxy_shape_factor_spectral_new (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_free (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);

void test_nc_galaxy_shape_factor_spectral_arb (void);
void test_nc_galaxy_shape_factor_spectral_unit_shear (void);
void test_nc_galaxy_shape_factor_spectral_large_shear (void);
void test_nc_galaxy_shape_factor_spectral_near_unit (void);
void test_nc_galaxy_shape_factor_spectral_near_unit_sides (void);
void test_nc_galaxy_shape_factor_spectral_rotation (void);
void test_nc_galaxy_shape_factor_spectral_reltol (void);
void test_nc_galaxy_shape_factor_spectral_panels (void);
void test_nc_galaxy_shape_factor_spectral_serialize (void);
void test_nc_galaxy_shape_factor_spectral_hooks (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_crit_hash (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_nodes (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_array (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_rebuild (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_data_prepare (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_obs (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch (void);
void test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch_subprocess (gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_g_max (void);
void test_nc_galaxy_shape_factor_spectral_restrict_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_stable (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_short (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_same_box (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_small_g_max (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_pop_order (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_pop_order_diverge (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_pop_order_off (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_pop_order_fallback_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_g_max_subprocess (void);
void test_nc_galaxy_shape_factor_spectral_invalid_noise (void);
void test_nc_galaxy_shape_factor_spectral_invalid_noise_subprocess (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add_func ("/nc/galaxy_shape_factor_spectral/arb", &test_nc_galaxy_shape_factor_spectral_arb);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/unit_shear", &test_nc_galaxy_shape_factor_spectral_unit_shear);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/large_shear", &test_nc_galaxy_shape_factor_spectral_large_shear);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/near_unit", &test_nc_galaxy_shape_factor_spectral_near_unit);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/near_unit/sides", &test_nc_galaxy_shape_factor_spectral_near_unit_sides);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/rotation", &test_nc_galaxy_shape_factor_spectral_rotation);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/reltol", &test_nc_galaxy_shape_factor_spectral_reltol);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/panels", &test_nc_galaxy_shape_factor_spectral_panels);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/serialize", &test_nc_galaxy_shape_factor_spectral_serialize);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_nodes", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_nodes,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_line", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_line,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/line", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_line,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/crit_hash", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_crit_hash,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/hooks", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_hooks,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_array", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_array,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/rebuild", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_rebuild,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/data_prepare", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_data_prepare,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/obs", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_obs,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/obs/mismatch", &test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch);
  {
    guint k;

    for (k = 0; k < G_N_ELEMENTS (test_spectral_obs_mismatch_cases); k++)
    {
      gchar *path = g_strdup_printf ("/nc/galaxy_shape_factor_spectral/interp/obs/mismatch/%s/subprocess", test_spectral_obs_mismatch_cases[k][0]);

      g_test_add_data_func (path, test_spectral_obs_mismatch_cases[k][0], &test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch_subprocess);
      g_free (path);
    }
  }
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/g_max", &test_nc_galaxy_shape_factor_spectral_interp_g_max);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/restrict_range", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/restrict_range/stable", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_stable,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/restrict_range/loaded_short", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_short,
              &test_nc_galaxy_shape_factor_spectral_free);

  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/restrict_range/loaded_same_box", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_same_box,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/small_g_max", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_small_g_max,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop_order", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_pop_order,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop_order/diverge", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_pop_order_diverge,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop_order/off", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_pop_order_off,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop_order/fallback_line", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_pop_order_fallback_line,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/g_max/subprocess", &test_nc_galaxy_shape_factor_spectral_interp_g_max_subprocess);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/invalid_noise", &test_nc_galaxy_shape_factor_spectral_invalid_noise);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/invalid_noise/subprocess", &test_nc_galaxy_shape_factor_spectral_invalid_noise_subprocess);

  g_test_run ();

  return 0;
}

/* A new population of kind pop with parameters p1 and p2, and its prepared data in *data */
static NcGalaxyShapePop *
_test_pop_new (const TestShapeArbPop pop, const gdouble p1, const gdouble p2, NcGalaxyShapePopData **data)
{
  NcGalaxyShapePop *gsp;

  if (pop == TEST_SHAPE_ARB_POP_GAUSS)
  {
    gsp = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
    ncm_model_param_set_by_name (NCM_MODEL (gsp), "sigma", p1, NULL);
  }
  else
  {
    gsp = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ());
    ncm_model_param_set_by_name (NCM_MODEL (gsp), "alpha", p1, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (gsp), "beta", p2, NULL);
  }

  *data = nc_galaxy_shape_pop_data_new (gsp);
  nc_galaxy_shape_pop_prepare (gsp, *data);

  return gsp;
}

void
test_nc_galaxy_shape_factor_spectral_new (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHICosmo *cosmo                  = NC_HICOSMO (nc_hicosmo_de_xcdm_new ());
  NcDistance *dist                  = nc_distance_new (100.0);
  NcHaloMassSummary *hms            = NC_HALO_MASS_SUMMARY (nc_halo_cm_param_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp          = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  NcHaloPosition *hp                = nc_halo_position_new (dist);
  NcWLSurfaceMassDensity *smd       = nc_wl_surface_mass_density_new (dist);
  NcGalaxyRedshiftPopLSSTSRD *pop_z = nc_galaxy_redshift_pop_lsst_srd_new_y1_source ();
  NcGalaxyRedshiftObsGauss *obs_z   = nc_galaxy_redshift_obs_gauss_new ();

  test->pop = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ());
  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "alpha", 1.57, NULL);
  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "beta", 1.65, NULL);
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
test_nc_galaxy_shape_factor_spectral_free (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  ncm_mset_clear (&test->mset);
  nc_galaxy_shape_pop_clear (&test->pop);
  nc_galaxy_position_factor_clear (&test->position_factor);
  nc_galaxy_redshift_factor_clear (&test->redshift_factor);
}

/* Every row of the certified table, at the default reltol */
void
test_nc_galaxy_shape_factor_spectral_arb (void)
{
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_shape_arb_cells); i++)
  {
    const TestShapeArbCell *c         = &test_shape_arb_cells[i];
    NcGalaxyShapePopData *pop_data    = NULL;
    NcGalaxyShapePop *pop             = _test_pop_new (c->pop, c->p1, c->p2, &pop_data);
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (c->conv);
    const gdouble lnP                 = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, c->g1, c->g2, c->e1, c->e2, c->std_noise);

    g_assert_cmpfloat_with_epsilon (lnP, c->lnP, TEST_SPECTRAL_TOL);

    nc_galaxy_shape_factor_spectral_free (gsfs);
    nc_galaxy_shape_pop_data_unref (pop_data);
    nc_galaxy_shape_pop_free (pop);
  }
}

/* At |g| = 1, ln P = ln N_2 (eps_obs - g), for both conventions */
void
test_nc_galaxy_shape_factor_spectral_unit_shear (void)
{
  const gdouble g[]              = {1.0, 0.0, 0.0, -1.0};
  const gdouble e[]              = {0.3, -0.2, 0.9, 0.4};
  const gdouble s                = 0.05;
  NcGalaxyShapePopData *pop_data = NULL;
  NcGalaxyShapePop *pop          = _test_pop_new (TEST_SHAPE_ARB_POP_BETA, 1.57, 1.65, &pop_data);
  guint j, k;

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);

    for (k = 0; k < G_N_ELEMENTS (g); k += 2)
    {
      const gdouble d2     = gsl_pow_2 (e[k] - g[k]) + gsl_pow_2 (e[k + 1] - g[k + 1]);
      const gdouble ln_ref = -d2 / (2.0 * s * s) - log (2.0 * M_PI * s * s);
      const gdouble lnP    = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g[k], g[k + 1], e[k], e[k + 1], s);

      g_assert_cmpfloat_with_epsilon (lnP, ln_ref, 1.0e-13 * fabs (ln_ref));
    }

    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/* For |g| > 1, ln P is the one at g / |g|^2 */
void
test_nc_galaxy_shape_factor_spectral_large_shear (void)
{
  const gdouble g[]              = {1.5, 0.0, -0.7, 2.9};
  const gdouble e1               = 0.3;
  const gdouble e2               = -0.4;
  const gdouble s                = 0.1;
  NcGalaxyShapePopData *pop_data = NULL;
  NcGalaxyShapePop *pop          = _test_pop_new (TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, &pop_data);
  guint j, k;

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);

    for (k = 0; k < G_N_ELEMENTS (g); k += 2)
    {
      const gdouble abs_g2 = gsl_pow_2 (g[k]) + gsl_pow_2 (g[k + 1]);
      const gdouble large  = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g[k], g[k + 1], e1, e2, s);
      const gdouble small  = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g[k] / abs_g2, g[k + 1] / abs_g2, e1, e2, s);

      g_assert_cmpfloat_with_epsilon (large, small, 1.0e-10);
    }

    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/* Every row of the table of shears next to |g| = 1, at the default reltol */
void
test_nc_galaxy_shape_factor_spectral_near_unit (void)
{
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_spectral_near_unit_cells); i++)
  {
    const TestShapeArbCell *c         = &test_spectral_near_unit_cells[i];
    NcGalaxyShapePopData *pop_data    = NULL;
    NcGalaxyShapePop *pop             = _test_pop_new (c->pop, c->p1, c->p2, &pop_data);
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (c->conv);
    const gdouble lnP                 = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, c->g1, c->g2, c->e1, c->e2, c->std_noise);

    g_assert_cmpfloat_with_epsilon (lnP, c->lnP, TEST_SPECTRAL_TOL);

    nc_galaxy_shape_factor_spectral_free (gsfs);
    nc_galaxy_shape_pop_data_unref (pop_data);
    nc_galaxy_shape_pop_free (pop);
  }
}

/*
 * Next to |g| = 1, for both conventions: ln P at g and at g / |g|^2 agree, and ln P is
 * finite on both sides down to the last doubles before 1. The two shears of the first
 * check are 1 + d and the double nearest to 1 / (1 + d), whose distance to 1 differs from
 * d / (1 + d) by up to 1.1e-16, a relative difference of up to 1.1e-16 / d that multiplies
 * the logarithmic slope of ln P; the tolerance covers it for the largest slope of the
 * populations of the test.
 */
void
test_nc_galaxy_shape_factor_spectral_near_unit_sides (void)
{
  const gdouble d[]              = {3.0e-6, 1.0e-7, 1.0e-9};
  const gdouble e[]              = {0.178728848695755, 0.3698130249977112};
  const gdouble s                = 0.04034872725605965;
  NcGalaxyShapePopData *pop_data = NULL;
  NcGalaxyShapePop *pop          = _test_pop_new (TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, &pop_data);
  guint j, k;

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);

    for (k = 0; k < G_N_ELEMENTS (d); k++)
    {
      const gdouble g_out   = 1.0 + d[k];
      const gdouble g_in    = 1.0 / g_out;
      const gdouble lnP_out = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g_out, 0.0, e[0], e[1], s);
      const gdouble lnP_in  = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g_in, 0.0, e[0], e[1], s);

      g_assert_cmpfloat_with_epsilon (lnP_out, lnP_in, 1.0e-9 + 5.0 * 1.2e-16 / d[k]);
    }

    for (k = 0; k < 2; k++)
    {
      const gdouble g   = (k == 0) ? nextafter (1.0, 0.0) : nextafter (1.0, 2.0);
      const gdouble lnP = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g, 0.0, e[0], e[1], s);

      g_assert_true (gsl_finite (lnP));
    }

    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/* ln P is the same when g and eps_obs are rotated by the same angle */
void
test_nc_galaxy_shape_factor_spectral_rotation (void)
{
  const gdouble angles[]         = {0.3, 1.9, -2.7};
  const gdouble g1               = 0.25;
  const gdouble g2               = -0.1;
  const gdouble e1               = 0.4;
  const gdouble e2               = 0.3;
  const gdouble s                = 0.03;
  NcGalaxyShapePopData *pop_data = NULL;
  NcGalaxyShapePop *pop          = _test_pop_new (TEST_SHAPE_ARB_POP_BETA, 1.92, 1.98, &pop_data);
  guint j, k;

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
    const gdouble ref                 = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g1, g2, e1, e2, s);

    for (k = 0; k < G_N_ELEMENTS (angles); k++)
    {
      const gdouble c   = cos (angles[k]);
      const gdouble sn  = sin (angles[k]);
      const gdouble rot = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data,
                                                                          g1 * c - g2 * sn, g2 * c + g1 * sn,
                                                                          e1 * c - e2 * sn, e2 * c + e1 * sn, s);

      g_assert_cmpfloat_with_epsilon (rot, ref, 1.0e-9);
    }

    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/* The default reltol, its setter, and ln P at reltol 1e-3 and 1e-7 within 1e-6 of each other */
void
test_nc_galaxy_shape_factor_spectral_reltol (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapePopData *pop_data    = NULL;
  NcGalaxyShapePop *pop             = _test_pop_new (TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, &pop_data);
  gdouble loose, tight;

  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_reltol (gsfs), ==, 1.0e-5);

  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-3);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_reltol (gsfs), ==, 1.0e-3);
  loose = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, 0.2, 0.1, 0.4, -0.3, 0.05);

  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-7);
  tight = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, 0.2, 0.1, 0.4, -0.3, 0.05);

  g_assert_cmpfloat_with_epsilon (loose, tight, 1.0e-6);

  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/*
 * The accepted panels cover the disk, ln P returned with them is the one of
 * eval_ln_direct, a rotation of the frame shifts their angles by the rotation angle (the
 * edges follow the peak found by the scan, which the rotation moves by rounding only), and
 * |g| = 1 has no panels.
 */
void
test_nc_galaxy_shape_factor_spectral_panels (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapePopData *pop_data    = NULL;
  NcGalaxyShapePop *pop             = _test_pop_new (TEST_SHAPE_ARB_POP_BETA, 1.57, 1.65, &pop_data);
  const gdouble g1                  = 0.2;
  const gdouble g2                  = 0.0;
  const gdouble e1                  = 0.4;
  const gdouble e2                  = -0.3;
  const gdouble s                   = 0.05;
  const gdouble angle               = 0.8;
  const gdouble c                   = cos (angle);
  const gdouble sn                  = sin (angle);
  gdouble area                      = 0.0;
  gdouble ln_P, ln_P_rot;
  NcmMatrix *panels     = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data, g1, g2, e1, e2, s, &ln_P);
  NcmMatrix *panels_rot = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data,
                                                                       g1 * c - g2 * sn, g2 * c + g1 * sn,
                                                                       e1 * c - e2 * sn, e2 * c + e1 * sn, s, &ln_P_rot);
  guint i;

  g_assert_nonnull (panels);
  g_assert_nonnull (panels_rot);
  g_assert_cmpuint (ncm_matrix_ncols (panels), ==, 8);
  g_assert_cmpuint (ncm_matrix_nrows (panels), ==, ncm_matrix_nrows (panels_rot));
  g_assert_cmpfloat (ln_P, ==, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, g1, g2, e1, e2, s));
  g_assert_cmpfloat_with_epsilon (ln_P_rot, ln_P, 1.0e-9);

  for (i = 0; i < ncm_matrix_nrows (panels); i++)
  {
    const gdouble ra = ncm_matrix_get (panels, i, 0);
    const gdouble rb = ncm_matrix_get (panels, i, 1);
    const gdouble ta = ncm_matrix_get (panels, i, 2);
    const gdouble tb = ncm_matrix_get (panels, i, 3);

    g_assert_cmpfloat (ra, >=, 0.0);
    g_assert_cmpfloat (rb, <=, 1.0);
    area += (rb - ra) * (tb - ta);

    g_assert_cmpfloat_with_epsilon (ncm_matrix_get (panels_rot, i, 0), ra, 1.0e-8);
    g_assert_cmpfloat_with_epsilon (ncm_matrix_get (panels_rot, i, 1), rb, 1.0e-8);
    g_assert_cmpfloat_with_epsilon (ncm_matrix_get (panels_rot, i, 2), ta + angle, 1.0e-8);
    g_assert_cmpfloat_with_epsilon (ncm_matrix_get (panels_rot, i, 3), tb + angle, 1.0e-8);
  }

  g_assert_cmpfloat_with_epsilon (area, 2.0 * M_PI, 1.0e-12);

  ncm_matrix_free (panels);
  ncm_matrix_free (panels_rot);

  panels = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data, 1.0, 0.0, e1, e2, s, &ln_P);
  g_assert_null (panels);
  g_assert_true (gsl_finite (ln_P));

  /* Below 1 - |g| = 1e-5 the panels are those of the integral at that distance */
  {
    NcmMatrix *panels_s = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data, 1.0 - 2.0e-5, 0.0, e1, e2, s, NULL);
    guint n_s;

    panels = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data, 1.0 - 1.0e-7, 0.0, e1, e2, s, &ln_P);
    g_assert_nonnull (panels);
    g_assert_cmpfloat (ln_P, ==, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, 1.0 - 1.0e-7, 0.0, e1, e2, s));
    n_s = ncm_matrix_nrows (panels);
    ncm_matrix_free (panels);

    panels = nc_galaxy_shape_factor_spectral_eval_panels (gsfs, pop, pop_data, 1.0 - 1.0e-9, 0.0, e1, e2, s, NULL);
    g_assert_nonnull (panels);
    g_assert_cmpuint (ncm_matrix_nrows (panels), ==, n_s);
    g_assert_cmpfloat (ncm_matrix_get (panels, 0, 6), !=, ncm_matrix_get (panels_s, 0, 6));

    ncm_matrix_free (panels);
    ncm_matrix_free (panels_s);
  }

  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_shape_pop_data_unref (pop_data);
  nc_galaxy_shape_pop_free (pop);
}

/* A copy made through NcmSerialize has the convention, the reltol and the interpolant settings of the original */
void
test_nc_galaxy_shape_factor_spectral_serialize (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcmSerialize *ser                 = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyShapeFactorSpectral *dup;

  g_assert_false (nc_galaxy_shape_factor_spectral_get_use_interp (gsfs));
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs), ==, 0.3);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs), ==, 1.0e-5);
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (gsfs), ==, 2);

  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-4);
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.2);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-4);
  nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, 4);
  dup = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (ncm_serialize_dup_obj (ser, G_OBJECT (gsfs)));

  g_assert_cmpint (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (dup)), ==, NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_reltol (dup), ==, 1.0e-4);
  g_assert_true (nc_galaxy_shape_factor_spectral_get_use_interp (dup));
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_g_max (dup), ==, 0.2);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_abstol (dup), ==, 1.0e-4);
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (dup), ==, 4);

  ncm_serialize_free (ser);
  nc_galaxy_shape_factor_spectral_free (dup);
  NCM_TEST_FREE (nc_galaxy_shape_factor_spectral_free, gsfs);
}

/*
 * nc_galaxy_shape_factor_eval_ln_marginal() and nc_galaxy_shape_factor_eval_marginal() on
 * the data of one galaxy give the value of eval_ln_direct with the population data and
 * the noise of that galaxy.
 */
void
test_nc_galaxy_shape_factor_spectral_hooks (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble s                    = 0.1;
  const gdouble g1                   = 0.15;
  const gdouble g2                   = -0.05;
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  NcGalaxyShapeFactorData *data;
  gdouble ln_direct;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);

  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s);

  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, ln_direct);
  /* With AVX2 the exponential is ncm_vmath_exp(), within 2 ulp of exp() */
  ncm_assert_cmpdouble_e (nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, exp (ln_direct), 1.0e-15, 0.0);

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/*
 * The interpolant of ln P in g against the integral, for both conventions and three
 * galaxies: one with the center g_c of the interpolant inside |g| <= interp-g-max, one
 * with it outside, and one with |eps_obs| > 1. At shears of the disk the two agree within
 * TEST_SPECTRAL_INTERP_TOL and the value does not change under the reflection of g on
 * the direction of eps_obs; outside the disk the value is the one of the integral.
 */
void
test_nc_galaxy_shape_factor_spectral_interp (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const gdouble gal[][3] = {
    {
      0.3, 0.2, 0.1
    }, {
      -0.55, 0.45, 0.05
    }, {
      1.05, 0.1, 0.2
    }
  };
  const gdouble abs_g[]              = {0.0, 0.05, 0.15, 0.29, 0.3 * (1.0 - 1.0e-12)};
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  guint j, k, l, m;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);

    for (k = 0; k < G_N_ELEMENTS (gal); k++)
    {
      const gdouble e1              = gal[k][0];
      const gdouble e2              = gal[k][1];
      const gdouble s               = gal[k][2];
      const gdouble e               = hypot (e1, e2);
      GPtrArray *data_array         = g_ptr_array_new ();
      NcGalaxyShapeFactorData *data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);

      nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
      g_ptr_array_add (data_array, data);
      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

      for (l = 0; l < G_N_ELEMENTS (abs_g); l++)
      {
        for (m = 0; m < 8; m++)
        {
          const gdouble angle  = 0.1 + 2.0 * M_PI * m / 8.0;
          const gdouble g1     = abs_g[l] * cos (angle);
          const gdouble g2     = abs_g[l] * sin (angle);
          const gdouble dot    = (g1 * e1 + g2 * e2) / (e * e);
          const gdouble g1_ref = 2.0 * dot * e1 - g1;
          const gdouble g2_ref = 2.0 * dot * e2 - g2;
          const gdouble ln_int = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s);
          const gdouble ln_P   = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);

          g_assert_cmpfloat_with_epsilon (ln_P, ln_int, TEST_SPECTRAL_INTERP_TOL);
          g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1_ref, g2_ref, e1, e2), ln_P, 1.0e-11);
          ncm_assert_cmpdouble_e (nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, exp (ln_P), 1.0e-15, 0.0);
        }
      }

      g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, 0.25, 0.2, e1, e2), ==,
                         nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, 0.25, 0.2, e1, e2, s));

      g_ptr_array_unref (data_array);
      nc_galaxy_shape_factor_data_unref (data);
    }

    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * The interpolant of a galaxy follows what it is evaluated for: after a change of the
 * population parameters, of eps_obs, of std_noise, of interp-abstol, of interp-g-max, of
 * reltol and of interp-min-radial-panels the value at one shear agrees with the integral for the new values, and it
 * differs from the value before the change where the change moves ln P by more than the
 * tolerance.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_rebuild (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble g1                   = 0.12;
  const gdouble g2                   = -0.07;
  gdouble e1                         = 0.3;
  gdouble e2                         = 0.2;
  gdouble s                          = 0.1;
  NcGalaxyShapeFactorData *data;
  gdouble ln_P, ln_P_prev;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);

  data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, ln_P);

  /* The population parameters */
  ln_P_prev = ln_P;
  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "beta", 2.2, NULL);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The observed ellipticity */
  ln_P_prev = ln_P;
  e1        = 0.25;
  e2        = -0.3;
  ln_P      = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The noise */
  ln_P_prev = ln_P;
  s         = 0.06;
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The tolerances and the range: a new interpolant, the same function */
  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 1.0e-4);

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-5);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, ln_P_prev);

  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.2);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-7);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, 4);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  /* Without the interpolant the value is the integral */
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, FALSE);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==,
                     nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s));

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/*
 * nc_galaxy_shape_factor_data_prepare() builds the interpolant: its value is the one of a
 * galaxy whose interpolant is built at the first evaluation. Without use-interp it builds
 * nothing and the value is the integral.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_data_prepare (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  GPtrArray *data_array              = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_factor_data_unref);
  const gdouble g1                   = 0.12;
  const gdouble g2                   = -0.07;
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble s                    = 0.1;
  NcGalaxyShapeFactorData *prepared, *lazy, *direct;
  gdouble ln_P;
  guint i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < 3; i++)
  {
    NcGalaxyShapeFactorData *data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);

    nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
  }

  prepared = g_ptr_array_index (data_array, 0);
  lazy     = g_ptr_array_index (data_array, 1);
  direct   = g_ptr_array_index (data_array, 2);

  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  /* Without use-interp nothing is built */
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, direct, 3.0);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, direct, g1, g2, e1, e2), ==,
                     nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, direct->pop_data, g1, g2, e1, e2, s));

  /* With use-interp; the second call finds the interpolant valid */
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, prepared, 3.0);
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, prepared, 3.0);

  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, prepared, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, ==, nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, lazy, g1, g2, e1, e2));
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, prepared->pop_data, g1, g2, e1, e2, s),
                                  TEST_SPECTRAL_INTERP_TOL);

  g_ptr_array_unref (data_array);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/* An interpolation range of 1/3 or more aborts */
void
test_nc_galaxy_shape_factor_spectral_interp_g_max (void)
{
  g_test_trap_subprocess ("/nc/galaxy_shape_factor_spectral/interp/g_max/subprocess", 0, 0);
  g_test_trap_assert_failed ();
}

void
test_nc_galaxy_shape_factor_spectral_interp_g_max_subprocess (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);

  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.34);
}

/* A noise of 1e-200, whose variance is zero in double precision, aborts with the point */
void
test_nc_galaxy_shape_factor_spectral_invalid_noise (void)
{
  g_test_trap_subprocess ("/nc/galaxy_shape_factor_spectral/invalid_noise/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_galaxy_shape_factor_spectral:*NcGalaxyShapePopGauss*std_noise*");
}

void
test_nc_galaxy_shape_factor_spectral_invalid_noise_subprocess (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcGalaxyShapePopData *pop_data    = NULL;
  NcGalaxyShapePop *pop             = _test_pop_new (TEST_SHAPE_ARB_POP_GAUSS, 0.3, 0.0, &pop_data);

  nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pop_data, 0.2, 0.1, 0.4, -0.3, 1.0e-200);
}

/* One shape factor data per galaxy of test_spectral_obs_gal, with its prepared fields */
static GPtrArray *
_test_obs_data_array_new (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactor *gsf)
{
  GPtrArray *data_array = g_ptr_array_new_with_free_func ((GDestroyNotify) nc_galaxy_shape_factor_data_unref);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_spectral_obs_gal); i++)
  {
    NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
    NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
    NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);

    p_data->ra  = 0.3;
    p_data->dec = 0.0;
    nc_galaxy_shape_factor_data_set (gsf, data, test_spectral_obs_gal[i][0], test_spectral_obs_gal[i][1], test_spectral_obs_gal[i][2],
                                     0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);

    nc_galaxy_position_factor_data_unref (p_data);
    nc_galaxy_redshift_factor_data_unref (z_data);
  }

  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  return data_array;
}

/*
 * A catalog of the galaxies of test_spectral_obs_gal prepared by gsfs (use-interp set) for
 * the population of test, after a serialization round trip; ln_P receives the value of each
 * galaxy at g.
 */
static NcGalaxyWLObs *
_test_obs_prepared_new (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactorSpectral *gsfs, const gdouble g1, const gdouble g2, gdouble *ln_P)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (gsfs);
  GPtrArray *data_array    = _test_obs_data_array_new (test, gsf);
  GList *columns           = nc_galaxy_shape_factor_data_required_columns (g_ptr_array_index (data_array, 0));
  GStrv col_names          = g_new0 (gchar *, g_list_length (columns) + 1);
  NcmSerialize *ser        = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyWLObs *obs, *obs_dup;
  GVariant *obs_ser;
  GList *l;
  guint i = 0;

  for (l = columns; l != NULL; l = l->next)
    col_names[i++] = g_strdup (l->data);

  obs = nc_galaxy_wl_obs_new (nc_galaxy_shape_factor_get_ellip_conv (gsf), NC_WL_ELLIPTICITY_FRAME_CELESTIAL, data_array->len, col_names);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

    nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 3.0);
    ln_P[i] = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2);
    nc_galaxy_shape_factor_data_write_row (data, obs, i);
  }

  nc_galaxy_shape_factor_spectral_set_obs_source (gsfs, test->pop, obs);

  obs_ser = ncm_serialize_to_variant (ser, G_OBJECT (obs));
  ncm_serialize_reset (ser, FALSE);
  obs_dup = NC_GALAXY_WL_OBS (ncm_serialize_from_variant (ser, obs_ser));

  g_variant_unref (obs_ser);
  nc_galaxy_wl_obs_free (obs);
  ncm_serialize_free (ser);
  g_strfreev (col_names);
  g_list_free_full (columns, g_free);
  g_ptr_array_unref (data_array);

  return obs_dup;
}

/* The rows of obs read into new data of gsf */
static GPtrArray *
_test_obs_read (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactor *gsf, NcGalaxyWLObs *obs)
{
  GPtrArray *data_array = _test_obs_data_array_new (test, gsf);
  guint i;

  for (i = 0; i < data_array->len; i++)
    nc_galaxy_shape_factor_data_read_row (g_ptr_array_index (data_array, i), obs, i);

  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  return data_array;
}

/*
 * A prepared catalog gives, after a serialization round trip, the values of the interpolants
 * it was prepared with; they come from the stored interpolants (a catalog whose row 0 holds
 * the interpolant of galaxy 1 gives another value). A catalog changed after the preparation
 * carries no interpolants, and the galaxies build theirs again.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_obs (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs_w = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs_r = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf_r          = NC_GALAXY_SHAPE_FACTOR (gsfs_r);
  const gdouble g1                    = 0.12;
  const gdouble g2                    = -0.07;
  gdouble ln_P[G_N_ELEMENTS (test_spectral_obs_gal)];
  NcGalaxyWLObs *obs;
  GPtrArray *data_array;
  guint i;

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs_w, TRUE);
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs_r, TRUE);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  obs = _test_obs_prepared_new (test, gsfs_w, g1, g2, ln_P);

  for (i = 0; i < G_N_ELEMENTS (test_spectral_obs_gal); i++)
    g_assert_true (NCM_IS_SPECTRAL_POLAR_PANELS (nc_galaxy_wl_obs_peek_shape_interp (obs, i)));

  data_array = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2), ==, ln_P[i]);
  }

  g_ptr_array_unref (data_array);

  /* Row 0 with the interpolant of galaxy 1 */
  nc_galaxy_wl_obs_set_shape_interp (obs, 0, nc_galaxy_wl_obs_peek_shape_interp (obs, 1));
  data_array = _test_obs_read (test, gsf_r, obs);
  nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, g_ptr_array_index (data_array, 0), 3.0);
  g_assert_cmpfloat (fabs (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, g_ptr_array_index (data_array, 0), g1, g2,
                                                                    test_spectral_obs_gal[0][0], test_spectral_obs_gal[0][1]) - ln_P[0]), >, 1.0e-2);
  g_ptr_array_unref (data_array);

  /* A change of a row: no interpolant is loaded, each galaxy builds its own */
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_M, 2, 0.0, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs, 0));
  data_array = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2), ==, ln_P[i]);
  }

  g_ptr_array_unref (data_array);
  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_shape_factor_spectral_free (gsfs_w);
  nc_galaxy_shape_factor_spectral_free (gsfs_r);
}

/* Every change of the fit after a catalog is prepared aborts with what differs */
void
test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch (void)
{
  guint k;

  for (k = 0; k < G_N_ELEMENTS (test_spectral_obs_mismatch_cases); k++)
  {
    gchar *path = g_strdup_printf ("/nc/galaxy_shape_factor_spectral/interp/obs/mismatch/%s/subprocess", test_spectral_obs_mismatch_cases[k][0]);

    g_test_trap_subprocess (path, 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (test_spectral_obs_mismatch_cases[k][1]);
    g_free (path);
  }
}

void
test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch_subprocess (gconstpointer pdata)
{
  const gchar *mismatch = pdata;
  TestNcGalaxyShapeFactorSpectral test;
  NcGalaxyShapeFactorSpectral *gsfs_w = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs_r;
  gdouble ln_P[G_N_ELEMENTS (test_spectral_obs_gal)];
  NcGalaxyWLObs *obs;
  GPtrArray *data_array;
  NcGalaxyShapeFactorData *data;

  test_nc_galaxy_shape_factor_spectral_new (&test, NULL);
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs_w, TRUE);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  obs = _test_obs_prepared_new (&test, gsfs_w, 0.12, -0.07, ln_P);

  gsfs_r = nc_galaxy_shape_factor_spectral_new (g_str_equal (mismatch, "conv") ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs_r, !g_str_equal (mismatch, "use_interp"));
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);

  if (g_str_equal (mismatch, "pop"))
    ncm_model_param_set_by_name (NCM_MODEL (test.pop), "beta", 1.65 * (1.0 + 1.0e-12), NULL);
  else if (g_str_equal (mismatch, "free"))
    ncm_model_param_set_ftype (NCM_MODEL (test.pop), 0, NCM_PARAM_TYPE_FREE);
  else if (g_str_equal (mismatch, "abstol"))
    nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, 1.0e-4);

  data_array = _test_obs_read (&test, NC_GALAXY_SHAPE_FACTOR (gsfs_r), obs);
  data       = g_ptr_array_index (data_array, 0);
  nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test.mset, data, 3.0);

  /* Passed the check: a change of the population afterwards is found by the evaluation */
  ncm_model_param_set_by_name (NCM_MODEL (test.pop), "beta", 2.0, NULL);
  nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test.pop, data, 0.12, -0.07, data->epsilon_obs_1, data->epsilon_obs_2);
}

/* nc_galaxy_shape_factor_eval_ln_marginal_array() gives, point by point, what
 * nc_galaxy_shape_factor_eval_ln_marginal() gives, by integration and from the interpolant,
 * for shears inside and outside interp-g-max and for two observed ellipticities. */
void
test_nc_galaxy_shape_factor_spectral_eval_array (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const guint n                      = 24;
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcmVector *g1                      = ncm_vector_new (n);
  NcmVector *g2                      = ncm_vector_new (n);
  NcmVector *e1                      = ncm_vector_new (n);
  NcmVector *e2                      = ncm_vector_new (n);
  NcmVector *out                     = ncm_vector_new (n);
  guint j, i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < n; i++)
  {
    const gdouble abs_g = 0.4 * i / (n - 1.0);
    const gdouble angle = 0.3 + 2.0 * M_PI * i / 7.0;

    ncm_vector_set (g1, i, abs_g * cos (angle));
    ncm_vector_set (g2, i, abs_g * sin (angle));
    ncm_vector_set (e1, i, (i < n / 2) ? 0.3 : -0.2);
    ncm_vector_set (e2, i, (i < n / 2) ? 0.2 : 0.4);
  }

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, j == 1);
    nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    nc_galaxy_shape_factor_eval_ln_marginal_array (gsf, test->pop, data, g1, g2, e1, e2, out);

    for (i = 0; i < n; i++)
    {
      const gdouble ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data,
                                                                    ncm_vector_get (g1, i), ncm_vector_get (g2, i),
                                                                    ncm_vector_get (e1, i), ncm_vector_get (e2, i));

      g_assert_true (gsl_finite (ln_P));
      g_assert_cmpfloat (ncm_vector_get (out, i), ==, ln_P);
    }

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  ncm_vector_free (g1);
  ncm_vector_free (g2);
  ncm_vector_free (e1);
  ncm_vector_free (e2);
  ncm_vector_free (out);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/* The critical-density hash follows the cosmology and the lens redshift but not the halo
 * center: moving ra and dec leaves it unchanged, moving z changes it, and setting z back
 * restores it. */
void
test_nc_galaxy_shape_factor_spectral_crit_hash (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE));
  NcmModel *hp             = ncm_mset_peek (test->mset, nc_halo_position_id ());
  NcmModel *cosmo          = ncm_mset_peek (test->mset, nc_hicosmo_id ());
  const gdouble z          = ncm_model_param_get (hp, NC_HALO_POSITION_Z);
  guint64 h0, h_center, h_z, h_back;

  nc_galaxy_shape_factor_prepare (gsf, test->mset);
  h0 = nc_galaxy_shape_factor_get_crit_hash (gsf);

  ncm_model_param_set (hp, NC_HALO_POSITION_RA, ncm_model_param_get (hp, NC_HALO_POSITION_RA) + 0.01);
  ncm_model_param_set (hp, NC_HALO_POSITION_DEC, ncm_model_param_get (hp, NC_HALO_POSITION_DEC) - 0.01);
  nc_galaxy_shape_factor_prepare (gsf, test->mset);
  h_center = nc_galaxy_shape_factor_get_crit_hash (gsf);
  g_assert_cmpuint (h_center, ==, h0);

  ncm_model_param_set (hp, NC_HALO_POSITION_Z, z + 0.05);
  nc_galaxy_shape_factor_prepare (gsf, test->mset);
  h_z = nc_galaxy_shape_factor_get_crit_hash (gsf);
  g_assert_cmpuint (h_z, !=, h0);

  ncm_model_param_set (hp, NC_HALO_POSITION_Z, z);
  nc_galaxy_shape_factor_prepare (gsf, test->mset);
  h_back = nc_galaxy_shape_factor_get_crit_hash (gsf);
  g_assert_cmpuint (h_back, ==, h0);

  ncm_model_param_set (cosmo, 0, ncm_model_param_get (cosmo, 0) * 1.01);
  nc_galaxy_shape_factor_prepare (gsf, test->mset);
  g_assert_cmpuint (nc_galaxy_shape_factor_get_crit_hash (gsf), !=, h0);

  nc_galaxy_shape_factor_free (gsf);
}

/* nc_galaxy_shape_factor_eval_nodes() gives, shear by shear, what the per-point functions
 * give, within rounding (TEST_SPECTRAL_ASSERT_LOOP), by integration and from the interpolant, inside and
 * outside interp-g-max, in the linear and the logarithmic forms. */
void
test_nc_galaxy_shape_factor_spectral_eval_nodes (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const guint n                      = 40;
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  gdouble g1[40], g2[40], out[40];
  guint i, use_interp, want_log;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < n; i++)
  {
    const gdouble abs_g = 0.45 * i / (n - 1.0);
    const gdouble angle = 0.2 + 2.0 * M_PI * i / 9.0;

    g1[i] = abs_g * cos (angle);
    g2[i] = abs_g * sin (angle);
  }

  for (use_interp = 0; use_interp < 2; use_interp++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, use_interp);
    nc_galaxy_shape_factor_data_set (gsf, data, -0.25, 0.35, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    for (want_log = 0; want_log < 2; want_log++)
    {
      nc_galaxy_shape_factor_eval_nodes (gsf, test->pop, data, n, g1, g2, -0.25, 0.35, want_log, out);

      for (i = 0; i < n; i++)
      {
        const gdouble v = want_log ?
                          nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1[i], g2[i], -0.25, 0.35) :
                          nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g1[i], g2[i], -0.25, 0.35);

        g_assert_true (gsl_finite (v));

        /* ln P within TEST_SPECTRAL_ASSERT_LOOP, and P relative within the same, since the
         * relative error of P is the error of ln P */
        if (want_log)
          TEST_SPECTRAL_ASSERT_LOOP (out[i], v);
        else
          g_assert_cmpfloat (fabs (out[i] / v - 1.0), <=, 1.0e-13 * MAX (1.0, fabs (log (v))));
      }
    }

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

#define TEST_LINE_N 60

/*
 * One case of eval_line: n positions on [0, t_max] in a scrambled order, on the line
 * b + t u. With exact, every value is that of eval_nodes at the same shear (to the rounding
 * of the shear); otherwise ln P is within tol of the integral at each shear (P within tol
 * relative).
 */
static void
_test_eval_line_case (NcGalaxyShapeFactor *gsf, NcGalaxyShapePop *pop, NcGalaxyShapeFactorData *data, const guint n, const gdouble t_max,
                      const gboolean exact, const gdouble tol)
{
  const gdouble u_1 = 1.05 * cos (0.7);
  const gdouble u_2 = 1.05 * sin (0.7);
  const gdouble b_1 = 0.004;
  const gdouble b_2 = -0.003;
  const gdouble e_1 = -0.25;
  const gdouble e_2 = 0.35;
  gdouble t[TEST_LINE_N], out[TEST_LINE_N], g1[TEST_LINE_N], g2[TEST_LINE_N], ref[TEST_LINE_N];
  guint want_log, j;

  g_assert_cmpuint (n, <=, TEST_LINE_N);

  for (j = 0; j < n; j++)
  {
    /* 37 is prime to 60: a permutation of the n equally spaced points */
    t[j]  = t_max * ((37 * j) % n) / (n - 1.0);
    g1[j] = b_1 + t[j] * u_1;
    g2[j] = b_2 + t[j] * u_2;
  }

  for (want_log = 0; want_log < 2; want_log++)
  {
    nc_galaxy_shape_factor_eval_line (gsf, pop, data, n, t, u_1, u_2, b_1, b_2, e_1, e_2, want_log, out);
    nc_galaxy_shape_factor_eval_nodes (gsf, pop, data, n, g1, g2, e_1, e_2, want_log, ref);

    for (j = 0; j < n; j++)
    {
      g_assert_true (gsl_finite (out[j]));

      if (exact) /* the shears formed here may differ by an ulp from those of the library */
        ncm_assert_cmpdouble_e (out[j], ==, ref[j], 1.0e-13, 0.0);
      else if (want_log)
        g_assert_cmpfloat (fabs (out[j] - ref[j]), <=, tol);
      else
        ncm_assert_cmpdouble_e (out[j], ==, ref[j], tol, 0.0);
    }
  }

  /* In place: out is t */
  {
    gdouble tt[TEST_LINE_N];

    memcpy (tt, t, sizeof (gdouble) * n);
    nc_galaxy_shape_factor_eval_line (gsf, pop, data, n, tt, u_1, u_2, b_1, b_2, e_1, e_2, TRUE, tt);
    nc_galaxy_shape_factor_eval_line (gsf, pop, data, n, t, u_1, u_2, b_1, b_2, e_1, e_2, TRUE, out);

    for (j = 0; j < n; j++)
      g_assert_cmpfloat (tt[j], ==, out[j]);
  }
}

/*
 * nc_galaxy_shape_factor_eval_line(): the table along the line agrees with the integral at
 * each shear within line-abstol, for both conventions and both populations; without it
 * (line-abstol zero, the interpolant, a line too short for a table, or one that reaches
 * |g| > 0.9) every value equals that of eval_nodes.
 */
void
test_nc_galaxy_shape_factor_spectral_eval_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data    = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data    = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapePopData *gauss_data      = NULL;
  NcGalaxyShapePop *gauss               = _test_pop_new (TEST_SHAPE_ARB_POP_GAUSS, 0.25, 0.0, &gauss_data);
  const NcGalaxyWLObsEllipConv convs[2] = {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET};
  guint ic;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (ic = 0; ic < 2; ic++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (convs[ic]);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    NcGalaxyShapePop *pops[2]         = {test->pop, gauss};
    guint ip;

    /* The default: reltol / 10; the cases below use 1e-8 */
    g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_line_abstol (gsfs), <, 0.0);
    nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 1.0e-8);
    nc_galaxy_shape_factor_data_set (gsf, data, -0.25, 0.35, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);

    for (ip = 0; ip < 2; ip++)
    {
      NcGalaxyShapePopData *pop_data = (ip == 0) ? data->pop_data : gauss_data;

      /* The population of the galaxy: eval_line passes data's pop_data to the integrals */
      if (ip == 1)
      {
        NcGalaxyShapePopData *own = data->pop_data;

        data->pop_data = pop_data;
        _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.08, FALSE, 1.0e-8);
        _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.3, FALSE, 1.0e-8);
        data->pop_data = own;
        continue;
      }

      nc_galaxy_shape_pop_prepare (pops[ip], pop_data);

      /* The table */
      _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.08, FALSE, 1.0e-8);
      _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.3, FALSE, 1.0e-8);

      /* Too short for a table: 9 nodes cannot hold two grids of 5 points */
      _test_eval_line_case (gsf, pops[ip], data, 9, 0.3, TRUE, 0.0);

      /* The shortest line with a table: 10 nodes, grids of 3 and 5 points */
      _test_eval_line_case (gsf, pops[ip], data, 10, 0.08, FALSE, 1.0e-8);

      /* A line that leaves |g| <= 0.9 */
      _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.92, TRUE, 0.0);

      /* No table */
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 0.0);
      _test_eval_line_case (gsf, pops[ip], data, TEST_LINE_N, 0.3, TRUE, 0.0);
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 1.0e-8);
    }

    /*
     * A table that cannot reach its tolerance: every node is integrated, with the values of
     * the integrals, at each evaluation alike.
     */
    {
      NcGalaxyShapeFactorData *data2 = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);

      nc_galaxy_shape_factor_data_set (gsf, data2, -0.25, 0.35, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
      nc_galaxy_shape_pop_prepare (test->pop, data2->pop_data);
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 1.0e-16);
      _test_eval_line_case (gsf, test->pop, data2, TEST_LINE_N, 0.3, TRUE, 0.0);
      nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 1.0e-8);
      nc_galaxy_shape_factor_data_unref (data2);
    }

    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  /* With the interpolant the default eval_line forms the shears for eval_nodes */
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
    nc_galaxy_shape_factor_data_set (gsf, data, -0.25, 0.35, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    _test_eval_line_case (gsf, test->pop, data, TEST_LINE_N, 0.3, TRUE, 0.0);

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  nc_galaxy_shape_pop_data_unref (gauss_data);
  nc_galaxy_shape_pop_free (gauss);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * nc_galaxy_shape_factor_line_new() and _eval(): the values of the linear shape integrand
 * at each redshift, from the table within line-abstol over [z_lo, z_hi] and exactly outside
 * it (a redshift beyond z_hi) and without the table (line-abstol zero, the interpolant).
 */
void
test_nc_galaxy_shape_factor_spectral_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data    = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data    = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  const NcGalaxyWLObsEllipConv convs[2] = {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET};
  const gdouble z_lo                    = 0.25;
  const gdouble z_hi                    = 3.0;
  guint ic, mode;

  p_data->ra  = 0.05;
  p_data->dec = 0.02;

  for (ic = 0; ic < 2; ic++)
  {
    /* mode 0: the table; 1: line-abstol zero; 2: the interpolant */
    for (mode = 0; mode < 3; mode++)
    {
      NcGalaxyShapeFactorSpectral *gsfs   = nc_galaxy_shape_factor_spectral_new (convs[ic]);
      NcGalaxyShapeFactor *gsf            = NC_GALAXY_SHAPE_FACTOR (gsfs);
      NcGalaxyShapeFactorData *data       = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
      GPtrArray *data_array               = g_ptr_array_new ();
      NcGalaxyShapeFactorIntegrand *integ = NULL;
      NcGalaxyShapeFactorLine *line;
      guint k;

      if (mode == 1)
        nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 0.0);

      nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, mode == 2);
      nc_galaxy_shape_factor_data_set (gsf, data, -0.25, 0.35, 0.1, 0.004, -0.003, 0.05, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
      g_ptr_array_add (data_array, data);
      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
      integ = nc_galaxy_shape_factor_integ (gsf, test->mset, FALSE);

      line = nc_galaxy_shape_factor_line_new (gsf, data, z_lo, z_hi);
      g_assert_true (nc_galaxy_shape_factor_line_has_table (line) == (mode == 0));

      for (k = 0; k <= 40; k++)
      {
        const gdouble z   = z_lo + (z_hi - z_lo) * k / 40.0;
        const gdouble ref = nc_galaxy_shape_factor_integrand_eval (integ, z, data);
        const gdouble val = nc_galaxy_shape_factor_line_eval (line, z);

        g_assert_true (gsl_finite (val) && (val > 0.0));

        if (mode == 0)
          ncm_assert_cmpdouble_e (val, ==, ref, 1.0e-8, 0.0);
        else
          ncm_assert_cmpdouble_e (val, ==, ref, 1.0e-14, 0.0);
      }

      /* Beyond the range of the table the integral is computed */
      ncm_assert_cmpdouble_e (nc_galaxy_shape_factor_line_eval (line, 4.0), ==, nc_galaxy_shape_factor_integrand_eval (integ, 4.0, data), 1.0e-14, 0.0);

      nc_galaxy_shape_factor_line_free (line);
      nc_galaxy_shape_factor_integrand_free (integ);
      g_ptr_array_unref (data_array);
      nc_galaxy_shape_factor_data_unref (data);
      nc_galaxy_shape_factor_free (gsf);
    }
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * restrict-range: with the mass bounded at 10^14 and the center fixed, a galaxy 0.3 degrees
 * from the halo reaches a shear far below 0.29, so its interpolant stops below it. A shear
 * there is the integral (bit for bit), a small shear is the interpolant (within the
 * tolerance of the integral, and not the integral itself), and eval_nodes agrees with
 * the per-point values within TEST_SPECTRAL_LOOP_ULP. With cDelta free over [1, 10] the same holds. A box of the center
 * that takes the halo to the galaxy (ra and dec free over +-0.4 degrees) gives back the
 * whole disk: 0.29 is interpolated again.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp           = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms                      = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcmModel *hp                       = ncm_mset_peek (test->mset, nc_halo_position_id ());
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble s                    = 0.1;
  guint c_free;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  ncm_model_param_set_upper_bound (hms, NC_HALO_CM_PARAM_LOG10M_DELTA, 14.0);

  for (c_free = 0; c_free < 3; c_free++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();
    const gboolean centre_box         = (c_free == 2);
    const gdouble g_far[2]            = {0.29 * cos (0.4), 0.29 * sin (0.4)};
    const gdouble g_near[2]           = {0.004 * cos (1.3), 0.004 * sin (1.3)};
    gdouble g1[2], g2[2], out[2];
    gdouble ln_far, ln_near, ln_far_int, ln_near_int;

    if (c_free == 1)
    {
      ncm_model_param_set_lower_bound (hms, NC_HALO_CM_PARAM_C_DELTA, 1.0);
      ncm_model_param_set_upper_bound (hms, NC_HALO_CM_PARAM_C_DELTA, 10.0);
      ncm_model_param_set_ftype (hms, NC_HALO_CM_PARAM_C_DELTA, NCM_PARAM_TYPE_FREE);
    }

    if (centre_box)
    {
      ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_RA, -0.4);
      ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_RA, 0.4);
      ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_DEC, -0.4);
      ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_DEC, 0.4);
      ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
      ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);
    }

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
    nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, TRUE);
    g_assert_true (nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs));

    nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
    /* The bound of the galaxy: the top of the redshift support of the fixture is 20 */
    nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);

    ln_far      = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_far[0], g_far[1], e1, e2);
    ln_near     = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_near[0], g_near[1], e1, e2);
    ln_far_int  = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g_far[0], g_far[1], e1, e2, s);
    ln_near_int = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g_near[0], g_near[1], e1, e2, s);

    g_assert_cmpfloat_with_epsilon (ln_near, ln_near_int, TEST_SPECTRAL_INTERP_TOL);
    g_assert_cmpfloat (ln_near, !=, ln_near_int);
    g_assert_cmpfloat_with_epsilon (ln_far, ln_far_int, TEST_SPECTRAL_INTERP_TOL);

    if (centre_box)
      g_assert_cmpfloat (ln_far, !=, ln_far_int);
    else
      g_assert_cmpfloat (ln_far, ==, ln_far_int);

    g1[0] = g_far[0];
    g2[0] = g_far[1];
    g1[1] = g_near[0];
    g2[1] = g_near[1];
    nc_galaxy_shape_factor_eval_nodes (gsf, test->pop, data, 2, g1, g2, e1, e2, TRUE, out);
    TEST_SPECTRAL_ASSERT_LOOP (out[0], ln_far);
    TEST_SPECTRAL_ASSERT_LOOP (out[1], ln_near);

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  ncm_model_param_set_ftype (hms, NC_HALO_CM_PARAM_C_DELTA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/* A one-row catalog of the galaxy of data, written with the interpolant it holds */
static NcGalaxyWLObs *
_test_restrict_obs_new (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapeFactorData *data)
{
  GList *columns  = nc_galaxy_shape_factor_data_required_columns (data);
  GStrv col_names = g_new0 (gchar *, g_list_length (columns) + 1);
  NcGalaxyWLObs *obs;
  GList *l;
  guint i = 0;

  for (l = columns; l != NULL; l = l->next)
    col_names[i++] = g_strdup (l->data);

  obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, 1, col_names);
  nc_galaxy_shape_factor_data_write_row (data, obs, 0);
  nc_galaxy_shape_factor_spectral_set_obs_source (gsfs, test->pop, obs);

  g_strfreev (col_names);
  g_list_free_full (columns, g_free);

  return obs;
}

/* The outer radius rho_max of the interpolant the galaxy of data holds, read from its catalog row */
static gdouble
_test_restrict_rho_max (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapeFactorData *data)
{
  NcGalaxyWLObs *obs = _test_restrict_obs_new (test, gsfs, data);
  gdouble ra, rb, ta, tb;

  g_assert_true (NCM_IS_SPECTRAL_POLAR_PANELS (nc_galaxy_wl_obs_peek_shape_interp (obs, 0)));
  ncm_spectral_polar_panels_get_range (nc_galaxy_wl_obs_peek_shape_interp (obs, 0), &ra, &rb, &ta, &tb);
  nc_galaxy_wl_obs_free (obs);

  return rb;
}

/*
 * restrict-range: the bound depends on the box only. With the center free over +-0.05
 * degrees, moving it inside the box leaves the interpolant of a galaxy 0.3 degrees away with
 * the same radius, below interp-g-max. With the redshift of the halo free the bound does not
 * hold and the interpolant covers interp-g-max.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_stable (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp           = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms                      = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcmModel *hp                       = ncm_mset_peek (test->mset, nc_halo_position_id ());
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble e                    = hypot (e1, e2);
  const gdouble g_c                  = e / (1.0 + sqrt (1.0 - e * e));
  const gdouble ra0                  = ncm_model_param_get (hp, NC_HALO_POSITION_RA);
  const gdouble dec0                 = ncm_model_param_get (hp, NC_HALO_POSITION_DEC);
  gdouble rho_0, rho_moved, rho_z;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  ncm_model_param_set_upper_bound (hms, NC_HALO_CM_PARAM_LOG10M_DELTA, 14.0);
  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_RA, ra0 - 0.05);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_RA, ra0 + 0.05);
  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_DEC, dec0 - 0.05);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_DEC, dec0 + 0.05);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, TRUE);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);

  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  rho_0 = _test_restrict_rho_max (test, gsfs, data);
  g_assert_cmpfloat (rho_0 - g_c, <, nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs));

  ncm_model_param_set (hp, NC_HALO_POSITION_RA, ra0 + 0.03);
  ncm_model_param_set (hp, NC_HALO_POSITION_DEC, dec0 - 0.02);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  rho_moved = _test_restrict_rho_max (test, gsfs, data);
  g_assert_cmpfloat (rho_moved, ==, rho_0);

  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_Z, NCM_PARAM_TYPE_FREE);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  rho_z = _test_restrict_rho_max (test, gsfs, data);
  ncm_assert_cmpdouble_e (rho_z - g_c, ==, nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs), 1.0e-14, 0.0);

  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_Z, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set (hp, NC_HALO_POSITION_RA, ra0);
  ncm_model_param_set (hp, NC_HALO_POSITION_DEC, dec0);

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * restrict-range: a catalog prepared with the center fixed, read by a fit whose box of the
 * center takes the halo to the galaxy, does not cover the shears the galaxy reaches: one
 * warning, the galaxy is counted, and a shear beyond the catalog is the integral.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_short (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp           = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms                      = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcmModel *hp                       = ncm_mset_peek (test->mset, nc_halo_position_id ());
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs2 = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, p_data, z_data);
  NcGalaxyShapeFactorData *data2     = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  GPtrArray *data_array2             = g_ptr_array_new ();
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble g_far[2]             = {0.29 * cos (0.4), 0.29 * sin (0.4)};
  const gdouble ra0                  = ncm_model_param_get (hp, NC_HALO_POSITION_RA);
  const gdouble dec0                 = ncm_model_param_get (hp, NC_HALO_POSITION_DEC);
  NcGalaxyWLObs *obs;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  ncm_model_param_set_upper_bound (hms, NC_HALO_CM_PARAM_LOG10M_DELTA, 14.0);

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, TRUE);
  nc_galaxy_shape_factor_data_set (NC_GALAXY_SHAPE_FACTOR (gsfs), data, e1, e2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, data, 20.0);
  obs = _test_restrict_obs_new (test, gsfs, data);

  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_RA, ra0 - 0.4);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_RA, ra0 + 0.4);
  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_DEC, dec0 - 0.4);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_DEC, dec0 + 0.4);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs2, TRUE);
  nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs2, TRUE);
  nc_galaxy_shape_factor_data_read_row (data2, obs, 0);
  g_ptr_array_add (data_array2, data2);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, data_array2, TRUE, TRUE));

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*do not cover the shears some galaxies reach in the box of this fit*");
  nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, data2, 20.0);
  g_test_assert_expected_messages ();
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_n_loaded_short (gsfs2), ==, 1);

  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->pop, data2, g_far[0], g_far[1], e1, e2),
                     ==,
                     nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs2, test->pop, data2->pop_data, g_far[0], g_far[1], e1, e2, 0.1));

  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);

  nc_galaxy_wl_obs_free (obs);
  g_ptr_array_unref (data_array);
  g_ptr_array_unref (data_array2);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_data_unref (data2);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_shape_factor_spectral_free (gsfs2);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * restrict-range: a catalog read by a fit with the box it was prepared for is never short.
 * The range of a loaded interpolant is (g_c + g_max) - g_c. With the smallest reach, 0.02,
 * and |eps_obs| so close to one that g_c + 0.02 passes one, the sum is rounded to the ulp
 * of [1, 2), well above the relative tolerance of the reach: without the allowance of a few
 * ulps, 10 of these 16 galaxies would be counted short of their own reach.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_loaded_same_box (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp           = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms                      = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  guint n_short                      = 0;
  guint i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  /* A reach below REACH_MIN: the interpolants are built over |g| <= 0.02 */
  ncm_model_param_set_upper_bound (hms, NC_HALO_CM_PARAM_LOG10M_DELTA, 13.0);

  for (i = 0; i < 16; i++)
  {
    NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactorSpectral *gsfs2 = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, p_data, z_data);
    NcGalaxyShapeFactorData *data2     = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, p_data, z_data);
    GPtrArray *data_array              = g_ptr_array_new ();
    GPtrArray *data_array2             = g_ptr_array_new ();
    const gdouble e                    = 0.9999 + 0.00009 * i / 15.0;
    NcGalaxyWLObs *obs;

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
    nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, TRUE);
    nc_galaxy_shape_factor_data_set (NC_GALAXY_SHAPE_FACTOR (gsfs), data, e, 0.0, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, data_array, TRUE, TRUE));
    nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, data, 20.0);
    obs = _test_restrict_obs_new (test, gsfs, data);

    nc_galaxy_shape_factor_spectral_set_use_interp (gsfs2, TRUE);
    nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs2, TRUE);
    nc_galaxy_shape_factor_data_read_row (data2, obs, 0);
    g_ptr_array_add (data_array2, data2);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, data_array2, TRUE, TRUE));
    nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs2), test->mset, data2, 20.0);
    n_short += nc_galaxy_shape_factor_spectral_get_n_loaded_short (gsfs2);

    nc_galaxy_wl_obs_free (obs);
    g_ptr_array_unref (data_array);
    g_ptr_array_unref (data_array2);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_data_unref (data2);
    nc_galaxy_shape_factor_spectral_free (gsfs);
    nc_galaxy_shape_factor_spectral_free (gsfs2);
  }

  g_assert_cmpuint (n_short, ==, 0);

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * An interpolant over a region narrower than std_noise / 8 (interp-g-max 0.025 at
 * std_noise 0.3, with g_c > interp-g-max) is built, and agrees with the integral.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_small_g_max (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble e1                   = 0.5;
  const gdouble e2                   = 0.2;
  const gdouble s                    = 0.3;
  guint m;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.025);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  for (m = 0; m < 8; m++)
  {
    const gdouble angle = 0.1 + 2.0 * M_PI * m / 8.0;
    const gdouble g1    = 0.02 * cos (angle);
    const gdouble g2    = 0.02 * sin (angle);

    g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2),
                                    nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s),
                                    TEST_SPECTRAL_INTERP_TOL);
  }

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * pop-order: with the free parameters of the population in a box, ln P at parameters inside
 * the box is that of the integral within the error of the expansion (Beta, alpha and beta
 * free, 7 nodes each; Gauss, sigma free), for both conventions; eval_nodes gives the values
 * of the per-point evaluation within TEST_SPECTRAL_LOOP_ULP; and parameters outside the box
 * are integrated, bit for bit.
 */
void
test_nc_galaxy_shape_factor_spectral_pop_order (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble s                    = 0.1;
  const gdouble g1[3]                = {0.004, -0.08, 0.17};
  const gdouble g2[3]                = {0.001, 0.05, -0.12};
  guint conv, kind;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (kind = 0; kind < 2; kind++)
  {
    NcGalaxyShapePop *pop = (kind == 0) ? NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ()) : NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
    const gchar *names[2] = {(kind == 0) ? "alpha" : "sigma", (kind == 0) ? "beta" : NULL};
    const gdouble lb[2]   = {(kind == 0) ? 1.45 : 0.24, 1.55};
    const gdouble ub[2]   = {(kind == 0) ? 1.65 : 0.26, 1.75};
    guint p;

    for (p = 0; p < 2; p++)
    {
      guint idx;

      if (names[p] == NULL)
        continue;

      g_assert_true (ncm_model_param_index_from_name (NCM_MODEL (pop), names[p], &idx, NULL));
      ncm_model_param_set_lower_bound (NCM_MODEL (pop), idx, lb[p]);
      ncm_model_param_set_upper_bound (NCM_MODEL (pop), idx, ub[p]);
      ncm_model_param_set (NCM_MODEL (pop), idx, 0.5 * (lb[p] + ub[p]));
      ncm_model_param_set_ftype (NCM_MODEL (pop), idx, NCM_PARAM_TYPE_FREE);
    }

    ncm_mset_set (test->mset, NCM_MODEL (pop), NULL);

    for (conv = 0; conv < 2; conv++)
    {
      NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((conv == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
      NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
      NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
      GPtrArray *data_array             = g_ptr_array_new ();
      guint step, j;

      nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
      nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 7);
      g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_order (gsfs), ==, 7);

      nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
      g_ptr_array_add (data_array, data);
      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
      nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
      g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (gsfs), ==, 0);

      for (step = 0; step < 3; step++)
      {
        /* Two points inside the box, then one outside */
        const gdouble u = (step == 0) ? 0.17 : ((step == 1) ? 0.83 : 1.3);
        NcGalaxyShapePopData *pd;
        gdouble out[3];

        for (p = 0; p < 2; p++)
        {
          if (names[p] != NULL)
            ncm_model_param_set_by_name (NCM_MODEL (pop), names[p], lb[p] + u * (ub[p] - lb[p]), NULL);
        }

        g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
        pd = nc_galaxy_shape_pop_data_new (pop);
        nc_galaxy_shape_pop_prepare (pop, pd);

        nc_galaxy_shape_factor_eval_nodes (gsf, pop, data, 3, g1, g2, e1, e2, TRUE, out);

        for (j = 0; j < 3; j++)
        {
          const gdouble ln_P   = nc_galaxy_shape_factor_eval_ln_marginal (gsf, pop, data, g1[j], g2[j], e1, e2);
          const gdouble ln_int = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, pd, g1[j], g2[j], e1, e2, s);

          TEST_SPECTRAL_ASSERT_LOOP (out[j], ln_P);

          if (step < 2)
            g_assert_cmpfloat_with_epsilon (ln_P, ln_int, 1.0e-6);
          else
            g_assert_cmpfloat (ln_P, ==, ln_int);
        }

        nc_galaxy_shape_pop_data_unref (pd);
      }

      g_ptr_array_unref (data_array);
      nc_galaxy_shape_factor_data_unref (data);
      nc_galaxy_shape_factor_spectral_free (gsfs);
    }

    ncm_mset_set (test->mset, NCM_MODEL (test->pop), NULL);
    nc_galaxy_shape_pop_free (pop);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/* A galaxy of data at the box [lb, ub] of alpha of test->pop, free at its center */
static void
_pop_order_free_alpha (TestNcGalaxyShapeFactorSpectral *test, const gdouble lb, const gdouble ub)
{
  guint idx;

  g_assert_true (ncm_model_param_index_from_name (NCM_MODEL (test->pop), "alpha", &idx, NULL));
  ncm_model_param_set_lower_bound (NCM_MODEL (test->pop), idx, lb);
  ncm_model_param_set_upper_bound (NCM_MODEL (test->pop), idx, ub);
  ncm_model_param_set (NCM_MODEL (test->pop), idx, 0.5 * (lb + ub));
  ncm_model_param_set_ftype (NCM_MODEL (test->pop), idx, NCM_PARAM_TYPE_FREE);
}

/* ln P of data at the shears g1, g2 against the integral, to within tol (0: bit for bit) */
static void
_pop_order_assert_integral (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapeFactorData *data,
                            const guint n, const gdouble *g1, const gdouble *g2, const gdouble tol)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapePopData *pd = nc_galaxy_shape_pop_data_new (test->pop);
  gdouble out[8];
  guint j;

  g_assert_cmpuint (n, <=, G_N_ELEMENTS (out));
  nc_galaxy_shape_pop_prepare (test->pop, pd);
  nc_galaxy_shape_factor_eval_nodes (gsf, test->pop, data, n, g1, g2, 0.3, 0.2, TRUE, out);

  for (j = 0; j < n; j++)
  {
    const gdouble ln_P   = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1[j], g2[j], 0.3, 0.2);
    const gdouble ln_int = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, pd, g1[j], g2[j], 0.3, 0.2, 0.1);

    TEST_SPECTRAL_ASSERT_LOOP (out[j], ln_P);

    if (tol == 0.0)
      g_assert_cmpfloat (ln_P, ==, ln_int);
    else
      g_assert_cmpfloat_with_epsilon (ln_P, ln_int, tol);
  }

  nc_galaxy_shape_pop_data_unref (pd);
}

/*
 * An expansion of order 3 over a wide box of alpha does not converge: the galaxy is counted,
 * warned about once, and integrated, bit for bit, inside the box and after a change of alpha.
 */
void
test_nc_galaxy_shape_factor_spectral_pop_order_diverge (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble g1[3]                = {0.004, -0.08, 0.17};
  const gdouble g2[3]                = {0.001, 0.05, -0.12};

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  _pop_order_free_alpha (test, 1.0, 3.0);

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 3);
  nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  g_test_expect_message ("NUMCOSMO", G_LOG_LEVEL_WARNING, "*expansion over the population parameters does not converge*");
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  g_test_assert_expected_messages ();
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (gsfs), ==, 1);

  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 0.0);

  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "alpha", 2.7, NULL);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 0.0);

  /* Not built again, so not counted again */
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (gsfs), ==, 1);

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * A galaxy whose interpolant was expanded over the population parameters is built again for
 * one population when the expansion stops: with pop-order set to 0, and with the parameters
 * of the population fixed; and expanded again when pop-order is set back.
 */
void
test_nc_galaxy_shape_factor_spectral_pop_order_off (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble g1[3]                = {0.004, -0.08, 0.17};
  const gdouble g2[3]                = {0.001, 0.05, -0.12};
  guint idx;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  _pop_order_free_alpha (test, 1.45, 1.65);
  g_assert_true (ncm_model_param_index_from_name (NCM_MODEL (test->pop), "alpha", &idx, NULL));

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 7);
  nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 1.0e-6);

  /* pop-order 0: the evaluation builds the interpolant of the current population */
  nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 0);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 1.0e-6);

  /* Expanded again, then the parameters fixed */
  nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 7);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 1.0e-6);

  ncm_model_param_set_ftype (NCM_MODEL (test->pop), idx, NCM_PARAM_TYPE_FIXED);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  _pop_order_assert_integral (test, gsfs, data, 3, g1, g2, 1.0e-6);

  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (gsfs), ==, 0);

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/* The counts of nc_galaxy_shape_factor_spectral_get_fallback_stats() for reason */
static void
_fallback_assert_counts (NcGalaxyShapeFactorSpectralFallback reason, const guint64 lines, const guint64 shears)
{
  guint64 l, sh;
  gdouble sec;

  nc_galaxy_shape_factor_spectral_get_fallback_stats (reason, &l, &sh, &sec);
  g_assert_cmpuint (l, ==, lines);
  g_assert_cmpuint (sh, ==, shears);
  g_assert_cmpfloat (sec, >=, 0.0);
}

/*
 * pop-order with the population outside the box of the expansion. Along a line, the galaxy
 * is tabulated as without use-interp, the same values bit for bit, and the line is counted;
 * inside the box the interpolant is used, the values of eval_nodes, and nothing is counted.
 * A shear evaluated alone outside the box is integrated and counted.
 */
void
test_nc_galaxy_shape_factor_spectral_pop_order_fallback_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyPositionFactorData *p_data    = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data    = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorSpectral *gsfs     = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs_int = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf              = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactor *gsf_int          = NC_GALAXY_SHAPE_FACTOR (gsfs_int);
  NcGalaxyShapeFactorData *data         = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  NcGalaxyShapeFactorData *data_int     = nc_galaxy_shape_factor_data_new (gsf_int, test->mset, p_data, z_data);
  GPtrArray *data_array                 = g_ptr_array_new ();
  GPtrArray *data_int_array             = g_ptr_array_new ();
  const gdouble u_1                     = 1.05 * cos (0.7);
  const gdouble u_2                     = 1.05 * sin (0.7);
  const gdouble b_1                     = 0.004;
  const gdouble b_2                     = -0.003;
  const gdouble e_1                     = -0.25;
  const gdouble e_2                     = 0.35;
  const guint n                         = 40;
  gdouble t[TEST_LINE_N], out[TEST_LINE_N], ref[TEST_LINE_N];
  guint j;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  _pop_order_free_alpha (test, 1.45, 1.65);

  nc_galaxy_shape_factor_spectral_set_use_interp (gsfs, TRUE);
  nc_galaxy_shape_factor_spectral_set_pop_order (gsfs, 5);
  nc_galaxy_shape_factor_data_set (gsf, data, e_1, e_2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  nc_galaxy_shape_factor_data_set (gsf_int, data_int, e_1, e_2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_ptr_array_add (data_int_array, data_int);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf_int, test->mset, data_int_array, TRUE, TRUE));
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, data, 20.0);
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_pop_n_unconverged (gsfs), ==, 0);

  for (j = 0; j < n; j++)
    t[j] = 0.25 * j / (n - 1.0);

  /* Inside the box: the interpolant, nothing counted */
  nc_galaxy_shape_factor_spectral_reset_fallback_stats ();
  _test_eval_line_case (gsf, test->pop, data, n, 0.25, TRUE, 0.0);
  _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION, 0, 0);
  _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_RANGE, 0, 0);
  _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_CAP, 0, 0);

  /* Outside it: the table of the line, as without use-interp */
  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "alpha", 1.8, NULL);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf_int, test->mset, data_int_array, TRUE, TRUE));

  nc_galaxy_shape_factor_eval_line (gsf, test->pop, data, n, t, u_1, u_2, b_1, b_2, e_1, e_2, TRUE, out);
  nc_galaxy_shape_factor_eval_line (gsf_int, test->pop, data_int, n, t, u_1, u_2, b_1, b_2, e_1, e_2, TRUE, ref);

  for (j = 0; j < n; j++)
    g_assert_cmpfloat (out[j], ==, ref[j]);

  _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION, 1, 0);

  /* One shear alone: the integral, counted */
  {
    NcGalaxyShapePopData *pd = nc_galaxy_shape_pop_data_new (test->pop);
    const gdouble g_1        = b_1 + 0.1 * u_1;
    const gdouble g_2        = b_2 + 0.1 * u_2;

    nc_galaxy_shape_pop_prepare (test->pop, pd);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1, g_2, e_1, e_2), ==,
                       nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, pd, g_1, g_2, e_1, e_2, 0.1));
    _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION, 1, 1);
    nc_galaxy_shape_pop_data_unref (pd);
  }

  nc_galaxy_shape_factor_spectral_reset_fallback_stats ();
  _fallback_assert_counts (NC_GALAXY_SHAPE_FACTOR_SPECTRAL_FALLBACK_POPULATION, 0, 0);

  g_ptr_array_unref (data_array);
  g_ptr_array_unref (data_int_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_data_unref (data_int);
  nc_galaxy_shape_factor_spectral_free (gsfs);
  nc_galaxy_shape_factor_spectral_free (gsfs_int);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

