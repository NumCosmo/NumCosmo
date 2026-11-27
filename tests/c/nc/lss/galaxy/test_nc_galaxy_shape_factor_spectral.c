/***************************************************************************
 *            test_nc_galaxy_shape_factor_spectral.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_factor_spectral.c
 * Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
 * Minor contributions:
 * Caio Lima de Oliveira <caiolimadeoliveira@pm.me> (2026)
 * - /nc/galaxy_shape_factor_spectral/crit_hash.
 */

/*
 * Tests of NcGalaxyShapeFactorSpectral: ln P against the values certified by the Arb
 * reference tests/tools/nc_galaxy_shape_marginal_arb.c, the |g| = 1 closed form, the
 * |g| > 1 reduction, shears next to |g| = 1, the rotation of the frame, the hooks of NcGalaxyShapeFactor, the
 * accepted panels, the interpolant of ln P in g and over the free population parameters, the
 * properties and the errors for an input that cannot be integrated and for an interpolation
 * range that is not accepted.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include "test_nc_galaxy_shape_marginal_arb_cases.h"

#include <math.h>
#include <glib.h>
#include <glib-object.h>

/*
 * Largest |ln P - ln P_Arb| accepted at the default reltol (1e-5). The largest difference
 * measured is 6.9e-9, on the 20,062 certified points of the certification grid with
 * (|eps_obs| - 1) / std_noise <= 8.
 */
#define TEST_SPECTRAL_TOL (1.0e-7)

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
 * Largest |ln P (table over the population) - ln P (integral)| accepted at interp-abstol
 * 1e-3 for the table over alpha and beta (/interp/pop) and for the table over sigma
 * (/interp/pop/gauss). The largest differences measured are 5.5e-8 and 4.9e-8.
 */
#define TEST_SPECTRAL_POP2_TOL (1.0e-6)

/* The range of beta, and of alpha when it is free, in the tests of the table over the population */
#define TEST_SPECTRAL_POP_BETA_LB (1.6)
#define TEST_SPECTRAL_POP_BETA_UB (1.7)
#define TEST_SPECTRAL_POP_ALPHA_LB (1.5)
#define TEST_SPECTRAL_POP_ALPHA_UB (1.6)

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
  {"pop",         "*built for other settings:*the population parameter beta (*"},
  {"free",        "*built for other settings:*the population parameter alpha is free, not a variable of the interpolants*"},
  {"pop_bounds",  "*built for other settings:*the bounds of beta*"},
  {"pop_level",   "*built for other settings:*interp-pop-level-cap (4, built with 3)*"},
  {"pop_one",     "*built for other settings:*interp-pop-one-integral (0, built with 1)*"},
  {"pop_compress", "*built for other settings:*interp-compress-frac*"},
  {"pop_outside", "*the population parameter beta =*is outside the range*of the interpolant*"},
  {"abstol",     "*built for other settings:*interp-abstol*"},
  {"trim_frac",  "*built for other settings:*interp-trim-frac*"},
  {"conv",       "*built for other settings:*the ellipticity convention*"},
  {"eval",       "*the interpolant of the galaxy does not hold for this evaluation*"},
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
static NcGalaxyShapeFactorSpectral *_restrict_range_build (TestNcGalaxyShapeFactorSpectral *test, const gboolean restrict_range, const gdouble ra, NcGalaxyShapeFactorData **data);
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
void test_nc_galaxy_shape_factor_spectral_interp (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_array (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_nodes (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_nodes_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside (void);
void test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside_subprocess (void);
void test_nc_galaxy_shape_factor_spectral_restrict_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_cdelta (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_centre_box (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_restrict_range_cm_relation (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_eval_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_rebuild (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_data_prepare (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_stale (void);
void test_nc_galaxy_shape_factor_spectral_interp_stale_subprocess (void);
void test_nc_galaxy_shape_factor_spectral_interp_obs_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_obs (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_pop (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_pop_obs (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_pop_gauss (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_pop_one_integral (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch (void);
void test_nc_galaxy_shape_factor_spectral_interp_obs_pop_fixed (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_obs_mismatch_subprocess (gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_g_max (void);
void test_nc_galaxy_shape_factor_spectral_interp_pop_compress (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_compress_frac (void);
void test_nc_galaxy_shape_factor_spectral_interp_compress_frac_subprocess (void);
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
  g_test_add ("/nc/galaxy_shape_factor_spectral/crit_hash", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_crit_hash,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/hooks", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_hooks,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_line", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_line,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_nodes", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_nodes,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/eval_nodes_range", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_eval_nodes_range,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/eval_nodes_range/outside", &test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/eval_nodes_range/outside/subprocess", &test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside_subprocess);
  g_test_add ("/nc/galaxy_shape_factor_spectral/restrict_range", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/restrict_range/cdelta", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_cdelta,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/restrict_range/centre_box", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_centre_box,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/restrict_range/cm_relation", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_restrict_range_cm_relation,
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
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_pop,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop/obs", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_pop_obs,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop/gauss", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_pop_gauss,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop/one_integral", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_pop_one_integral,
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
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/obs/pop_fixed", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_obs_pop_fixed,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/stale", &test_nc_galaxy_shape_factor_spectral_interp_stale);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/stale/subprocess", &test_nc_galaxy_shape_factor_spectral_interp_stale_subprocess);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/obs/range", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_obs_range,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/pop/compress", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_pop_compress,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/compress_frac", &test_nc_galaxy_shape_factor_spectral_interp_compress_frac);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/compress_frac/subprocess", &test_nc_galaxy_shape_factor_spectral_interp_compress_frac_subprocess);
  g_test_add_func ("/nc/galaxy_shape_factor_spectral/interp/g_max", &test_nc_galaxy_shape_factor_spectral_interp_g_max);
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

/* Frees the parameter idx of pop over [lb, ub] */
static void
_test_pop_free (NcGalaxyShapePop *pop, const guint idx, const gdouble lb, const gdouble ub)
{
  ncm_model_param_set_lower_bound (NCM_MODEL (pop), idx, lb);
  ncm_model_param_set_upper_bound (NCM_MODEL (pop), idx, ub);
  ncm_model_param_set_ftype (NCM_MODEL (pop), idx, NCM_PARAM_TYPE_FREE);
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

  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_g_max (gsfs), ==, 0.3);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs), ==, 1.0e-5);
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (gsfs), ==, 2);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_trim_frac (gsfs), ==, 1.0e-4);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_compress_frac (gsfs), ==, 1.0e-2);

  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-4);
  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.2);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-4);
  nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, 4);
  nc_galaxy_shape_factor_spectral_set_interp_trim_frac (gsfs, 1.0e-3);
  nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, 0.25);
  dup = NC_GALAXY_SHAPE_FACTOR_SPECTRAL (ncm_serialize_dup_obj (ser, G_OBJECT (gsfs)));

  g_assert_cmpint (nc_galaxy_shape_factor_get_ellip_conv (NC_GALAXY_SHAPE_FACTOR (dup)), ==, NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_reltol (dup), ==, 1.0e-4);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_g_max (dup), ==, 0.2);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_abstol (dup), ==, 1.0e-4);
  g_assert_cmpuint (nc_galaxy_shape_factor_spectral_get_interp_min_radial_panels (dup), ==, 4);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_trim_frac (dup), ==, 1.0e-3);
  g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_interp_compress_frac (dup), ==, 0.25);

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
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, exp (ln_direct));

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
      nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

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
          g_assert_cmpfloat (nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, exp (ln_P));
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
 * A build follows the galaxy and the settings: after a change of the population parameters,
 * of eps_obs, of std_noise, of interp-abstol, of interp-g-max, of reltol, of
 * interp-min-radial-panels and of interp-trim-frac, a new build gives at one shear the
 * integral for the new values within the tolerance, and a value that differs from the one
 * before the change where the change moves ln P by more than the tolerance. A galaxy without
 * an interpolant takes the integral.
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
  NcGalaxyShapeFactorData *data, *bare;
  gdouble ln_P, ln_P_prev;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  bare = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  nc_galaxy_shape_factor_data_set (gsf, bare, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_ptr_array_add (data_array, bare);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, ln_P);

  /* Without an interpolant the value is the integral */
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, bare, g1, g2, e1, e2), ==,
                     nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, bare->pop_data, g1, g2, e1, e2, s));

  /* The population parameters */
  ln_P_prev = ln_P;
  ncm_model_param_set_by_name (NCM_MODEL (test->pop), "beta", 2.2, NULL);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The observed ellipticity */
  ln_P_prev = ln_P;
  e1        = 0.25;
  e2        = -0.3;
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The noise */
  ln_P_prev = ln_P;
  s         = 0.06;
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g1, g2, e1, e2, s), TEST_SPECTRAL_INTERP_TOL);
  g_assert_cmpfloat (fabs (ln_P - ln_P_prev), >, 1.0e-3);

  /* The tolerances and the range: a new interpolant, the same function */
  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 1.0e-4);

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-5);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2), ==, ln_P_prev);

  nc_galaxy_shape_factor_spectral_set_interp_g_max (gsfs, 0.2);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_reltol (gsfs, 1.0e-7);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_interp_min_radial_panels (gsfs, 4);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 2.0 * TEST_SPECTRAL_INTERP_TOL);

  /* The trim: the stored series differs from the full one by at most trim-frac interp-abstol */
  ln_P_prev = ln_P;
  nc_galaxy_shape_factor_spectral_set_interp_trim_frac (gsfs, 0.1);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g1, g2, e1, e2);
  g_assert_cmpfloat (ln_P, !=, ln_P_prev);
  g_assert_cmpfloat_with_epsilon (ln_P, ln_P_prev, 0.1 * nc_galaxy_shape_factor_spectral_get_interp_abstol (gsfs));

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_data_unref (bare);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/*
 * nc_galaxy_shape_factor_data_prepare() builds nothing: a galaxy without an interpolant
 * takes the integral after it. nc_galaxy_shape_factor_spectral_build_interp() builds one,
 * which data_prepare() then accepts, and its value agrees with the integral within
 * TEST_SPECTRAL_INTERP_TOL.
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
  NcGalaxyShapeFactorData *built, *direct;
  guint i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < 2; i++)
  {
    NcGalaxyShapeFactorData *data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);

    nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
  }

  built  = g_ptr_array_index (data_array, 0);
  direct = g_ptr_array_index (data_array, 1);

  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, direct, 3.0);
  g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, direct, g1, g2, e1, e2), ==,
                     nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, direct->pop_data, g1, g2, e1, e2, s));

  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, built, 3.0);
  nc_galaxy_shape_factor_data_prepare (gsf, test->mset, built, 3.0);
  g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, built, g1, g2, e1, e2),
                                  nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, built->pop_data, g1, g2, e1, e2, s),
                                  TEST_SPECTRAL_INTERP_TOL);

  g_ptr_array_unref (data_array);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/* An interpolant evaluated after a change of a setting, without a new build, aborts */
void
test_nc_galaxy_shape_factor_spectral_interp_stale (void)
{
  g_test_trap_subprocess ("/nc/galaxy_shape_factor_spectral/interp/stale/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the interpolant of the galaxy does not hold for this evaluation*");
}

void
test_nc_galaxy_shape_factor_spectral_interp_stale_subprocess (void)
{
  TestNcGalaxyShapeFactorSpectral test;
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
  GPtrArray *data_array             = g_ptr_array_new ();
  NcGalaxyPositionFactorData *p_data;
  NcGalaxyRedshiftFactorData *z_data;
  NcGalaxyShapeFactorData *data;

  test_nc_galaxy_shape_factor_spectral_new (&test, NULL);
  p_data      = nc_galaxy_position_factor_data_new (test.position_factor, test.mset);
  z_data      = nc_galaxy_redshift_factor_data_new (test.redshift_factor, test.mset);
  p_data->ra  = 0.3;
  p_data->dec = 0.0;
  data        = nc_galaxy_shape_factor_data_new (gsf, test.mset, p_data, z_data);

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test.mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test.mset, data, 3.0);

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);
  nc_galaxy_shape_factor_eval_ln_marginal (gsf, test.pop, data, 0.12, -0.07, 0.3, 0.2);
}

/*
 * A catalog prepared with restrict-range, read by a shape factor without it: each loaded
 * interpolant keeps the range it was built with, so a shear above that range takes the
 * integral, bit for bit, also through eval_nodes, and a shear inside it agrees with the
 * integral within interp-abstol.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_obs_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  const gdouble g_far[2]   = {0.29 * cos (0.4), 0.29 * sin (0.4)};
  const gdouble g_near[2]  = {0.001, 0.0};
  NcGalaxyShapeFactorSpectral *gsfs_r, *gsfs_f;
  NcGalaxyShapeFactorData *data_r, *data_f;
  GPtrArray *data_array;
  NcGalaxyWLObs *obs;
  guint log10M_idx = 0;
  gdouble ln_direct, ln_f, ln_n;

  g_assert_true (ncm_model_param_index_from_name (hms, "log10MDelta", &log10M_idx, NULL));
  ncm_model_param_set_upper_bound (hms, log10M_idx, 14.5);

  gsfs_r = _restrict_range_build (test, TRUE, 0.8, &data_r);
  gsfs_f = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  g_assert_true (nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs_f));
  nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs_f, FALSE);

  {
    GList *columns = nc_galaxy_shape_factor_data_required_columns (data_r);
    GStrv names    = g_new0 (gchar *, g_list_length (columns) + 1);
    GList *l;
    guint i = 0;

    for (l = columns; l != NULL; l = l->next)
      names[i++] = g_strdup (l->data);

    obs = nc_galaxy_wl_obs_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, NC_WL_ELLIPTICITY_FRAME_CELESTIAL, 1, names);
    nc_galaxy_shape_factor_data_write_row (data_r, obs, 0);
    nc_galaxy_shape_factor_spectral_set_obs_source (gsfs_r, test->pop, obs);
    g_strfreev (names);
    g_list_free_full (columns, g_free);
  }

  {
    NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
    NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);

    data_f     = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->mset, p_data, z_data);
    data_array = g_ptr_array_new ();
    nc_galaxy_shape_factor_data_read_row (data_f, obs, 0);
    g_ptr_array_add (data_array, data_f);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->mset, data_array, TRUE, TRUE));
    nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->mset, data_f, 3.0);
    nc_galaxy_position_factor_data_unref (p_data);
    nc_galaxy_redshift_factor_data_unref (z_data);
  }

  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_f, test->pop, data_f->pop_data, g_far[0], g_far[1], 0.3, 0.2, 0.1);
  ln_f      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->pop, data_f, g_far[0], g_far[1], 0.3, 0.2);
  g_assert_cmpfloat (ln_f, ==, ln_direct);

  NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsfs_f)->eval_nodes (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->pop, data_f, 1, &g_far[0], &g_far[1], 0.3, 0.2, TRUE, &ln_n);
  g_assert_cmpfloat (ln_n, ==, ln_direct);

  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_f, test->pop, data_f->pop_data, g_near[0], g_near[1], 0.3, 0.2, 0.1);
  ln_f      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->pop, data_f, g_near[0], g_near[1], 0.3, 0.2);
  g_assert_cmpfloat (fabs (ln_f - ln_direct), <=, 1.0e-5);

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data_r);
  nc_galaxy_shape_factor_data_unref (data_f);
  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_shape_factor_spectral_free (gsfs_r);
  nc_galaxy_shape_factor_spectral_free (gsfs_f);
}

/* An interpolation range of 1/3 or more aborts */

/* The numbers a table stores, and how many of its panels are TUCKER */
static gsize
_test_table_stored (NcmSpectralTensor *st, guint *n_tucker)
{
  gsize stored = 0;
  guint p;

  *n_tucker = 0;

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
  {
    if (ncm_spectral_tensor_get_panel_format (st, p) == NCM_SPECTRAL_TENSOR_PANEL_FORMAT_TUCKER)
    {
      GArray *ranks;
      GArray *d = ncm_spectral_tensor_get_panel_tucker (st, p, &ranks);

      stored += d->len;
      (*n_tucker)++;
      g_array_unref (ranks);
      g_array_unref (d);
    }
    else
    {
      GArray *lower, *upper, *K;
      GArray *c = ncm_spectral_tensor_get_panel (st, p, &lower, &upper, &K);

      stored += c->len;
      g_array_unref (lower);
      g_array_unref (upper);
      g_array_unref (K);
      g_array_unref (c);
    }
  }

  return stored;
}

/*
 * The compression of a table over alpha and beta, at interp-abstol 1e-3: built with
 * interp-compress-frac 0 the table has DENSE panels only; built with the default fraction it
 * has TUCKER panels and stores fewer numbers, and ln P of
 * nc_galaxy_shape_factor_eval_ln_marginal() and of eval_nodes differs from the DENSE values
 * by at most the fraction times interp-abstol, at the ends and inside the ranges of both
 * parameters. nc_galaxy_shape_factor_spectral_data_take_interp() gives each table and
 * leaves the galaxy without one.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_pop_compress (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const gdouble shears[][2] = {
    {
      0.0, 0.0
    }, {
      0.12, -0.07
    }, {
      -0.2, 0.15
    }, {
      0.05, 0.28
    }
  };
  const gdouble betas[]              = {TEST_SPECTRAL_POP_BETA_LB, 1.63, TEST_SPECTRAL_POP_BETA_UB};
  const gdouble alphas[]             = {TEST_SPECTRAL_POP_ALPHA_LB, 1.55, TEST_SPECTRAL_POP_ALPHA_UB};
  const gdouble abstol               = 1.0e-3;
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  const gdouble frac                 = nc_galaxy_shape_factor_spectral_get_interp_compress_frac (gsfs);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorClass *klass    = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  gdouble g_1[G_N_ELEMENTS (shears)], g_2[G_N_ELEMENTS (shears)], out[G_N_ELEMENTS (shears)];
  gdouble ln_dense[G_N_ELEMENTS (alphas)][G_N_ELEMENTS (betas)][G_N_ELEMENTS (shears)];
  NcmSpectralTensor *st_dense, *st_comp;
  guint n_tucker_dense, n_tucker_comp;
  gsize stored_dense, stored_comp;
  guint pass, a, b, l;

  g_assert_cmpfloat (frac, >, 0.0);

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (l = 0; l < G_N_ELEMENTS (shears); l++)
  {
    g_1[l] = shears[l][0];
    g_2[l] = shears[l][1];
  }

  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, 3);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, abstol);
  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);
  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_ALPHA, TEST_SPECTRAL_POP_ALPHA_LB, TEST_SPECTRAL_POP_ALPHA_UB);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, 0.2, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);

  st_dense = st_comp = NULL;

  /* Pass 0: no compression; pass 1: the default fraction */
  for (pass = 0; pass < 2; pass++)
  {
    nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, (pass == 0) ? 0.0 : frac);

    for (a = 0; a < G_N_ELEMENTS (alphas); a++)
    {
      for (b = 0; b < G_N_ELEMENTS (betas); b++)
      {
        ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_ALPHA, alphas[a]);
        ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, betas[b]);
        g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

        if ((a == 0) && (b == 0))
          nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

        klass->eval_nodes (gsf, test->pop, data, G_N_ELEMENTS (shears), g_1, g_2, e1, e2, TRUE, out);

        for (l = 0; l < G_N_ELEMENTS (shears); l++)
        {
          const gdouble ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1[l], g_2[l], e1, e2);

          g_assert_cmpfloat_with_epsilon (out[l], ln_P, 1.0e-12);

          if (pass == 0)
          {
            ln_dense[a][b][l] = ln_P;
          }
          else
          {
            g_assert_cmpfloat (fabs (ln_P - ln_dense[a][b][l]), <=, frac * abstol);
            g_assert_cmpfloat (fabs (out[l] - ln_dense[a][b][l]), <=, frac * abstol);
          }
        }
      }
    }

    if (pass == 0)
      st_dense = nc_galaxy_shape_factor_spectral_data_take_interp (data);
    else
      st_comp = nc_galaxy_shape_factor_spectral_data_take_interp (data);

    g_assert_null (nc_galaxy_shape_factor_spectral_data_take_interp (data));
  }

  g_assert_true (NCM_IS_SPECTRAL_TENSOR (st_dense));
  g_assert_true (NCM_IS_SPECTRAL_TENSOR (st_comp));

  stored_dense = _test_table_stored (st_dense, &n_tucker_dense);
  stored_comp  = _test_table_stored (st_comp, &n_tucker_comp);

  g_assert_cmpuint (n_tucker_dense, ==, 0);
  g_assert_cmpuint (n_tucker_comp, >, 0);
  g_assert_cmpuint (stored_comp, <, stored_dense);
  g_test_message ("stored numbers: %zu DENSE, %zu compressed, %u TUCKER panels of %u", stored_dense, stored_comp, n_tucker_comp,
                  ncm_spectral_tensor_get_n_panels (st_comp));

  ncm_spectral_tensor_free (st_dense);
  ncm_spectral_tensor_free (st_comp);
  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

void
test_nc_galaxy_shape_factor_spectral_interp_compress_frac (void)
{
  g_test_trap_subprocess ("/nc/galaxy_shape_factor_spectral/interp/compress_frac/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*nc_galaxy_shape_factor_spectral_set_interp_compress_frac: 1.5 is outside*");
}

void
test_nc_galaxy_shape_factor_spectral_interp_compress_frac_subprocess (void)
{
  NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);

  nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, 1.5);
}

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
 * A catalog of the galaxies of test_spectral_obs_gal prepared by gsfs for
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

    nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
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
 * carries no interpolants, and the galaxies take the integral.
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

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  obs = _test_obs_prepared_new (test, gsfs_w, g1, g2, ln_P);

  for (i = 0; i < G_N_ELEMENTS (test_spectral_obs_gal); i++)
    g_assert_true (NCM_IS_SPECTRAL_TENSOR (nc_galaxy_wl_obs_peek_shape_interp (obs, i)));

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

  /* A change of a row: no interpolant is loaded, and no galaxy builds one: each takes the integral */
  nc_galaxy_wl_obs_set (obs, NC_GALAXY_SHAPE_FACTOR_COL_M, 2, 0.0, NULL);
  g_assert_null (nc_galaxy_wl_obs_peek_shape_interp (obs, 0));
  data_array = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2), ==,
                       nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data->pop_data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2,
                                                                       data->std_noise));
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

  const gboolean pop_table = g_str_has_prefix (mismatch, "pop_");

  test_nc_galaxy_shape_factor_spectral_new (&test, NULL);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_w, 3);

  if (pop_table)
    _test_pop_free (test.pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);

  obs = _test_obs_prepared_new (&test, gsfs_w, 0.12, -0.07, ln_P);

  gsfs_r = nc_galaxy_shape_factor_spectral_new (g_str_equal (mismatch, "conv") ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_r, 3);

  if (g_str_equal (mismatch, "pop"))
    ncm_model_param_set_by_name (NCM_MODEL (test.pop), "beta", 1.65 * (1.0 + 1.0e-12), NULL);
  else if (g_str_equal (mismatch, "free"))
    ncm_model_param_set_ftype (NCM_MODEL (test.pop), 0, NCM_PARAM_TYPE_FREE);
  else if (g_str_equal (mismatch, "abstol"))
    nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, 1.0e-4);
  else if (g_str_equal (mismatch, "trim_frac"))
    nc_galaxy_shape_factor_spectral_set_interp_trim_frac (gsfs_r, 0.0);
  else if (g_str_equal (mismatch, "pop_bounds"))
    ncm_model_param_set_upper_bound (NCM_MODEL (test.pop), NC_GALAXY_SHAPE_POP_BETA_BETA, 1.8);
  else if (g_str_equal (mismatch, "pop_level"))
    nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_r, 4);
  else if (g_str_equal (mismatch, "pop_one"))
    nc_galaxy_shape_factor_spectral_set_interp_pop_one_integral (gsfs_r, FALSE);
  else if (g_str_equal (mismatch, "pop_compress"))
    nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs_r, 0.5);

  data_array = _test_obs_read (&test, NC_GALAXY_SHAPE_FACTOR (gsfs_r), obs);
  data       = g_ptr_array_index (data_array, 0);
  nc_galaxy_shape_factor_data_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test.mset, data, 3.0);

  /* Passed the check: a change of the population afterwards is found by the evaluation (beta outside the table for pop_outside) */
  ncm_model_param_set_by_name (NCM_MODEL (test.pop), "beta", 2.0, NULL);
  nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test.pop, data, 0.12, -0.07, data->epsilon_obs_1, data->epsilon_obs_2);
}

/*
 * A catalog prepared with beta a variable on [TEST_SPECTRAL_POP_BETA_LB,
 * TEST_SPECTRAL_POP_BETA_UB], read by a fit that holds beta fixed at 1.63: the check accepts
 * it, and ln P of each galaxy equals, bit for bit, the value of a fit with beta free at 1.63
 * (both evaluate the same table at the same point) and agrees with the integral within
 * TEST_SPECTRAL_OBS_INTERP_ABSTOL.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_obs_pop_fixed (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs_w = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs_r = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf_r          = NC_GALAXY_SHAPE_FACTOR (gsfs_r);
  const gdouble g1                    = 0.12;
  const gdouble g2                    = -0.07;
  gdouble ln_P[G_N_ELEMENTS (test_spectral_obs_gal)];
  gdouble ln_fixed[G_N_ELEMENTS (test_spectral_obs_gal)];
  GPtrArray *data_fixed, *data_free;
  NcGalaxyWLObs *obs;
  guint i;

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_w, 3);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_r, 3);

  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);
  obs = _test_obs_prepared_new (test, gsfs_w, g1, g2, ln_P);

  ncm_model_param_set_ftype (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, 1.63);
  data_fixed = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_fixed->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_fixed, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    ln_fixed[i] = nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2);
    g_assert_cmpfloat_with_epsilon (ln_fixed[i],
                                    nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data->pop_data, g1, g2,
                                                                                    data->epsilon_obs_1, data->epsilon_obs_2, data->std_noise),
                                    TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  }

  ncm_model_param_set_ftype (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, NCM_PARAM_TYPE_FREE);
  data_free = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_free->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_free, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2), ==, ln_fixed[i]);
  }

  g_ptr_array_unref (data_fixed);
  g_ptr_array_unref (data_free);
  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_shape_factor_spectral_free (gsfs_w);
  nc_galaxy_shape_factor_spectral_free (gsfs_r);
}

/*
 * The table over the population: with beta free on [1.6, 1.7] the interpolant of a galaxy
 * agrees with the integral within TEST_SPECTRAL_INTERP_TOL at the default interp-abstol, at
 * both ends of the range and inside it, by eval_ln_marginal and by eval_nodes (which agree
 * to rounding). A change of the fixed alpha builds it again for the new value. With alpha
 * free on [1.5, 1.6] too, at interp-abstol 1e-3, it agrees within TEST_SPECTRAL_POP2_TOL at
 * the nine pairs of the ends and middles.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_pop (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const gdouble shears[][2] = {
    {
      0.0, 0.0
    }, {
      0.12, -0.07
    }, {
      -0.2, 0.15
    }, {
      0.05, 0.28
    }
  };
  const gdouble betas[]              = {TEST_SPECTRAL_POP_BETA_LB, 1.63, TEST_SPECTRAL_POP_BETA_UB};
  const gdouble alphas[]             = {TEST_SPECTRAL_POP_ALPHA_LB, 1.55, TEST_SPECTRAL_POP_ALPHA_UB};
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyShapeFactorClass *klass    = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorData *data      = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  GPtrArray *data_array              = g_ptr_array_new ();
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  gdouble g_1[G_N_ELEMENTS (shears)], g_2[G_N_ELEMENTS (shears)], out[G_N_ELEMENTS (shears)];
  guint a, b, l;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (l = 0; l < G_N_ELEMENTS (shears); l++)
  {
    g_1[l] = shears[l][0];
    g_2[l] = shears[l][1];
  }

  /* The build against the integral: no compression, tested in /interp/pop/compress */
  nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, 0.0);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, 3);
  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);

  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, data);

  for (b = 0; b < G_N_ELEMENTS (betas); b++)
  {
    ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, betas[b]);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    if (b == 0)
      nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

    klass->eval_nodes (gsf, test->pop, data, G_N_ELEMENTS (shears), g_1, g_2, e1, e2, TRUE, out);

    for (l = 0; l < G_N_ELEMENTS (shears); l++)
    {
      const gdouble ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1[l], g_2[l], e1, e2);

      g_assert_cmpfloat_with_epsilon (ln_P, nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g_1[l], g_2[l], e1, e2, 0.1),
                                      TEST_SPECTRAL_INTERP_TOL);
      g_assert_cmpfloat_with_epsilon (out[l], ln_P, 1.0e-12);
    }
  }

  /* A fixed parameter: a new table for its new value */
  ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_ALPHA, 1.52);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
  g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, 0.12, -0.07, e1, e2),
                                  nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, 0.12, -0.07, e1, e2, 0.1),
                                  TEST_SPECTRAL_INTERP_TOL);

  /* alpha and beta free */
  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_ALPHA, TEST_SPECTRAL_POP_ALPHA_LB, TEST_SPECTRAL_POP_ALPHA_UB);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);
  nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, 0.2, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);

  for (a = 0; a < G_N_ELEMENTS (alphas); a++)
  {
    for (b = 0; b < G_N_ELEMENTS (betas); b++)
    {
      ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_ALPHA, alphas[a]);
      ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, betas[b]);
      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

      if ((a == 0) && (b == 0))
        nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

      for (l = 0; l < G_N_ELEMENTS (shears); l++)
        g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1[l], g_2[l], e1, e2),
                                        nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g_1[l], g_2[l], e1, e2, 0.2),
                                        TEST_SPECTRAL_POP2_TOL);
    }
  }

  g_ptr_array_unref (data_array);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_factor_spectral_free (gsfs);
}

/*
 * A catalog prepared with beta free on [1.6, 1.7] stores tables with three variables and
 * gives, after a serialization round trip, the values it was prepared with. A loaded table is
 * never rebuilt, so the evaluations at other values of beta inside the range, which agree
 * with the integral within TEST_SPECTRAL_OBS_INTERP_ABSTOL, show that it holds for them.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_pop_obs (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcGalaxyShapeFactorSpectral *gsfs_w = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactorSpectral *gsfs_r = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf_r          = NC_GALAXY_SHAPE_FACTOR (gsfs_r);
  const gdouble betas[]               = {TEST_SPECTRAL_POP_BETA_LB, 1.68, TEST_SPECTRAL_POP_BETA_UB};
  const gdouble g1                    = 0.12;
  const gdouble g2                    = -0.07;
  gdouble ln_P[G_N_ELEMENTS (test_spectral_obs_gal)];
  NcGalaxyWLObs *obs;
  GPtrArray *data_array;
  guint i, b;

  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_w, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs_r, TEST_SPECTRAL_OBS_INTERP_ABSTOL);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_w, 3);
  nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs_r, 3);
  _test_pop_free (test->pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);

  obs = _test_obs_prepared_new (test, gsfs_w, g1, g2, ln_P);

  for (i = 0; i < G_N_ELEMENTS (test_spectral_obs_gal); i++)
    g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (nc_galaxy_wl_obs_peek_shape_interp (obs, i)), ==, 3);

  data_array = _test_obs_read (test, gsf_r, obs);

  for (i = 0; i < data_array->len; i++)
  {
    NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

    nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
    g_assert_cmpfloat (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2), ==, ln_P[i]);
  }

  for (b = 0; b < G_N_ELEMENTS (betas); b++)
  {
    ncm_model_param_set (NCM_MODEL (test->pop), NC_GALAXY_SHAPE_POP_BETA_BETA, betas[b]);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf_r, test->mset, data_array, TRUE, TRUE));

    for (i = 0; i < data_array->len; i++)
    {
      NcGalaxyShapeFactorData *data = g_ptr_array_index (data_array, i);

      nc_galaxy_shape_factor_data_prepare (gsf_r, test->mset, data, 3.0);
      g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf_r, test->pop, data, g1, g2, data->epsilon_obs_1, data->epsilon_obs_2),
                                      nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data->pop_data, g1, g2,
                                                                                      data->epsilon_obs_1, data->epsilon_obs_2, data->std_noise),
                                      TEST_SPECTRAL_OBS_INTERP_ABSTOL);
    }
  }

  g_ptr_array_unref (data_array);
  nc_galaxy_wl_obs_free (obs);
  nc_galaxy_shape_factor_spectral_free (gsfs_w);
  nc_galaxy_shape_factor_spectral_free (gsfs_r);
}

/*
 * The table over sigma of a truncated Gaussian population, which replaces the Beta population
 * of the fixture: for TRACE_DET with sigma free on its default bounds, [0.2, 1], and for TRACE
 * with sigma free on [0.3, 0.6], at interp-abstol 1e-3, the interpolant agrees with the
 * integral within TEST_SPECTRAL_POP2_TOL at both ends of the range and inside it.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_pop_gauss (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const gdouble shears[][2] = {
    {
      0.0, 0.0
    }, {
      0.12, -0.07
    }, {
      -0.2, 0.15
    }, {
      0.05, 0.28
    }
  };
  NcGalaxyShapePop *pop              = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  const gdouble s                    = 0.2;
  guint j, k, l;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  g_assert_cmpfloat (ncm_model_param_get_lower_bound (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA), ==, 0.2);
  g_assert_cmpfloat (ncm_model_param_get_upper_bound (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA), ==, 1.0);
  ncm_model_param_set_ftype (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, NCM_PARAM_TYPE_FREE);
  ncm_mset_set (test->mset, NCM_MODEL (pop), NULL);

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();
    const gdouble sigmas[]            = {(j == 0) ? 0.2 : 0.3, 0.45, (j == 0) ? 1.0 : 0.6};

    if (j == 1)
      _test_pop_free (pop, NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.3, 0.6);

    /* The build against the integral: no compression, tested in /interp/pop/compress */
    nc_galaxy_shape_factor_spectral_set_interp_compress_frac (gsfs, 0.0);
    nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);
    nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, 3);
    nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);

    for (k = 0; k < G_N_ELEMENTS (sigmas); k++)
    {
      ncm_model_param_set (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, sigmas[k]);
      g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

      if (k == 0)
        nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

      for (l = 0; l < G_N_ELEMENTS (shears); l++)
        g_assert_cmpfloat_with_epsilon (nc_galaxy_shape_factor_eval_ln_marginal (gsf, pop, data, shears[l][0], shears[l][1], e1, e2),
                                        nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, pop, data->pop_data, shears[l][0], shears[l][1], e1, e2, s),
                                        TEST_SPECTRAL_POP2_TOL);
    }

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_pop_free (pop);
}

/*
 * interp-pop-one-integral, TRUE by default: the interpolants built with one integral per
 * shear for the populations of a call agree with those built with one integral per
 * population, at the corners and the center of the range of the free parameters and at four
 * shears. The cases: beta free (TRACE, default interp-abstol), alpha and beta free
 * (TRACE_DET, interp-abstol 1e-3) and sigma free of the truncated Gaussian (TRACE,
 * interp-abstol 1e-3). The largest difference measured is 1.8e-10.
 */
void
test_nc_galaxy_shape_factor_spectral_interp_pop_one_integral (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const gdouble shears[][2] = {
    {
      0.0, 0.0
    }, {
      0.12, -0.07
    }, {
      -0.2, 0.15
    }, {
      0.05, 0.28
    }
  };
  const gdouble corner[]             = {0.0, 1.0, 0.5};
  NcGalaxyShapePop *gauss            = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  const gdouble e1                   = 0.3;
  const gdouble e2                   = 0.2;
  guint c, one, a, b, l;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (c = 0; c < 3; c++)
  {
    const NcGalaxyWLObsEllipConv conv = (c == 1) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE;
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (conv);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapePop *pop             = (c == 2) ? gauss : test->pop;
    const gdouble s                   = (c == 0) ? 0.1 : 0.2;
    const guint n_a                   = (c == 1) ? G_N_ELEMENTS (corner) : 1;
    GPtrArray *data_array             = g_ptr_array_new ();
    NcGalaxyShapeFactorData *data;
    gdouble vals[2][G_N_ELEMENTS (corner) * G_N_ELEMENTS (corner) * G_N_ELEMENTS (shears)];

    g_assert_true (nc_galaxy_shape_factor_spectral_get_interp_pop_one_integral (gsfs));

    if (c == 0)
    {
      _test_pop_free (pop, NC_GALAXY_SHAPE_POP_BETA_BETA, TEST_SPECTRAL_POP_BETA_LB, TEST_SPECTRAL_POP_BETA_UB);
    }
    else if (c == 1)
    {
      _test_pop_free (pop, NC_GALAXY_SHAPE_POP_BETA_ALPHA, TEST_SPECTRAL_POP_ALPHA_LB, TEST_SPECTRAL_POP_ALPHA_UB);
    }
    else
    {
      _test_pop_free (pop, NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.3, 0.6);
      ncm_mset_set (test->mset, NCM_MODEL (pop), NULL);
    }

    data = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    nc_galaxy_shape_factor_spectral_set_interp_pop_level_cap (gsfs, 3);

    if (c > 0)
      nc_galaxy_shape_factor_spectral_set_interp_abstol (gsfs, 1.0e-3);

    nc_galaxy_shape_factor_data_set (gsf, data, e1, e2, s, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);

    for (one = 0; one < 2; one++)
    {
      guint q = 0;

      nc_galaxy_shape_factor_spectral_set_interp_pop_one_integral (gsfs, one == 0);

      for (a = 0; a < n_a; a++)
      {
        for (b = 0; b < G_N_ELEMENTS (corner); b++)
        {
          if (c == 2)
          {
            ncm_model_param_set (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_GAUSS_SIGMA, 0.3 + 0.3 * corner[b]);
          }
          else
          {
            if (c == 1)
              ncm_model_param_set (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_BETA_ALPHA,
                                   TEST_SPECTRAL_POP_ALPHA_LB + (TEST_SPECTRAL_POP_ALPHA_UB - TEST_SPECTRAL_POP_ALPHA_LB) * corner[a]);

            ncm_model_param_set (NCM_MODEL (pop), NC_GALAXY_SHAPE_POP_BETA_BETA,
                                 TEST_SPECTRAL_POP_BETA_LB + (TEST_SPECTRAL_POP_BETA_UB - TEST_SPECTRAL_POP_BETA_LB) * corner[b]);
          }

          g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

          if (q == 0)
            nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

          for (l = 0; l < G_N_ELEMENTS (shears); l++)
            vals[one][q++] = nc_galaxy_shape_factor_eval_ln_marginal (gsf, pop, data, shears[l][0], shears[l][1], e1, e2);
        }
      }

      if (one == 1)
      {
        guint i;

        for (i = 0; i < q; i++)
          g_assert_cmpfloat_with_epsilon (vals[0][i], vals[1][i], 1.0e-8);
      }
    }

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_spectral_free (gsfs);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_shape_pop_free (gauss);
}

/*
 * The eval_line of the class at 50 shears b + t u, for two observed ellipticities, with
 * line-abstol 1e-6, for a galaxy without an interpolant. For t from 0.002 to 0.102 it
 * returns TRUE, and ln P and P agree with the integral at each shear within line-abstol
 * (largest difference measured: 1.2e-10 and 2.9e-11). For t up to 0.292 and
 * e_obs = (-0.25, 0.35) the series does not reach the tolerance within the grids allowed, and it returns
 * FALSE with out unchanged. It also returns FALSE with line-abstol 0, for a galaxy with an
 * interpolant, for 5 shears, and for a line that reaches |g| > 0.9.
 */
void
test_nc_galaxy_shape_factor_spectral_eval_line (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const guint n             = 50;
  const gdouble e_obs[2][2] = {
    {
      -0.25, 0.35
    }, {
      0.6, -0.3
    }
  };
  const gdouble u_1                  = cos (0.7);
  const gdouble u_2                  = sin (0.7);
  const gdouble b_1                  = 0.003;
  const gdouble b_2                  = -0.002;
  const gdouble line_abstol          = 1.0e-6;
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  gdouble t[50], t_long[50], t_far[50], ln_out[50], lin_out[50];
  guint c, i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < n; i++)
  {
    t[i]      = 0.002 + 0.1 * i / (n - 1.0);
    t_long[i] = 0.002 + 0.29 * i / (n - 1.0);
    t_far[i]  = 0.002 + 0.95 * i / (n - 1.0);
  }

  for (c = 0; c < 2; c++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorClass *klass   = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();

    nc_galaxy_shape_factor_data_set (gsf, data, e_obs[c][0], e_obs[c][1], 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    /* Off by default */
    g_assert_cmpfloat (nc_galaxy_shape_factor_spectral_get_line_abstol (gsfs), ==, 0.0);
    ln_out[0] = 123.0;
    g_assert_false (klass->eval_line (gsf, test->pop, data, n, t, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));
    g_assert_cmpfloat (ln_out[0], ==, 123.0);

    nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, line_abstol);
    g_assert_true (klass->eval_line (gsf, test->pop, data, n, t, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));
    g_assert_true (klass->eval_line (gsf, test->pop, data, n, t, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], FALSE, lin_out));

    for (i = 0; i < n; i++)
    {
      const gdouble g_1  = b_1 + t[i] * u_1;
      const gdouble g_2  = b_2 + t[i] * u_2;
      const gdouble ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1, g_2, e_obs[c][0], e_obs[c][1]);

      g_assert_cmpfloat (fabs (ln_out[i] - ln_P), <=, line_abstol);
      g_assert_cmpfloat (fabs (lin_out[i] / exp (ln_P) - 1.0), <=, 2.0 * line_abstol);
    }

    /* Declined: not converged, 5 shears, a line beyond |g| = 0.9, a galaxy with an interpolant */
    ln_out[0] = 123.0;

    if (c == 0)
      g_assert_false (klass->eval_line (gsf, test->pop, data, n, t_long, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));

    g_assert_false (klass->eval_line (gsf, test->pop, data, 5, t, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));
    g_assert_false (klass->eval_line (gsf, test->pop, data, n, t_far, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));
    nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
    g_assert_false (klass->eval_line (gsf, test->pop, data, n, t, u_1, u_2, b_1, b_2, e_obs[c][0], e_obs[c][1], TRUE, ln_out));
    g_assert_cmpfloat (ln_out[0], ==, 123.0);

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/*
 * The eval_nodes of the class gives, for one galaxy, what eval_ln_marginal and eval_marginal
 * give shear by shear, at 40 shears on a line and around it, inside and outside
 * interp-g-max (0.3). ln P by integration is the same value. ln P from the interpolant agrees
 * to rounding (the array evaluation of NcmSpectralTensor sums in another order): measured
 * 4.4e-16, bound 1e-14. P from eval_nodes is the exponential of ln P in a loop that calls the
 * vector exp of libmvec, within 4 ulp of exp (): measured 4.4e-16 of P, bound 1e-14 of P.
 */
void
test_nc_galaxy_shape_factor_spectral_eval_nodes (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const guint n                      = 40;
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  gdouble g_1[40], g_2[40], ln_out[40], lin_out[40];
  guint j, i;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;

  for (i = 0; i < n; i++)
  {
    const gdouble t = 0.4 * i / (n - 1.0);

    g_1[i] = 0.002 + t * 0.8 + ((i % 5 == 0) ? 0.01 : 0.0);
    g_2[i] = -t * 0.6;
  }

  for (j = 0; j < 2; j++)
  {
    NcGalaxyShapeFactorSpectral *gsfs = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
    NcGalaxyShapeFactor *gsf          = NC_GALAXY_SHAPE_FACTOR (gsfs);
    NcGalaxyShapeFactorClass *klass   = NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsf);
    NcGalaxyShapeFactorData *data     = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
    GPtrArray *data_array             = g_ptr_array_new ();

    nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    if (j == 1)
      nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);

    klass->eval_nodes (gsf, test->pop, data, n, g_1, g_2, data->epsilon_obs_1, data->epsilon_obs_2, TRUE, ln_out);
    klass->eval_nodes (gsf, test->pop, data, n, g_1, g_2, data->epsilon_obs_1, data->epsilon_obs_2, FALSE, lin_out);

    for (i = 0; i < n; i++)
    {
      const gdouble ln_P = nc_galaxy_shape_factor_eval_ln_marginal (gsf, test->pop, data, g_1[i], g_2[i], data->epsilon_obs_1, data->epsilon_obs_2);
      const gdouble P    = nc_galaxy_shape_factor_eval_marginal (gsf, test->pop, data, g_1[i], g_2[i], data->epsilon_obs_1, data->epsilon_obs_2);

      g_assert_true (gsl_finite (ln_P));

      if ((j == 0) || (g_1[i] * g_1[i] + g_2[i] * g_2[i] > 0.09))
      {
        g_assert_cmpfloat (ln_out[i], ==, ln_P);
        g_assert_cmpfloat (fabs (lin_out[i] / P - 1.0), <=, 1.0e-14);
      }
      else
      {
        g_assert_cmpfloat (fabs (ln_out[i] - ln_P), <=, 1.0e-14);
        g_assert_cmpfloat (fabs (lin_out[i] / P - 1.0), <=, 1.0e-14);
      }
    }

    g_ptr_array_unref (data_array);
    nc_galaxy_shape_factor_data_unref (data);
    nc_galaxy_shape_factor_free (gsf);
  }

  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);
}

/* A shape factor with restrict-range, and the prepared data of one galaxy at (ra, dec)
 * [deg] with its interpolant built */
static NcGalaxyShapeFactorSpectral *
_restrict_range_build_at (TestNcGalaxyShapeFactorSpectral *test, const gboolean restrict_range, const gdouble ra, const gdouble dec, NcGalaxyShapeFactorData **data)
{
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  GPtrArray *data_array              = g_ptr_array_new ();

  nc_galaxy_shape_factor_spectral_set_restrict_range (gsfs, restrict_range);
  p_data->ra  = ra;
  p_data->dec = dec;

  *data = nc_galaxy_shape_factor_data_new (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (NC_GALAXY_SHAPE_FACTOR (gsfs), *data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
  g_ptr_array_add (data_array, *data);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array (NC_GALAXY_SHAPE_FACTOR (gsfs), test->mset, data_array, TRUE, TRUE));
  nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, *data, 3.0);

  g_ptr_array_unref (data_array);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_redshift_factor_data_unref (z_data);

  return gsfs;
}

/* _restrict_range_build_at() for a galaxy at right ascension ra (degrees from the halo
 * center) on the equator */
static NcGalaxyShapeFactorSpectral *
_restrict_range_build (TestNcGalaxyShapeFactorSpectral *test, const gboolean restrict_range, const gdouble ra, NcGalaxyShapeFactorData **data)
{
  return _restrict_range_build_at (test, restrict_range, ra, 0.0, data);
}

/*
 * restrict-range on a galaxy 0.8 degrees from the halo center, log10MDelta below 14.5: the
 * interpolant covers less than interp-g-max (0.3). At |g| = 0.29 the restricted factor
 * gives the integral, also through eval_nodes, and the unrestricted one the interpolant; at
 * |g| = 0.001 both agree with the integral within interp-abstol. With Omega_c free the
 * bound does not hold, and a new build of the restricted factor covers interp-g-max again
 * and gives the unrestricted value at |g| = 0.29; fixed again, the integral.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  const gdouble g_far[2]   = {0.29 * cos (0.4), 0.29 * sin (0.4)};
  const gdouble g_near[2]  = {0.001, 0.0};
  NcGalaxyShapeFactorData *data_r, *data_f;
  NcGalaxyShapeFactorSpectral *gsfs_r, *gsfs_f;
  guint log10M_idx = 0;
  gdouble ln_direct, ln_r, ln_f;

  g_assert_true (ncm_model_param_index_from_name (hms, "log10MDelta", &log10M_idx, NULL));
  ncm_model_param_set_upper_bound (hms, log10M_idx, 14.5);

  gsfs_r = _restrict_range_build (test, TRUE, 0.8, &data_r);
  gsfs_f = _restrict_range_build (test, FALSE, 0.8, &data_f);
  g_assert_true (nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs_r));
  g_assert_false (nc_galaxy_shape_factor_spectral_get_restrict_range (gsfs_f));

  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data_r->pop_data, g_far[0], g_far[1], 0.3, 0.2, 0.1);
  ln_r      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->pop, data_r, g_far[0], g_far[1], 0.3, 0.2);
  ln_f      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_f), test->pop, data_f, g_far[0], g_far[1], 0.3, 0.2);
  g_assert_cmpfloat (ln_r, ==, ln_direct);
  g_assert_cmpfloat (ln_f, !=, ln_direct);

  /* eval_nodes follows the same range */
  {
    gdouble ln_n;

    NC_GALAXY_SHAPE_FACTOR_GET_CLASS (gsfs_r)->eval_nodes (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->pop, data_r, 1, &g_far[0], &g_far[1], 0.3, 0.2, TRUE, &ln_n);
    g_assert_cmpfloat (ln_n, ==, ln_direct);
  }
  g_assert_cmpfloat (fabs (ln_f - ln_direct), <=, 1.0e-5);

  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data_r->pop_data, g_near[0], g_near[1], 0.3, 0.2, 0.1);
  ln_r      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->pop, data_r, g_near[0], g_near[1], 0.3, 0.2);
  g_assert_cmpfloat (fabs (ln_r - ln_direct), <=, 1.0e-5);

  ncm_mset_param_set_ftype (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->mset);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs_r, test->mset, data_r, 3.0);
  ln_r = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->pop, data_r, g_far[0], g_far[1], 0.3, 0.2);
  g_assert_cmpfloat (ln_r, ==, ln_f);

  ncm_mset_param_set_ftype (test->mset, nc_hicosmo_id (), NC_HICOSMO_DE_OMEGA_C, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
  nc_galaxy_shape_factor_prepare (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->mset);
  nc_galaxy_shape_factor_spectral_build_interp (gsfs_r, test->mset, data_r, 3.0);
  ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs_r, test->pop, data_r->pop_data, g_far[0], g_far[1], 0.3, 0.2, 0.1);
  ln_r      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs_r), test->pop, data_r, g_far[0], g_far[1], 0.3, 0.2);
  g_assert_cmpfloat (ln_r, ==, ln_direct);

  nc_galaxy_shape_factor_data_unref (data_r);
  nc_galaxy_shape_factor_data_unref (data_f);
  nc_galaxy_shape_factor_spectral_free (gsfs_r);
  nc_galaxy_shape_factor_spectral_free (gsfs_f);
}

/* The |g| the interpolant of data covers along the direction of g_far, found by bisection: the
 * interpolant differs from the integral, the integral it falls back to beyond its range is
 * the same bit for bit. */
static gdouble
_restrict_range_covered (TestNcGalaxyShapeFactorSpectral *test, NcGalaxyShapeFactorSpectral *gsfs, NcGalaxyShapeFactorData *data)
{
  gdouble lo = 1.0e-4, hi = 0.299;
  guint i;

  for (i = 0; i < 24; i++)
  {
    const gdouble g         = 0.5 * (lo + hi);
    const gdouble g_1       = g * cos (0.4);
    const gdouble g_2       = g * sin (0.4);
    const gdouble ln_direct = nc_galaxy_shape_factor_spectral_eval_ln_direct (gsfs, test->pop, data->pop_data, g_1, g_2, 0.3, 0.2, 0.1);
    const gdouble ln_r      = nc_galaxy_shape_factor_eval_ln_marginal (NC_GALAXY_SHAPE_FACTOR (gsfs), test->pop, data, g_1, g_2, 0.3, 0.2);

    if (ln_r != ln_direct)
      lo = g;
    else
      hi = g;
  }

  return lo;
}

static gdouble _restrict_range_covered_at (TestNcGalaxyShapeFactorSpectral *test, const gdouble ra, const gdouble dec);

/*
 * restrict-range with a free concentration: the shear is not monotonic in cDelta, so the
 * range of a galaxy is the largest over the box of the concentration. With cDelta free over
 * [1, 10] the interpolant covers at least what it covers at every fixed concentration of the
 * box (a grid of them), and not much more: a restriction given up would cover interp-g-max.
 * The galaxy at 0.8 degrees has its largest shear at the smallest concentration, the one at
 * 0.05 degrees at a larger one. A galaxy at 0.002 degrees with log10MDelta up to 16 reaches
 * kappa >= 1 at the corner for some concentration of the box, and covers interp-g-max.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_cdelta (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  const gdouble ra[2]      = {0.8, 0.05};
  guint log10M_idx         = 0, c_idx = 0;
  guint j, k;

  g_assert_true (ncm_model_param_index_from_name (hms, "log10MDelta", &log10M_idx, NULL));
  g_assert_true (ncm_model_param_index_from_name (hms, "cDelta", &c_idx, NULL));
  ncm_model_param_set_upper_bound (hms, log10M_idx, 14.5);
  ncm_model_param_set_lower_bound (hms, c_idx, 1.0);
  ncm_model_param_set_upper_bound (hms, c_idx, 10.0);

  for (j = 0; j < 2; j++)
  {
    gdouble covered_max = 0.0, covered_free;

    for (k = 0; k < 10; k++)
    {
      ncm_model_param_set (hms, c_idx, 1.0 + k);
      covered_max = MAX (covered_max, _restrict_range_covered_at (test, ra[j], 0.0));
    }

    ncm_model_param_set (hms, c_idx, 4.0);
    ncm_model_param_set_ftype (hms, c_idx, NCM_PARAM_TYPE_FREE);
    ncm_mset_prepare_fparam_map (test->mset);
    covered_free = _restrict_range_covered_at (test, ra[j], 0.0);
    ncm_model_param_set_ftype (hms, c_idx, NCM_PARAM_TYPE_FIXED);
    ncm_mset_prepare_fparam_map (test->mset);

    g_assert_cmpfloat (covered_max, >, 1.0e-3);
    g_assert_cmpfloat (covered_free, >=, covered_max * (1.0 - 1.0e-3));
    g_assert_cmpfloat (covered_free, <, 0.299);
    g_assert_cmpfloat (covered_free, <=, covered_max * 1.1);
  }

  ncm_model_param_set_upper_bound (hms, log10M_idx, 16.0);
  ncm_model_param_set_ftype (hms, c_idx, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);
  g_assert_cmpfloat (_restrict_range_covered_at (test, 0.002, 0.0), >, 0.298);

  ncm_model_param_set_ftype (hms, c_idx, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
}

/* The |g| the interpolant of a galaxy at (ra, dec) covers, built with the models of the mset */
static gdouble
_restrict_range_covered_at (TestNcGalaxyShapeFactorSpectral *test, const gdouble ra, const gdouble dec)
{
  NcGalaxyShapeFactorData *data;
  NcGalaxyShapeFactorSpectral *gsfs = _restrict_range_build_at (test, TRUE, ra, dec, &data);
  const gdouble covered             = _restrict_range_covered (test, gsfs, data);

  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_shape_factor_spectral_free (gsfs);

  return covered;
}

/*
 * restrict-range with a free halo center over the box [-0.5, 0.5]^2 [deg], and a galaxy at
 * (0.6, 0), 0.1 degrees from the box but closer to the center of the box than the largest
 * displacement in it: the interpolant covers less than interp-g-max, and what it covers with
 * the halo center fixed at the nearest point of the box, (0.5, 0). The galaxy at
 * (359.4, 0), across the wrap of the right ascension, and the one at (0, 0.6), inside the
 * right ascensions of the box and 0.1 degrees above it, cover the same.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_centre_box (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp));
  NcmModel *hp             = ncm_mset_peek (test->mset, nc_halo_position_id ());
  gdouble covered_nearest, covered_free, covered_wrap, covered_dec;
  guint log10M_idx = 0;

  g_assert_true (ncm_model_param_index_from_name (hms, "log10MDelta", &log10M_idx, NULL));
  ncm_model_param_set_upper_bound (hms, log10M_idx, 14.5);

  ncm_model_param_set (hp, NC_HALO_POSITION_RA, 0.5);
  ncm_model_param_set (hp, NC_HALO_POSITION_DEC, 0.0);
  covered_nearest = _restrict_range_covered_at (test, 0.6, 0.0);

  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_RA, -0.5);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_RA, 0.5);
  ncm_model_param_set_lower_bound (hp, NC_HALO_POSITION_DEC, -0.5);
  ncm_model_param_set_upper_bound (hp, NC_HALO_POSITION_DEC, 0.5);
  ncm_model_param_set (hp, NC_HALO_POSITION_RA, 0.0);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FREE);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FREE);
  ncm_mset_prepare_fparam_map (test->mset);

  covered_free = _restrict_range_covered_at (test, 0.6, 0.0);
  covered_wrap = _restrict_range_covered_at (test, 359.4, 0.0);
  covered_dec  = _restrict_range_covered_at (test, 0.0, 0.6);

  g_assert_cmpfloat (covered_nearest, >, 1.0e-3);
  g_assert_cmpfloat (covered_free, <, 0.299);
  g_assert_cmpfloat (fabs (covered_free / covered_nearest - 1.0), <, 1.0e-3);
  g_assert_cmpfloat (fabs (covered_wrap / covered_free - 1.0), <, 1.0e-3);
  g_assert_cmpfloat (fabs (covered_dec / covered_free - 1.0), <, 1.0e-3);

  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_RA, NCM_PARAM_TYPE_FIXED);
  ncm_model_param_set_ftype (hp, NC_HALO_POSITION_DEC, NCM_PARAM_TYPE_FIXED);
  ncm_mset_prepare_fparam_map (test->mset);
}

/*
 * restrict-range with a concentration-mass relation (Duffy08, no cDelta): the bound is
 * taken at the concentration of the relation at the top of the mass box, so the interpolant
 * of the galaxy at 0.8 degrees covers what it covers with the parametric profile at that
 * concentration, less than interp-g-max.
 */
void
test_nc_galaxy_shape_factor_spectral_restrict_range_cm_relation (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  NcHaloDensityProfile *dp_param = NC_HALO_DENSITY_PROFILE (ncm_mset_peek (test->mset, nc_halo_density_profile_id ()));
  NcmModel *hms_param            = NCM_MODEL (nc_halo_density_profile_peek_mass_summary (dp_param));
  NcHICosmo *cosmo               = NC_HICOSMO (ncm_mset_peek (test->mset, nc_hicosmo_id ()));
  NcHaloPosition *hp             = NC_HALO_POSITION (ncm_mset_peek (test->mset, nc_halo_position_id ()));
  NcHaloMassSummary *hms         = NC_HALO_MASS_SUMMARY (nc_halo_cm_duffy08_new (NC_HALO_MASS_SUMMARY_MASS_DEF_MEAN, 200.0));
  NcHaloDensityProfile *dp       = NC_HALO_DENSITY_PROFILE (nc_halo_density_profile_nfw_new (hms));
  guint idx                      = 0;
  gdouble covered_param, covered, c_top;

  nc_halo_density_profile_ref (dp_param);
  g_assert_false (ncm_model_param_index_from_name (NCM_MODEL (hms), "cDelta", &idx, NULL));
  g_assert_true (ncm_model_param_index_from_name (NCM_MODEL (hms), "log10MDelta", &idx, NULL));
  ncm_model_param_set (NCM_MODEL (hms), idx, 14.5);
  c_top = nc_halo_mass_summary_concentration (hms, cosmo, nc_halo_position_get_redshift (hp));
  ncm_model_param_set (NCM_MODEL (hms), idx, 14.0);
  ncm_model_param_set_upper_bound (NCM_MODEL (hms), idx, 14.5);

  g_assert_true (ncm_model_param_index_from_name (hms_param, "log10MDelta", &idx, NULL));
  ncm_model_param_set_upper_bound (hms_param, idx, 14.5);
  g_assert_true (ncm_model_param_index_from_name (hms_param, "cDelta", &idx, NULL));
  ncm_model_param_set (hms_param, idx, c_top);
  covered_param = _restrict_range_covered_at (test, 0.8, 0.0);

  ncm_mset_set (test->mset, NCM_MODEL (dp), NULL);
  covered = _restrict_range_covered_at (test, 0.8, 0.0);

  g_assert_cmpfloat (covered_param, >, 2.0e-4);
  g_assert_cmpfloat (covered, <, 0.299);
  g_assert_cmpfloat (fabs (covered / covered_param - 1.0), <, 1.0e-3);

  ncm_mset_set (test->mset, NCM_MODEL (dp_param), NULL);
  nc_halo_density_profile_free (dp_param);
  nc_halo_mass_summary_free (hms);
  nc_halo_density_profile_free (dp);
}

/* nc_galaxy_shape_factor_eval_ln_marginal_array() gives, point by point, what
 * nc_galaxy_shape_factor_eval_ln_marginal() gives, by integration for two observed
 * ellipticities and from the interpolant for the one of the galaxy, for shears inside and
 * outside interp-g-max. */
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

    nc_galaxy_shape_factor_data_set (gsf, data, 0.3, 0.2, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);
    g_ptr_array_add (data_array, data);
    g_assert_true (nc_galaxy_shape_factor_prepare_data_array (gsf, test->mset, data_array, TRUE, TRUE));

    if (j == 1)
    {
      nc_galaxy_shape_factor_spectral_build_interp (gsfs, test->mset, data, 3.0);
      ncm_vector_set_all (e1, 0.3);
      ncm_vector_set_all (e2, 0.2);
    }

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

/*
 * A galaxy with 21 nodes from 0.05 to 1.5, some in front of the lens: four ranges that tile
 * the nodes give the values of the full evaluation exactly, linear and logarithmic, with the
 * line table off; one range call leaves the other entries unchanged. With line-abstol 1e-6
 * each range builds its own line table, so the ranges agree with the full call within twice
 * that tolerance.
 */
void
test_nc_galaxy_shape_factor_spectral_eval_nodes_range (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata)
{
  const guint n                      = 21;
  const guint cuts[5]                = {0, 1, 5, 13, 21};
  NcGalaxyShapeFactorSpectral *gsfs  = nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE);
  NcGalaxyShapeFactor *gsf           = NC_GALAXY_SHAPE_FACTOR (gsfs);
  NcGalaxyPositionFactorData *p_data = nc_galaxy_position_factor_data_new (test->position_factor, test->mset);
  NcGalaxyRedshiftFactorData *z_data = nc_galaxy_redshift_factor_data_new (test->redshift_factor, test->mset);
  NcGalaxyShapeFactorData *data;
  NcmVector *z_nodes   = ncm_vector_new (n);
  NcmVector *full      = ncm_vector_new (n);
  NcmVector *part      = ncm_vector_new (n);
  GPtrArray *data_arr  = g_ptr_array_new ();
  GPtrArray *nodes_arr = g_ptr_array_new ();
  guint i, k, l;

  p_data->ra  = 0.3;
  p_data->dec = 0.0;
  data        = nc_galaxy_shape_factor_data_new (gsf, test->mset, p_data, z_data);
  nc_galaxy_shape_factor_data_set (gsf, data, 0.2, -0.1, 0.1, 0.0, 0.0, 0.0, NC_WL_ELLIPTICITY_FRAME_CELESTIAL);

  for (i = 0; i < n; i++)
    ncm_vector_set (z_nodes, i, 0.05 + 1.45 * i / (n - 1.0));

  g_ptr_array_add (data_arr, data);
  g_ptr_array_add (nodes_arr, z_nodes);
  g_assert_true (nc_galaxy_shape_factor_prepare_data_array_at_nodes (gsf, test->mset, data_arr, nodes_arr, TRUE, TRUE, TRUE));

  for (l = 0; l < 2; l++)
  {
    if (l == 0)
      nc_galaxy_shape_factor_eval_at_nodes (gsf, test->mset, data, z_nodes, full);
    else
      nc_galaxy_shape_factor_eval_ln_at_nodes (gsf, test->mset, data, z_nodes, full);

    ncm_vector_set_all (part, GSL_NAN);

    for (k = 0; k < 4; k++)
    {
      if (l == 0)
        nc_galaxy_shape_factor_eval_at_nodes_range (gsf, test->mset, data, z_nodes, cuts[k], cuts[k + 1], part);
      else
        nc_galaxy_shape_factor_eval_ln_at_nodes_range (gsf, test->mset, data, z_nodes, cuts[k], cuts[k + 1], part);

      for (i = 0; i < n; i++)
      {
        if (i < cuts[k + 1])
          g_assert_cmpfloat (ncm_vector_get (part, i), ==, ncm_vector_get (full, i));
        else
          g_assert_true (gsl_isnan (ncm_vector_get (part, i)));
      }
    }
  }

  nc_galaxy_shape_factor_spectral_set_line_abstol (gsfs, 1.0e-6);
  nc_galaxy_shape_factor_eval_ln_at_nodes (gsf, test->mset, data, z_nodes, full);

  for (k = 0; k < 4; k++)
    nc_galaxy_shape_factor_eval_ln_at_nodes_range (gsf, test->mset, data, z_nodes, cuts[k], cuts[k + 1], part);

  for (i = 0; i < n; i++)
    g_assert_cmpfloat (fabs (ncm_vector_get (part, i) - ncm_vector_get (full, i)), <=, 2.0e-6);

  g_ptr_array_unref (data_arr);
  g_ptr_array_unref (nodes_arr);
  ncm_vector_free (z_nodes);
  ncm_vector_free (full);
  ncm_vector_free (part);
  nc_galaxy_shape_factor_data_unref (data);
  nc_galaxy_redshift_factor_data_unref (z_data);
  nc_galaxy_position_factor_data_unref (p_data);
  nc_galaxy_shape_factor_free (gsf);
}

void
test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside (void)
{
  g_test_trap_subprocess ("/nc/galaxy_shape_factor_spectral/eval_nodes_range/outside/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the range [1, 5) is not inside the 3 nodes*");
}

void
test_nc_galaxy_shape_factor_spectral_eval_nodes_range_outside_subprocess (void)
{
  NcGalaxyShapeFactor *gsf = NC_GALAXY_SHAPE_FACTOR (nc_galaxy_shape_factor_spectral_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE));
  NcmVector *z_nodes       = ncm_vector_new (3);
  NcmVector *out           = ncm_vector_new (3);

  nc_galaxy_shape_factor_eval_at_nodes_range (gsf, NULL, NULL, z_nodes, 1, 5, out);
}

