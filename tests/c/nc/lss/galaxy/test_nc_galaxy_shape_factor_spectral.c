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
#include <glib.h>
#include <glib-object.h>

/*
 * Largest |ln P - ln P_Arb| accepted at the default reltol (1e-5). The largest difference
 * measured is 6.9e-9, on the 20,055 certified points of the certification grid with
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
 * Points with 1 - |g| from 1e-5 to 1e-8, certified by the Arb reference with --coords=uv
 * (radius of ln P below 1e-15): HWL16a-002 rows 475 and 7139 with the Beta population
 * fitted to that catalog. Below 1 - |g| = 1e-5 the class forms the value from integrals at
 * that distance. The largest difference of the class measured on them is 1.9e-9.
 */
static const TestShapeArbCell test_spectral_near_unit_cells[] = {
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591169317498968, 0.28297475387508314, 0.691024972903032, 0.0, -18.15200444413578233, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591264271024746, 0.28297755535316127, 0.691024972903032, 0.0, -18.15200449764241692, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591169317498968, 0.28297475387508314, 0.691024972903032, 0.0, -18.14997346866772457, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.0587911494076252, 0.9591264271024746, 0.28297755535316127, 0.691024972903032, 0.0, -18.15199473594852941, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.999999, 0.0, 0.178728848695755, 0.3698130249977112, -59.93454784752150307, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.99999999, 0.0, 0.178728848695755, 0.3698130249977112, -79.77604456776660732, FALSE, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_SHAPE_ARB_POP_BETA, 1.5709751684077675, 1.6542629053559708, 0.04034872725605965, 0.99999999, 0.0, 0.178728848695755, 0.3698130249977112, -47.75380931482856171, FALSE, FALSE},
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
void test_nc_galaxy_shape_factor_spectral_interp (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_rebuild (TestNcGalaxyShapeFactorSpectral *test, gconstpointer pdata);
void test_nc_galaxy_shape_factor_spectral_interp_g_max (void);
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
  g_test_add ("/nc/galaxy_shape_factor_spectral/hooks", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_hooks,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp,
              &test_nc_galaxy_shape_factor_spectral_free);
  g_test_add ("/nc/galaxy_shape_factor_spectral/interp/rebuild", TestNcGalaxyShapeFactorSpectral, NULL,
              &test_nc_galaxy_shape_factor_spectral_new,
              &test_nc_galaxy_shape_factor_spectral_interp_rebuild,
              &test_nc_galaxy_shape_factor_spectral_free);
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

