/***************************************************************************
 *            test_nc_galaxy_shape_factor_quad_common.c
 *
 *  Tue Nov 10 2026
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
 * Tests of NcGalaxyShapeFactorQuad, see test_nc_galaxy_shape_factor_quad_common.h: ln P of
 * the default path and of NcGalaxyShapeFactorQuad:use-lp against values certified by the
 * Arb reference tests/tools/nc_galaxy_shape_marginal_arb.c, the |g| = 1 closed form on both
 * paths, and the properties' defaults and serialization round trip.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */

#include "test_nc_galaxy_shape_factor_quad_common.h"

#include <math.h>
#include <glib-object.h>

/*
 * Largest |ln P - ln P_Arb| accepted at the default reltol (1e-7), about 100 times the
 * largest difference measured on the whole table: 1.2e-8 (default path) and 5.0e-9
 * (use-lp).
 */
#define TEST_QUAD_TOL_PSI (1.0e-6)
#define TEST_QUAD_TOL_LP (5.0e-7)

typedef enum _TestQuadPop
{
  TEST_QUAD_POP_GAUSS,
  TEST_QUAD_POP_BETA,
} TestQuadPop;

/*
 * One certified value: the shear convention, the population (Gauss: sigma = p1; Beta:
 * alpha = p1, beta = p2), std_noise, g, eps_obs, ln P, and whether the unit tier checks it.
 */
typedef struct _TestQuadCell
{
  NcGalaxyWLObsEllipConv conv;
  TestQuadPop pop;
  gdouble p1;
  gdouble p2;
  gdouble std_noise;
  gdouble g1;
  gdouble g2;
  gdouble e1;
  gdouble e2;
  gdouble lnP;
  gboolean unit;
} TestQuadCell;

/* Certified values of tests/tools/make_shape_marginal_arb_cases.py */
static const TestQuadCell test_quad_cells[] = {
  /* BEGIN make_shape_marginal_arb_cases.py */
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_GAUSS, 0.3, 0.0, 0.1, 0.3, 0.0, 0.0, 0.5, -1.402431117754455, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_GAUSS, 0.3, 0.0, 0.01, 0.6, 0.0, 0.5, 0.0, -2.536716442503229, TRUE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_GAUSS, 0.05, 0.0, 0.005, 0.1, 0.0, 0.05, 0.0, -0.31255850401115876, TRUE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_GAUSS, 0.05, 0.0, 0.03, 0.6, 0.0, -0.9, 0.0, -182.8521316909928, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_GAUSS, 0.2, 0.0, 0.05, 0.2, -0.4, 0.3, 0.25, -3.978049030486719, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_GAUSS, 0.3, 0.0, 1.0, 3.0, 0.0, 2.0, 0.0, -2.872477807137722, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_GAUSS, 0.3, 0.0, 0.1, 1.5, 0.0, 1.2, 0.0, -4.056650475069435, TRUE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_BETA, 1.0, 1.0, 0.03, 0.6, 0.0, 0.0, 0.5, -4.003703139471936, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_BETA, 4.0, 1.0, 0.1, 0.6, 0.0, 0.9, 0.0, 0.7997452547544199, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_BETA, 1.0, 1.0, 0.01, 0.3, 0.0, 0.0, 1.2, -206.2067527885707, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_BETA, 4.0, 1.0, 0.3, 1.5, 0.0, -0.5, 0.0, -4.309400768872423, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_BETA, 1.57, 1.65, 0.1, 0.3, 0.0, 0.5, 0.0, 0.5334663856834809, TRUE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_BETA, 1.92, 1.98, 0.03, 0.6, 0.0, 0.0, 0.5, -2.5176649373938598, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_BETA, 1.57, 1.65, 0.05, -0.1, 0.25, 0.2, -0.3, -1.532207192020498, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_GAUSS, 0.05, 0.0, 0.1, 0.005111126056639749, 0.0, 0.0, 0.9, -29.85714990728496, TRUE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, TEST_QUAD_POP_BETA, 1.57, 1.65, 0.1, 0.01, 0.0, 0.0, 0.02, 0.2155439914200594, FALSE},
  {NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE, TEST_QUAD_POP_BETA, 4.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -2.1641169160898808, TRUE},
  /* END make_shape_marginal_arb_cases.py */
};

static TestNcGalaxyShapeFactorQuadMode test_quad_mode;

void test_nc_galaxy_shape_factor_quad_arb_psi (void);
void test_nc_galaxy_shape_factor_quad_arb_lp (void);
void test_nc_galaxy_shape_factor_quad_unit_shear (void);
void test_nc_galaxy_shape_factor_quad_properties (void);

gint
test_nc_galaxy_shape_factor_quad_main (gint argc, gchar *argv[], TestNcGalaxyShapeFactorQuadMode mode)
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  test_quad_mode = mode;

  g_test_add_func ("/nc/galaxy_shape_factor_quad/arb/psi", &test_nc_galaxy_shape_factor_quad_arb_psi);
  g_test_add_func ("/nc/galaxy_shape_factor_quad/arb/lp", &test_nc_galaxy_shape_factor_quad_arb_lp);

  if (mode == TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_UNIT)
  {
    g_test_add_func ("/nc/galaxy_shape_factor_quad/unit_shear", &test_nc_galaxy_shape_factor_quad_unit_shear);
    g_test_add_func ("/nc/galaxy_shape_factor_quad/properties", &test_nc_galaxy_shape_factor_quad_properties);
  }

  g_test_run ();

  return 0;
}

/* Returns a new population of kind @pop with parameters @p1, @p2. */
static NcGalaxyShapePop *
_test_quad_pop_new (const TestQuadPop pop, const gdouble p1, const gdouble p2)
{
  if (pop == TEST_QUAD_POP_GAUSS)
  {
    NcGalaxyShapePopGauss *gauss = nc_galaxy_shape_pop_gauss_new ();

    ncm_model_param_set_by_name (NCM_MODEL (gauss), "sigma", p1, NULL);

    return NC_GALAXY_SHAPE_POP (gauss);
  }
  else
  {
    NcGalaxyShapePopBeta *beta = nc_galaxy_shape_pop_beta_new ();

    ncm_model_param_set_by_name (NCM_MODEL (beta), "alpha", p1, NULL);
    ncm_model_param_set_by_name (NCM_MODEL (beta), "beta", p2, NULL);

    return NC_GALAXY_SHAPE_POP (beta);
  }
}

/* Returns a new Quad for @conv, integrating in log-polar coordinates when @use_lp. */
static NcGalaxyShapeFactorQuad *
_test_quad_new (const NcGalaxyWLObsEllipConv conv, const gboolean use_lp)
{
  return g_object_new (NC_TYPE_GALAXY_SHAPE_FACTOR_QUAD, "ellip-conv", conv, "use-lp", use_lp, NULL);
}

/* Compares ln P with the certified table, the rows marked unit only in the unit tier. */
static void
_test_quad_arb (const gboolean use_lp, const gdouble tol)
{
  guint i;

  for (i = 0; i < G_N_ELEMENTS (test_quad_cells); i++)
  {
    const TestQuadCell *c = &test_quad_cells[i];

    if ((test_quad_mode == TEST_NC_GALAXY_SHAPE_FACTOR_QUAD_UNIT) && !c->unit)
      continue;

    {
      NcGalaxyShapePop *pop        = _test_quad_pop_new (c->pop, c->p1, c->p2);
      NcGalaxyShapeFactorQuad *gsf = _test_quad_new (c->conv, use_lp);
      const gdouble lnP            = log (nc_galaxy_shape_factor_quad_eval_direct (gsf, pop, c->g1, c->g2, c->e1, c->e2, c->std_noise));

      g_assert_cmpfloat_with_epsilon (lnP, c->lnP, tol);

      nc_galaxy_shape_factor_quad_free (gsf);
      nc_galaxy_shape_pop_free (pop);
    }
  }
}

void
test_nc_galaxy_shape_factor_quad_arb_psi (void)
{
  _test_quad_arb (FALSE, TEST_QUAD_TOL_PSI);
}

void
test_nc_galaxy_shape_factor_quad_arb_lp (void)
{
  _test_quad_arb (TRUE, TEST_QUAD_TOL_LP);
}

/*
 * At |g| = 1 both paths return ln N_2 (eps_obs - g) = -|eps_obs - g|^2 / 2 s^2 - ln (2 pi s^2)
 * for any population.
 */
void
test_nc_galaxy_shape_factor_quad_unit_shear (void)
{
  const gdouble g[]     = {1.0, 0.0, 0.0, -1.0};
  const gdouble e[]     = {0.3, -0.2, 0.9, 0.4};
  const gdouble s       = 0.05;
  NcGalaxyShapePop *pop = _test_quad_pop_new (TEST_QUAD_POP_BETA, 1.57, 1.65);
  guint i, j, k;

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
    {
      NcGalaxyShapeFactorQuad *gsf = _test_quad_new ((j == 0) ? NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE : NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET, i == 1);

      for (k = 0; k < G_N_ELEMENTS (g); k += 2)
      {
        const gdouble d2     = gsl_pow_2 (e[k] - g[k]) + gsl_pow_2 (e[k + 1] - g[k + 1]);
        const gdouble ln_ref = -d2 / (2.0 * s * s) - log (2.0 * M_PI * s * s);
        const gdouble lnP    = log (nc_galaxy_shape_factor_quad_eval_direct (gsf, pop, g[k], g[k + 1], e[k], e[k + 1], s));

        g_assert_cmpfloat_with_epsilon (lnP, ln_ref, 1.0e-13 * fabs (ln_ref));
      }

      nc_galaxy_shape_factor_quad_free (gsf);
    }
  }

  nc_galaxy_shape_pop_free (pop);
}

/* The properties' defaults, and their values after a serialization round trip. */
void
test_nc_galaxy_shape_factor_quad_properties (void)
{
  NcGalaxyShapeFactorQuad *gsf = nc_galaxy_shape_factor_quad_new (NC_GALAXY_WL_OBS_ELLIP_CONV_TRACE_DET);
  NcmSerialize *ser            = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcGalaxyShapeFactorQuad *dup;
  gboolean use_lp;
  gdouble reltol;

  g_object_get (gsf, "use-lp", &use_lp, "reltol", &reltol, NULL);
  g_assert_false (use_lp);
  g_assert_cmpfloat (reltol, ==, 1.0e-7);

  g_object_set (gsf, "use-lp", TRUE, "reltol", 1.0e-9, NULL);
  dup = NC_GALAXY_SHAPE_FACTOR_QUAD (ncm_serialize_dup_obj (ser, G_OBJECT (gsf)));

  g_object_get (dup, "use-lp", &use_lp, "reltol", &reltol, NULL);
  g_assert_true (use_lp);
  g_assert_cmpfloat (reltol, ==, 1.0e-9);

  nc_galaxy_shape_factor_quad_free (dup);
  nc_galaxy_shape_factor_quad_free (gsf);
  ncm_serialize_free (ser);
}

