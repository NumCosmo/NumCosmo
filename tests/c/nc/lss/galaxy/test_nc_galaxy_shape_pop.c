/***************************************************************************
 *            test_nc_galaxy_shape_pop.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_nc_galaxy_shape_pop.c
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#undef GSL_RANGE_CHECK_OFF
#endif /* HAVE_CONFIG_H */
#include <numcosmo/numcosmo.h>

#include <math.h>
#include <glib.h>
#include <glib-object.h>
#include <gsl/gsl_integration.h>

typedef enum _TestPopKind
{
  TEST_POP_GAUSS,
  TEST_POP_GAUSS_LOCAL,
  TEST_POP_BETA,
} TestPopKind;

/* A population: sigma (Gauss), e_rms (GaussLocal) or alpha and beta (Beta) */
typedef struct _TestPopCase
{
  const gchar *name;
  TestPopKind kind;
  gdouble p1;
  gdouble p2;
} TestPopCase;

typedef struct _TestNcGalaxyShapePop
{
  NcGalaxyShapePop *pop;
  NcGalaxyShapePopData *data;
} TestNcGalaxyShapePop;

static const TestPopCase test_cases[] = {
  {"gauss/0.05", TEST_POP_GAUSS, 0.05, 0.0},
  {"gauss/0.3", TEST_POP_GAUSS, 0.3, 0.0},
  {"gauss_local/0.25", TEST_POP_GAUSS_LOCAL, 0.25, 0.0},
  {"beta/1_1", TEST_POP_BETA, 1.0, 1.0},
  {"beta/1.57_1.65", TEST_POP_BETA, 1.57, 1.65},
  {"beta/4_1", TEST_POP_BETA, 4.0, 1.0},
  {"beta/0.5_0.5", TEST_POP_BETA, 0.5, 0.5},
};

void test_nc_galaxy_shape_pop_new (TestNcGalaxyShapePop *test, gconstpointer pdata);
void test_nc_galaxy_shape_pop_free (TestNcGalaxyShapePop *test, gconstpointer pdata);

void test_nc_galaxy_shape_pop_exponent_at_one (TestNcGalaxyShapePop *test, gconstpointer pdata);
void test_nc_galaxy_shape_pop_reduced_interior (TestNcGalaxyShapePop *test, gconstpointer pdata);
void test_nc_galaxy_shape_pop_reduced_ends (TestNcGalaxyShapePop *test, gconstpointer pdata);
void test_nc_galaxy_shape_pop_reduced_norm (TestNcGalaxyShapePop *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  guint i;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  for (i = 0; i < G_N_ELEMENTS (test_cases); i++)
  {
    const TestPopCase *tc = &test_cases[i];
    gchar *path;

    path = g_strdup_printf ("/nc/galaxy_shape_pop/%s/exponent_at_one", tc->name);
    g_test_add (path, TestNcGalaxyShapePop, tc, &test_nc_galaxy_shape_pop_new, &test_nc_galaxy_shape_pop_exponent_at_one, &test_nc_galaxy_shape_pop_free);
    g_free (path);

    path = g_strdup_printf ("/nc/galaxy_shape_pop/%s/reduced/interior", tc->name);
    g_test_add (path, TestNcGalaxyShapePop, tc, &test_nc_galaxy_shape_pop_new, &test_nc_galaxy_shape_pop_reduced_interior, &test_nc_galaxy_shape_pop_free);
    g_free (path);

    path = g_strdup_printf ("/nc/galaxy_shape_pop/%s/reduced/ends", tc->name);
    g_test_add (path, TestNcGalaxyShapePop, tc, &test_nc_galaxy_shape_pop_new, &test_nc_galaxy_shape_pop_reduced_ends, &test_nc_galaxy_shape_pop_free);
    g_free (path);

    path = g_strdup_printf ("/nc/galaxy_shape_pop/%s/reduced/norm", tc->name);
    g_test_add (path, TestNcGalaxyShapePop, tc, &test_nc_galaxy_shape_pop_new, &test_nc_galaxy_shape_pop_reduced_norm, &test_nc_galaxy_shape_pop_free);
    g_free (path);
  }

  g_test_run ();

  return 0;
}

/* The integrand of the QAWS rule: the reduced density, without its two power laws */
static gdouble
_reduced (gdouble r, gpointer user_data)
{
  TestNcGalaxyShapePop *test = user_data;

  return nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, r);
}

void
test_nc_galaxy_shape_pop_new (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  const TestPopCase *tc = pdata;

  switch (tc->kind)
  {
    case TEST_POP_GAUSS:
      test->pop  = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_new ());
      test->data = nc_galaxy_shape_pop_data_new (test->pop);
      ncm_model_param_set_by_name (NCM_MODEL (test->pop), "sigma", tc->p1, NULL);
      break;
    case TEST_POP_GAUSS_LOCAL:
      test->pop  = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_gauss_local_new ());
      test->data = nc_galaxy_shape_pop_data_new (test->pop);
      nc_galaxy_shape_pop_gauss_local_data_set (NC_GALAXY_SHAPE_POP_GAUSS_LOCAL (test->pop), test->data, tc->p1);
      break;
    case TEST_POP_BETA:
      test->pop  = NC_GALAXY_SHAPE_POP (nc_galaxy_shape_pop_beta_new ());
      test->data = nc_galaxy_shape_pop_data_new (test->pop);
      ncm_model_param_set_by_name (NCM_MODEL (test->pop), "alpha", tc->p1, NULL);
      ncm_model_param_set_by_name (NCM_MODEL (test->pop), "beta", tc->p2, NULL);
      break;
    default:
      g_assert_not_reached ();
      break;
  }

  nc_galaxy_shape_pop_prepare (test->pop, test->data);
}

void
test_nc_galaxy_shape_pop_free (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  nc_galaxy_shape_pop_data_unref (test->data);
  NCM_TEST_FREE (nc_galaxy_shape_pop_free, test->pop);
}

/* 0 for the truncated Gaussians, beta - 1 for Beta */
void
test_nc_galaxy_shape_pop_exponent_at_one (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  const TestPopCase *tc  = pdata;
  const gdouble expected = (tc->kind == TEST_POP_BETA) ? tc->p2 - 1.0 : 0.0;

  g_assert_cmpfloat (nc_galaxy_shape_pop_exponent_at_one (test->pop), ==, expected);
}

/* r^alpha_o (1 - r)^alpha_1 times the reduced density is the density */
void
test_nc_galaxy_shape_pop_reduced_interior (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  const gdouble radii[] = {1.0e-6, 0.01, 0.3, 0.7, 0.99, 1.0 - 1.0e-6};
  const gdouble alpha_o = nc_galaxy_shape_pop_exponent_at_origin (test->pop);
  const gdouble alpha_1 = nc_galaxy_shape_pop_exponent_at_one (test->pop);
  guint i;

  for (i = 0; i < G_N_ELEMENTS (radii); i++)
  {
    const gdouble r       = radii[i];
    const gdouble p       = nc_galaxy_shape_pop_eval_p (test->pop, test->data, r);
    const gdouble reduced = nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, r);

    ncm_assert_cmpdouble_e (pow (r, alpha_o) * pow (1.0 - r, alpha_1) * reduced, ==, p, 1.0e-13, 0.0);
  }
}

/* Finite and positive at r = 0 and r = 1, and continuous there */
void
test_nc_galaxy_shape_pop_reduced_ends (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  const gdouble at_0 = nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, 0.0);
  const gdouble at_1 = nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, 1.0);

  g_assert_true (gsl_finite (at_0));
  g_assert_true (gsl_finite (at_1));
  g_assert_cmpfloat (at_0, >, 0.0);
  g_assert_cmpfloat (at_1, >, 0.0);

  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, 1.0e-8), ==, at_0, 1.0e-12, 0.0);
  ncm_assert_cmpdouble_e (nc_galaxy_shape_pop_eval_p_reduced (test->pop, test->data, 1.0 - 1.0e-10), ==, at_1, 1.0e-6, 0.0);
}

/*
 * The reduced density against its two power laws integrates to 1, by the QAWS rule of GSL
 * with the exponents of exponent_at_origin and exponent_at_one.
 */
void
test_nc_galaxy_shape_pop_reduced_norm (TestNcGalaxyShapePop *test, gconstpointer pdata)
{
  const gdouble alpha_o           = nc_galaxy_shape_pop_exponent_at_origin (test->pop);
  const gdouble alpha_1           = nc_galaxy_shape_pop_exponent_at_one (test->pop);
  gsl_integration_workspace *ws   = gsl_integration_workspace_alloc (1000);
  gsl_integration_qaws_table *tab = gsl_integration_qaws_table_alloc (alpha_o, alpha_1, 0, 0);
  gsl_function F                  = { &_reduced, test };
  gdouble res, err;

  gsl_integration_qaws (&F, 0.0, 1.0, tab, 0.0, 1.0e-12, 1000, ws, &res, &err);

  ncm_assert_cmpdouble_e (res, ==, 1.0, 1.0e-11, 0.0);

  gsl_integration_qaws_table_free (tab);
  gsl_integration_workspace_free (ws);
}

