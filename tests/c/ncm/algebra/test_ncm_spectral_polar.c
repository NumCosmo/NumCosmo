/***************************************************************************
 *            test_ncm_spectral_polar.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spectral_polar.c
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
#include <gsl/gsl_sf_gamma.h>

typedef struct _TestNcmSpectralPolar
{
  NcmSpectralPolar *sp;
  GArray *r_edges;
  GArray *theta_edges;
} TestNcmSpectralPolar;

/*
 * A Gaussian of width s centered at (cx, cy), normalized to 1 on the plane, and the number
 * of values the integrator asked for.
 */
typedef struct _TestGauss
{
  gdouble cx;
  gdouble cy;
  gdouble s;
  guint64 n_values;
} TestGauss;

/* 1 + r cos (3 theta) / 2 plus a Gaussian of amplitude A and width s at (cx, cy) */
typedef struct _TestSpike
{
  gdouble cx;
  gdouble cy;
  gdouble s;
  gdouble A;
} TestSpike;

/* exp (-r^2 / 2 sigma^2) (1 + cos (theta) / 2) */
typedef struct _TestRadial
{
  gdouble sigma;
  gdouble a;
  gdouble b;
} TestRadial;

typedef struct _TestThread
{
  NcmSpectralPolar *sp;
  gdouble s;
  gdouble result;
} TestThread;

void test_ncm_spectral_polar_new (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_free (TestNcmSpectralPolar *test, gconstpointer pdata);

void test_ncm_spectral_polar_defaults (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_serialize (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_weight (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak_origin (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak_wrong_seed (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak_edges (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_power_laws (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_retest (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_bisect (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_level_cap (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_threads (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_max_depth (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_max_depth_no_error (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_max_depth_no_error_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_levels (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_levels_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_exponents (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_exponents_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_edges_order (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_edges_order_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_edge_ratio (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_edge_ratio_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak_edges_args (TestNcmSpectralPolar *test, gconstpointer pdata);
void test_ncm_spectral_polar_peak_edges_args_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata);

#define TEST_NCM_SPECTRAL_POLAR_ADD(name, func)                              \
        g_test_add ("/ncm/spectral_polar/" name, TestNcmSpectralPolar, NULL, \
                    &test_ncm_spectral_polar_new, &func, &test_ncm_spectral_polar_free)

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  TEST_NCM_SPECTRAL_POLAR_ADD ("defaults", test_ncm_spectral_polar_defaults);
  TEST_NCM_SPECTRAL_POLAR_ADD ("serialize", test_ncm_spectral_polar_serialize);
  TEST_NCM_SPECTRAL_POLAR_ADD ("weight", test_ncm_spectral_polar_weight);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak", test_ncm_spectral_polar_peak);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak/origin", test_ncm_spectral_polar_peak_origin);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak/wrong_seed", test_ncm_spectral_polar_peak_wrong_seed);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak_edges", test_ncm_spectral_polar_peak_edges);
  TEST_NCM_SPECTRAL_POLAR_ADD ("power_laws", test_ncm_spectral_polar_power_laws);
  TEST_NCM_SPECTRAL_POLAR_ADD ("retest", test_ncm_spectral_polar_retest);
  TEST_NCM_SPECTRAL_POLAR_ADD ("panels", test_ncm_spectral_polar_panels);
  TEST_NCM_SPECTRAL_POLAR_ADD ("bisect", test_ncm_spectral_polar_bisect);
  TEST_NCM_SPECTRAL_POLAR_ADD ("level_cap", test_ncm_spectral_polar_level_cap);
  TEST_NCM_SPECTRAL_POLAR_ADD ("threads", test_ncm_spectral_polar_threads);
  TEST_NCM_SPECTRAL_POLAR_ADD ("max_depth", test_ncm_spectral_polar_max_depth);
  TEST_NCM_SPECTRAL_POLAR_ADD ("max_depth/no_error", test_ncm_spectral_polar_max_depth_no_error);
  TEST_NCM_SPECTRAL_POLAR_ADD ("max_depth/no_error/subprocess", test_ncm_spectral_polar_max_depth_no_error_subprocess);
  TEST_NCM_SPECTRAL_POLAR_ADD ("levels", test_ncm_spectral_polar_levels);
  TEST_NCM_SPECTRAL_POLAR_ADD ("levels/subprocess", test_ncm_spectral_polar_levels_subprocess);
  TEST_NCM_SPECTRAL_POLAR_ADD ("exponents", test_ncm_spectral_polar_exponents);
  TEST_NCM_SPECTRAL_POLAR_ADD ("exponents/subprocess", test_ncm_spectral_polar_exponents_subprocess);
  TEST_NCM_SPECTRAL_POLAR_ADD ("edges_order", test_ncm_spectral_polar_edges_order);
  TEST_NCM_SPECTRAL_POLAR_ADD ("edges_order/subprocess", test_ncm_spectral_polar_edges_order_subprocess);
  TEST_NCM_SPECTRAL_POLAR_ADD ("edge_ratio", test_ncm_spectral_polar_edge_ratio);
  TEST_NCM_SPECTRAL_POLAR_ADD ("edge_ratio/subprocess", test_ncm_spectral_polar_edge_ratio_subprocess);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak_edges_args", test_ncm_spectral_polar_peak_edges_args);
  TEST_NCM_SPECTRAL_POLAR_ADD ("peak_edges_args/subprocess", test_ncm_spectral_polar_peak_edges_args_subprocess);

  g_test_run ();

  return 0;
}

/* Replaces the entries of edges by the n values of v */
static void
_set_edges (GArray *edges, const gdouble *v, guint n)
{
  g_array_set_size (edges, 0);
  g_array_append_vals (edges, v, n);
}

/* G = 1 */
static void
_one (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = 1.0;
}

/* The Gaussian of a TestGauss; checks the cosines and sines it is given */
static void
_gauss (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  TestGauss *g = user_data;
  guint j;

  g->n_values += n;

  for (j = 0; j < n; j++)
  {
    const gdouble dx = r * cos_theta[j] - g->cx;
    const gdouble dy = r * sin_theta[j] - g->cy;

    g_assert_cmpfloat (cos_theta[j], ==, cos (theta[j]));
    g_assert_cmpfloat (sin_theta[j], ==, sin (theta[j]));

    fval[j] = exp (-0.5 * (dx * dx + dy * dy) / (g->s * g->s)) / (2.0 * M_PI * g->s * g->s);
  }
}

static void
_spike (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  TestSpike *sk = user_data;
  guint j;

  for (j = 0; j < n; j++)
  {
    const gdouble dx = r * cos_theta[j] - sk->cx;
    const gdouble dy = r * sin_theta[j] - sk->cy;

    fval[j] = 1.0 + 0.5 * r * cos (3.0 * theta[j]) + sk->A * exp (-0.5 * (dx * dx + dy * dy) / (sk->s * sk->s));
  }
}

static void
_radial (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  TestRadial *rd  = user_data;
  const gdouble g = exp (-0.5 * r * r / (rd->sigma * rd->sigma));
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = g * (1.0 + 0.5 * cos_theta[j]);
}

/* The radial factor of _radial, for the one-dimensional reference */
static gdouble
_radial_1d (gdouble r, gpointer user_data)
{
  TestRadial *rd = user_data;

  return exp (-0.5 * r * r / (rd->sigma * rd->sigma));
}

/* 1 below theta = 1 and 2 above: a jump at an angle that no bisection of [0, 2 pi] has as an edge */
static void
_step (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  guint j;

  for (j = 0; j < n; j++)
    fval[j] = (theta[j] < 1.0) ? 1.0 : 2.0;
}

/*
 * The integral of the Gaussian on the unit disk, with weight r, from the edges of
 * ncm_spectral_polar_peak_edges() for a peak at (r_peak, theta_peak) of half widths L_r
 * and L_theta.
 */
static gdouble
_gauss_integral (TestNcmSpectralPolar *test, TestGauss *g, gdouble r_peak, gdouble theta_peak, gdouble L_r, gdouble L_theta, gdouble reltol, guint64 *n_eval)
{
  GError *error = NULL;
  gdouble res;

  ncm_spectral_polar_peak_edges (test->sp, 1.0, r_peak, theta_peak, L_r, L_theta, test->r_edges, test->theta_edges);
  res = ncm_spectral_polar_integrate (test->sp, &_gauss, g, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                      reltol, 0.0, n_eval, NULL, &error);
  g_assert_no_error (error);

  return res;
}

static gpointer
_thread_integrate (gpointer data)
{
  TestThread *tt            = data;
  TestGauss g               = { 0.4 * cos (0.7), 0.4 * sin (0.7), tt->s, 0 };
  TestNcmSpectralPolar test = { tt->sp, g_array_new (FALSE, FALSE, sizeof (gdouble)), g_array_new (FALSE, FALSE, sizeof (gdouble)) };

  tt->result = _gauss_integral (&test, &g, 0.4, 0.7, tt->s, tt->s / 0.4, 1.0e-6, NULL);

  g_array_unref (test.r_edges);
  g_array_unref (test.theta_edges);

  return NULL;
}

void
test_ncm_spectral_polar_new (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  test->sp          = ncm_spectral_polar_new ();
  test->r_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->theta_edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
}

void
test_ncm_spectral_polar_free (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_array_unref (test->r_edges);
  g_array_unref (test->theta_edges);
  NCM_TEST_FREE (ncm_spectral_polar_free, test->sp);
}

void
test_ncm_spectral_polar_defaults (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  NcmSpectralPolar *sp = test->sp;

  g_assert_cmpuint (ncm_spectral_polar_get_level_min (sp), ==, 3);
  g_assert_cmpuint (ncm_spectral_polar_get_level_cap (sp), ==, 7);
  g_assert_cmpuint (ncm_spectral_polar_get_max_depth (sp), ==, 30);
  g_assert_cmpfloat (ncm_spectral_polar_get_edge_ratio (sp), ==, 4.0);

  ncm_spectral_polar_set_level_min (sp, 2);
  ncm_spectral_polar_set_level_cap (sp, 6);
  ncm_spectral_polar_set_max_depth (sp, 12);
  ncm_spectral_polar_set_edge_ratio (sp, 2.5);
  g_assert_cmpuint (ncm_spectral_polar_get_level_min (sp), ==, 2);
  g_assert_cmpuint (ncm_spectral_polar_get_level_cap (sp), ==, 6);
  g_assert_cmpuint (ncm_spectral_polar_get_max_depth (sp), ==, 12);
  g_assert_cmpfloat (ncm_spectral_polar_get_edge_ratio (sp), ==, 2.5);
}

/*
 * A copy made through NcmSerialize has the four properties of the original and gives the
 * same integral and the same number of evaluations.
 */
void
test_ncm_spectral_polar_serialize (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  NcmSerialize *ser         = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  TestGauss g               = { 0.4 * cos (0.7), 0.4 * sin (0.7), 0.02, 0 };
  TestNcmSpectralPolar copy = { NULL, test->r_edges, test->theta_edges };
  guint64 n_eval, n_eval_dup;
  gdouble res, res_dup;

  ncm_spectral_polar_set_level_min (test->sp, 2);
  ncm_spectral_polar_set_level_cap (test->sp, 6);
  ncm_spectral_polar_set_max_depth (test->sp, 12);
  ncm_spectral_polar_set_edge_ratio (test->sp, 2.5);

  copy.sp = NCM_SPECTRAL_POLAR (ncm_serialize_dup_obj (ser, G_OBJECT (test->sp)));

  g_assert_cmpuint (ncm_spectral_polar_get_level_min (copy.sp), ==, 2);
  g_assert_cmpuint (ncm_spectral_polar_get_level_cap (copy.sp), ==, 6);
  g_assert_cmpuint (ncm_spectral_polar_get_max_depth (copy.sp), ==, 12);
  g_assert_cmpfloat (ncm_spectral_polar_get_edge_ratio (copy.sp), ==, 2.5);

  res     = _gauss_integral (test, &g, 0.4, 0.7, 0.02, 0.05, 1.0e-6, &n_eval);
  res_dup = _gauss_integral (&copy, &g, 0.4, 0.7, 0.02, 0.05, 1.0e-6, &n_eval_dup);

  g_assert_cmpfloat (res_dup, ==, res);
  g_assert_cmpuint (n_eval_dup, ==, n_eval);

  ncm_spectral_polar_free (copy.sp);
  ncm_serialize_free (ser);
}

/* G = 1: 2 pi R^(a + b + 1) B (a + 1, b + 1), with both power laws, each alone and none */
void
test_ncm_spectral_polar_weight (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble exponents[][2] = {
    {
      0.0, 0.0
    }, {
      1.0, 0.0
    }, {
      -0.43, 0.65
    }, {
      0.5, -0.5
    }, {
      3.0, 0.98
    }, {
      0.0, -0.9
    }, {
      -0.9, 0.0
    }
  };
  const gdouble radii[] = {1.0, 2.5};
  const gdouble reltol  = 1.0e-10;
  guint i, j;

  for (i = 0; i < G_N_ELEMENTS (exponents); i++)
  {
    for (j = 0; j < G_N_ELEMENTS (radii); j++)
    {
      const gdouble a         = exponents[i][0];
      const gdouble b         = exponents[i][1];
      const gdouble R         = radii[j];
      const gdouble r_edges[] = {0.0, 0.5 * R, R};
      const gdouble t_edges[] = {0.0, 0.5 * M_PI, M_PI, 1.5 * M_PI, 2.0 * M_PI};
      const gdouble exact     = 2.0 * M_PI * pow (R, a + b + 1.0) * gsl_sf_beta (a + 1.0, b + 1.0);
      GError *error           = NULL;
      gdouble res;

      _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
      _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

      res = ncm_spectral_polar_integrate (test->sp, &_one, NULL, R, a, b, test->r_edges, test->theta_edges,
                                          reltol, 0.0, NULL, NULL, &error);
      g_assert_no_error (error);
      g_assert_cmpfloat (fabs (res / exact - 1.0), <=, reltol);
    }
  }
}

/* A Gaussian inside the disk, of widths from 0.005 to 0.05, with the edges of its own peak */
void
test_ncm_spectral_polar_peak (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble widths[] = {0.005, 0.02, 0.05};
  const gdouble reltol   = 1.0e-6;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (widths); i++)
  {
    TestGauss g = { 0.4 * cos (0.7), 0.4 * sin (0.7), widths[i], 0 };
    guint64 n_eval;
    const gdouble res = _gauss_integral (test, &g, 0.4, 0.7, widths[i], widths[i] / 0.4, reltol, &n_eval);

    g_assert_cmpfloat (fabs (res - 1.0), <=, reltol);
    g_assert_cmpuint (n_eval, ==, g.n_values);
  }
}

/* A Gaussian at the origin: angular half width above pi / 4, four equal angular panels */
void
test_ncm_spectral_polar_peak_origin (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-6;
  TestGauss g          = { 0.0, 0.0, 0.05, 0 };
  const gdouble res    = _gauss_integral (test, &g, 0.0, 0.0, 0.05, G_MAXDOUBLE, reltol, NULL);

  g_assert_cmpuint (test->theta_edges->len, ==, 5);
  g_assert_cmpfloat (fabs (res - 1.0), <=, reltol);
}

/* Edges made for a peak 15 widths away and 10 times wider: same accuracy, more evaluations */
void
test_ncm_spectral_polar_peak_wrong_seed (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-6;
  const gdouble s      = 0.01;
  TestGauss g          = { 0.4 * cos (0.7), 0.4 * sin (0.7), s, 0 };
  guint64 n_right, n_wrong;
  const gdouble right = _gauss_integral (test, &g, 0.4, 0.7, s, s / 0.4, reltol, &n_right);
  const gdouble wrong = _gauss_integral (test, &g, 0.55, 0.2, 10.0 * s, 10.0 * s / 0.4, reltol, &n_wrong);

  g_assert_cmpfloat (fabs (right - 1.0), <=, reltol);
  g_assert_cmpfloat (fabs (wrong - 1.0), <=, reltol);
  g_assert_cmpuint (n_wrong, >, n_right);
}

/*
 * The edges cover [0, R] and [theta_peak - pi, theta_peak + pi], ascending, with the peak
 * at the center of a panel in each direction.
 */
void
test_ncm_spectral_polar_peak_edges (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble R = 2.0, r_peak = 0.7, theta_peak = -1.3, L_r = 0.01, L_theta = 0.02;
  gboolean found = FALSE;
  guint i;

  ncm_spectral_polar_peak_edges (test->sp, R, r_peak, theta_peak, L_r, L_theta, test->r_edges, test->theta_edges);

  g_assert_cmpfloat (g_array_index (test->r_edges, gdouble, 0), ==, 0.0);
  g_assert_cmpfloat (g_array_index (test->r_edges, gdouble, test->r_edges->len - 1), ==, R);
  g_assert_cmpfloat (g_array_index (test->theta_edges, gdouble, 0), ==, theta_peak - M_PI);
  g_assert_cmpfloat (g_array_index (test->theta_edges, gdouble, test->theta_edges->len - 1), ==, theta_peak + M_PI);

  for (i = 0; i + 1 < test->r_edges->len; i++)
  {
    const gdouble e0 = g_array_index (test->r_edges, gdouble, i);
    const gdouble e1 = g_array_index (test->r_edges, gdouble, i + 1);

    g_assert_cmpfloat (e0, <, e1);

    if ((e0 == r_peak - L_r) && (e1 == r_peak + L_r))
      found = TRUE;
  }

  g_assert_true (found);
  found = FALSE;

  for (i = 0; i + 1 < test->theta_edges->len; i++)
  {
    const gdouble e0 = g_array_index (test->theta_edges, gdouble, i);
    const gdouble e1 = g_array_index (test->theta_edges, gdouble, i + 1);

    g_assert_cmpfloat (e0, <, e1);

    if ((e0 == theta_peak - L_theta) && (e1 == theta_peak + L_theta))
      found = TRUE;
  }

  g_assert_true (found);
}

/*
 * Both power laws with a function of r and theta, against the QAWS rule of GSL for
 * 2 pi int_0^1 r^a (1 - r)^b exp (-r^2 / 2 sigma^2) dr.
 */
void
test_ncm_spectral_polar_power_laws (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble exponents[][2] = {
    {
      -0.43, 0.65
    }, {
      0.92, 0.98
    }, {
      3.0, 0.0
    }, {
      -0.5, -0.5
    }
  };
  const gdouble r_edges[]       = {0.0, 0.25, 0.5, 1.0};
  const gdouble t_edges[]       = {-M_PI, -0.5 * M_PI, 0.0, 0.5 * M_PI, M_PI};
  const gdouble reltol          = 1.0e-9;
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (1000);
  guint i;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  for (i = 0; i < G_N_ELEMENTS (exponents); i++)
  {
    TestRadial rd                   = { 0.3, exponents[i][0], exponents[i][1] };
    gsl_integration_qaws_table *tab = gsl_integration_qaws_table_alloc (rd.a, rd.b, 0, 0);
    gsl_function F                  = { &_radial_1d, &rd };
    GError *error                   = NULL;
    gdouble ref, ref_err, res;

    gsl_integration_qaws (&F, 0.0, 1.0, tab, 0.0, 1.0e-13, 1000, ws, &ref, &ref_err);
    ref *= 2.0 * M_PI;

    res = ncm_spectral_polar_integrate (test->sp, &_radial, &rd, 1.0, rd.a, rd.b, test->r_edges, test->theta_edges,
                                        reltol, 0.0, NULL, NULL, &error);
    g_assert_no_error (error);
    g_assert_cmpfloat (fabs (res / ref - 1.0), <=, reltol);

    gsl_integration_qaws_table_free (tab);
  }

  gsl_integration_workspace_free (ws);
}

/*
 * A background of integral pi plus a Gaussian of amplitude 1e5 and width 0.001, of integral
 * 0.2 pi, at a corner shared by four panels. The starting grids sample the Gaussian at its
 * maximum, so their sum is many times the integral, and panels of the background accepted
 * with that sum are expanded again with the final one.
 */
void
test_ncm_spectral_polar_retest (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 0.4, 0.7, 1.0};
  const gdouble t_edges[] = {-M_PI, -0.5 * M_PI, 0.0, 0.5 * M_PI, M_PI};
  const gdouble reltol    = 1.0e-6;
  TestSpike sk            = { 0.4, 0.0, 1.0e-3, 1.0e5 };
  const gdouble exact     = M_PI + sk.A * 2.0 * M_PI * sk.s * sk.s;
  GError *error           = NULL;
  gdouble res;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  res = ncm_spectral_polar_integrate (test->sp, &_spike, &sk, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                      reltol, 0.0, NULL, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (fabs (res / exact - 1.0), <=, reltol);
}

/*
 * The accepted panels of the integral of the retest case: their values add up to the
 * integral, their areas to the area of the domain, their levels are between level-min and
 * level-cap, and the same call without the output gives the same integral.
 */
void
test_ncm_spectral_polar_panels (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 0.4, 0.7, 1.0};
  const gdouble t_edges[] = {-M_PI, -0.5 * M_PI, 0.0, 0.5 * M_PI, M_PI};
  const gdouble reltol    = 1.0e-6;
  TestSpike sk            = { 0.4, 0.0, 1.0e-3, 1.0e5 };
  NcmMatrix *panels       = NULL;
  GError *error           = NULL;
  gdouble sum             = 0.0;
  gdouble area            = 0.0;
  gdouble res, res_no_panels;
  guint i;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  res = ncm_spectral_polar_integrate (test->sp, &_spike, &sk, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                      reltol, 0.0, NULL, &panels, &error);
  g_assert_no_error (error);
  g_assert_nonnull (panels);
  g_assert_cmpuint (ncm_matrix_ncols (panels), ==, 8);
  g_assert_cmpuint (ncm_matrix_nrows (panels), >, 12);

  for (i = 0; i < ncm_matrix_nrows (panels); i++)
  {
    const gdouble ra = ncm_matrix_get (panels, i, 0);
    const gdouble rb = ncm_matrix_get (panels, i, 1);
    const gdouble ta = ncm_matrix_get (panels, i, 2);
    const gdouble tb = ncm_matrix_get (panels, i, 3);

    g_assert_cmpfloat (ra, <, rb);
    g_assert_cmpfloat (ta, <, tb);
    g_assert_cmpfloat (ncm_matrix_get (panels, i, 4), >=, ncm_spectral_polar_get_level_min (test->sp));
    g_assert_cmpfloat (ncm_matrix_get (panels, i, 4), <=, ncm_spectral_polar_get_level_cap (test->sp));
    g_assert_cmpfloat (ncm_matrix_get (panels, i, 5), >=, ncm_spectral_polar_get_level_min (test->sp));
    g_assert_cmpfloat (ncm_matrix_get (panels, i, 5), <=, ncm_spectral_polar_get_level_cap (test->sp));
    g_assert_cmpfloat (ncm_matrix_get (panels, i, 7), >=, 0.0);

    sum  += ncm_matrix_get (panels, i, 6);
    area += (rb - ra) * (tb - ta);
  }

  ncm_assert_cmpdouble_e (sum, ==, res, 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (area, ==, 2.0 * M_PI, 1.0e-13, 0.0);

  res_no_panels = ncm_spectral_polar_integrate (test->sp, &_spike, &sk, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                                reltol, 0.0, NULL, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (res_no_panels, ==, res);

  ncm_matrix_free (panels);
}

/* Two panels for a Gaussian of width 0.01 and level-cap 4: bisections in r and in theta */
void
test_ncm_spectral_polar_bisect (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {-M_PI, 0.0, M_PI};
  const gdouble reltol    = 1.0e-6;
  TestGauss g             = { 0.4 * cos (0.7), 0.4 * sin (0.7), 0.01, 0 };
  GError *error           = NULL;
  gdouble res;

  ncm_spectral_polar_set_level_cap (test->sp, 4);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  res = ncm_spectral_polar_integrate (test->sp, &_gauss, &g, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                      reltol, 0.0, NULL, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (fabs (res - 1.0), <=, reltol);
}

/* The same integral with level-cap 7, 5 and 9 and level-min 3 and 2 on one object */
void
test_ncm_spectral_polar_level_cap (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const guint caps[]   = {7, 5, 9, 4};
  const guint mins[]   = {3, 3, 2, 2};
  const gdouble reltol = 1.0e-6;
  guint i;

  for (i = 0; i < G_N_ELEMENTS (caps); i++)
  {
    TestGauss g = { 0.4 * cos (0.7), 0.4 * sin (0.7), 0.02, 0 };
    gdouble res;

    ncm_spectral_polar_set_level_cap (test->sp, caps[i]);
    ncm_spectral_polar_set_level_min (test->sp, mins[i]);

    res = _gauss_integral (test, &g, 0.4, 0.7, 0.02, 0.05, reltol, NULL);
    g_assert_cmpfloat (fabs (res - 1.0), <=, reltol);
  }
}

/* Four concurrent calls on one object give the values of the serial calls */
void
test_ncm_spectral_polar_threads (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble widths[4] = {0.005, 0.01, 0.02, 0.05};
  TestThread tt[4];
  GThread *threads[4];
  gdouble serial[4];
  guint t;

  for (t = 0; t < 4; t++)
  {
    tt[t].sp     = test->sp;
    tt[t].s      = widths[t];
    tt[t].result = 0.0;

    _thread_integrate (&tt[t]);
    serial[t] = tt[t].result;
  }

  for (t = 0; t < 4; t++)
    threads[t] = g_thread_new ("integrate", _thread_integrate, &tt[t]);

  for (t = 0; t < 4; t++)
  {
    g_thread_join (threads[t]);
    g_assert_cmpfloat (tt[t].result, ==, serial[t]);
  }
}

/* A jump in theta with max-depth 3: the error, NaN, and no panels */
void
test_ncm_spectral_polar_max_depth (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, M_PI, 2.0 * M_PI};
  NcmMatrix *panels       = (NcmMatrix *) test; /* any non-NULL value: the call must reset it */
  GError *error           = NULL;
  gdouble res;

  ncm_spectral_polar_set_max_depth (test->sp, 3);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  res = ncm_spectral_polar_integrate (test->sp, &_step, NULL, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                      1.0e-8, 0.0, NULL, &panels, &error);
  g_assert_error (error, NCM_SPECTRAL_POLAR_ERROR, NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH);
  g_assert_null (panels);
  g_assert_cmpstr (g_quark_to_string (NCM_SPECTRAL_POLAR_ERROR), ==, "ncm-spectral-polar-error");
  g_assert_true (gsl_isnan (res));
  g_clear_error (&error);
}

void
test_ncm_spectral_polar_max_depth_no_error (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/max_depth/no_error/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not pass the test after 3 bisections*");
}

void
test_ncm_spectral_polar_max_depth_no_error_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, M_PI, 2.0 * M_PI};

  ncm_spectral_polar_set_max_depth (test->sp, 3);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_integrate (test->sp, &_step, NULL, 1.0, 1.0, 0.0, test->r_edges, test->theta_edges,
                                1.0e-8, 0.0, NULL, NULL, NULL);
}

void
test_ncm_spectral_polar_levels (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/levels/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*not below level-cap*");
}

void
test_ncm_spectral_polar_levels_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble edges[] = {0.0, 1.0};

  ncm_spectral_polar_set_level_min (test->sp, 5);
  ncm_spectral_polar_set_level_cap (test->sp, 5);

  _set_edges (test->r_edges, edges, G_N_ELEMENTS (edges));
  _set_edges (test->theta_edges, edges, G_N_ELEMENTS (edges));

  ncm_spectral_polar_integrate (test->sp, &_one, NULL, 1.0, 0.0, 0.0, test->r_edges, test->theta_edges,
                                1.0e-8, 0.0, NULL, NULL, NULL);
}

void
test_ncm_spectral_polar_exponents (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/exponents/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*must be above -1*");
}

void
test_ncm_spectral_polar_exponents_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble edges[] = {0.0, 1.0};

  _set_edges (test->r_edges, edges, G_N_ELEMENTS (edges));
  _set_edges (test->theta_edges, edges, G_N_ELEMENTS (edges));

  ncm_spectral_polar_integrate (test->sp, &_one, NULL, 1.0, -1.0, 0.0, test->r_edges, test->theta_edges,
                                1.0e-8, 0.0, NULL, NULL, NULL);
}

void
test_ncm_spectral_polar_edges_order (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/edges_order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*radial edges are not ascending*");
}

void
test_ncm_spectral_polar_edges_order_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 0.6, 0.5, 1.0};
  const gdouble t_edges[] = {0.0, 1.0};

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_integrate (test->sp, &_one, NULL, 1.0, 0.0, 0.0, test->r_edges, test->theta_edges,
                                1.0e-8, 0.0, NULL, NULL, NULL);
}

void
test_ncm_spectral_polar_edge_ratio (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/edge_ratio/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is not above 1*");
}

void
test_ncm_spectral_polar_edge_ratio_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  ncm_spectral_polar_set_edge_ratio (test->sp, 1.0);
}

void
test_ncm_spectral_polar_peak_edges_args (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar/peak_edges_args/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is not inside*");
}

void
test_ncm_spectral_polar_peak_edges_args_subprocess (TestNcmSpectralPolar *test, gconstpointer pdata)
{
  ncm_spectral_polar_peak_edges (test->sp, 1.0, 1.5, 0.0, 0.1, 0.1, test->r_edges, test->theta_edges);
}

