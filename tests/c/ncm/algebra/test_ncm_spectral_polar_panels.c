/***************************************************************************
 *            test_ncm_spectral_polar_panels.c
 *
 *  Wed Nov 11 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spectral_polar_panels.c
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

#define TEST_N_POINTS (2000)
#define TEST_N_THREADS (4)

typedef struct _TestNcmSpectralPolarPanels
{
  NcmSpectralPolar *sp;
  NcmSpectralPolarPanels *spp;
  GArray *r_edges;
  GArray *theta_edges;
  NcmRNG *rng;
} TestNcmSpectralPolarPanels;

/*
 * The functions expanded by the tests: kind 0 is exp (r cos (theta) / 2) sin (r sin (theta)
 * + 0.3), kind 1 is 1 + r cos (3 theta) / 2 plus a Gaussian of amplitude A and width s at
 * (cx, cy), kind 2 is 1 + r + r theta / 2, kind 3 is 1 for theta < 1 and 2 otherwise, each
 * times scale. n_values counts the values the build asked for.
 *
 * The Gaussian of kind 1 is used with edges whose starting grids have nodes where it is far
 * above the absolute tolerance. One that is below the tolerance at every node changes no
 * coefficient, and the test accepts the panel without it.
 */
typedef struct _TestFunc
{
  guint kind;
  gdouble scale;
  gdouble cx;
  gdouble cy;
  gdouble s;
  gdouble A;
  guint64 n_values;
} TestFunc;

typedef struct _TestThread
{
  NcmSpectralPolar *sp;
  GArray *r_edges;
  GArray *theta_edges;
  gdouble s;
  gdouble value;
} TestThread;

void test_ncm_spectral_polar_panels_new (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_free (TestNcmSpectralPolarPanels *test, gconstpointer pdata);

void test_ncm_spectral_polar_panels_empty (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_submatrix (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_smooth (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_peak (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_reltol (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_bisect (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_reuse (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_threads (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_outside (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_serialize (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_max_depth (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_max_depth_no_error (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_max_depth_no_error_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_tolerances (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_tolerances_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_negative_radius (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_build_negative_radius_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_few (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_few_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_order (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_order_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_order_theta (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_edges_order_theta_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_no_edges (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_no_edges_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_empty (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_empty_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_size (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_size_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_level (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_level_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_across (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_across_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_order (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_append_order_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_variant_size (TestNcmSpectralPolarPanels *test, gconstpointer pdata);
void test_ncm_spectral_polar_panels_variant_size_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata);

#define TEST_NCM_SPECTRAL_POLAR_PANELS_ADD(name, func)                                    \
        g_test_add ("/ncm/spectral_polar_panels/" name, TestNcmSpectralPolarPanels, NULL, \
                    &test_ncm_spectral_polar_panels_new, &func, &test_ncm_spectral_polar_panels_free)

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("empty", test_ncm_spectral_polar_panels_empty);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append", test_ncm_spectral_polar_panels_append);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/submatrix", test_ncm_spectral_polar_panels_append_submatrix);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/smooth", test_ncm_spectral_polar_panels_build_smooth);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/peak", test_ncm_spectral_polar_panels_build_peak);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/reltol", test_ncm_spectral_polar_panels_build_reltol);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/bisect", test_ncm_spectral_polar_panels_build_bisect);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/reuse", test_ncm_spectral_polar_panels_build_reuse);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/threads", test_ncm_spectral_polar_panels_build_threads);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("outside", test_ncm_spectral_polar_panels_outside);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("serialize", test_ncm_spectral_polar_panels_serialize);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("max_depth", test_ncm_spectral_polar_panels_max_depth);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("max_depth/no_error", test_ncm_spectral_polar_panels_max_depth_no_error);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("max_depth/no_error/subprocess", test_ncm_spectral_polar_panels_max_depth_no_error_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/tolerances", test_ncm_spectral_polar_panels_build_tolerances);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/tolerances/subprocess", test_ncm_spectral_polar_panels_build_tolerances_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/negative_radius", test_ncm_spectral_polar_panels_build_negative_radius);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("build/negative_radius/subprocess", test_ncm_spectral_polar_panels_build_negative_radius_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/few", test_ncm_spectral_polar_panels_edges_few);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/few/subprocess", test_ncm_spectral_polar_panels_edges_few_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/order", test_ncm_spectral_polar_panels_edges_order);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/order/subprocess", test_ncm_spectral_polar_panels_edges_order_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/order_theta", test_ncm_spectral_polar_panels_edges_order_theta);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("edges/order_theta/subprocess", test_ncm_spectral_polar_panels_edges_order_theta_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/no_edges", test_ncm_spectral_polar_panels_append_no_edges);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/no_edges/subprocess", test_ncm_spectral_polar_panels_append_no_edges_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/empty", test_ncm_spectral_polar_panels_append_empty);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/empty/subprocess", test_ncm_spectral_polar_panels_append_empty_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/size", test_ncm_spectral_polar_panels_append_size);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/size/subprocess", test_ncm_spectral_polar_panels_append_size_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/level", test_ncm_spectral_polar_panels_append_level);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/level/subprocess", test_ncm_spectral_polar_panels_append_level_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/across", test_ncm_spectral_polar_panels_append_across);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/across/subprocess", test_ncm_spectral_polar_panels_append_across_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/order", test_ncm_spectral_polar_panels_append_order);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("append/order/subprocess", test_ncm_spectral_polar_panels_append_order_subprocess);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("variant/size", test_ncm_spectral_polar_panels_variant_size);
  TEST_NCM_SPECTRAL_POLAR_PANELS_ADD ("variant/size/subprocess", test_ncm_spectral_polar_panels_variant_size_subprocess);

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

/* The function f at (r, theta) */
static gdouble
_func_at (const TestFunc *f, gdouble r, gdouble theta)
{
  const gdouble x = r * cos (theta);
  const gdouble y = r * sin (theta);

  switch (f->kind)
  {
    case 0:

      return f->scale * exp (0.5 * x) * sin (y + 0.3);

    case 1:
    {
      const gdouble dx = x - f->cx;
      const gdouble dy = y - f->cy;

      return f->scale * (1.0 + 0.5 * r * cos (3.0 * theta) + f->A * exp (-0.5 * (dx * dx + dy * dy) / (f->s * f->s)));
    }
    case 2:

      return f->scale * (1.0 + r + 0.5 * r * theta);

    default:

      return f->scale * ((theta < 1.0) ? 1.0 : 2.0);
  }
}

/* The NcmSpectralPolarF of a TestFunc */
static void
_func (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  TestFunc *f = user_data;
  guint j;

  f->n_values += n;

  for (j = 0; j < n; j++)
  {
    g_assert_cmpfloat (cos_theta[j], ==, cos (theta[j]));
    g_assert_cmpfloat (sin_theta[j], ==, sin (theta[j]));

    fval[j] = _func_at (f, r, theta[j]);
  }
}

/* The largest |spp - f| on TEST_N_POINTS random points of the range of spp */
static gdouble
_max_error (TestNcmSpectralPolarPanels *test, NcmSpectralPolarPanels *spp, const TestFunc *f)
{
  gdouble ra, rb, ta, tb;
  gdouble err = 0.0;
  guint i;

  ncm_spectral_polar_panels_get_range (spp, &ra, &rb, &ta, &tb);

  for (i = 0; i < TEST_N_POINTS; i++)
  {
    const gdouble r     = ncm_rng_uniform_gen (test->rng, ra, rb);
    const gdouble theta = ncm_rng_uniform_gen (test->rng, ta, tb);
    const gdouble d     = fabs (ncm_spectral_polar_panels_eval (spp, r, theta) - _func_at (f, r, theta));

    err = MAX (err, d);
  }

  return err;
}

/* T_k (s) */
static gdouble
_cheb (guint k, gdouble s)
{
  return cos (k * acos (s));
}

/* A new (n_rows) by (n_cols) matrix with entry (l, k) = 1 / (1 + k + 2 l) */
static NcmMatrix *
_coeffs_new (guint n_rows, guint n_cols)
{
  NcmMatrix *coeffs = ncm_matrix_new (n_rows, n_cols);
  guint k, l;

  for (l = 0; l < n_rows; l++)
  {
    for (k = 0; k < n_cols; k++)
      ncm_matrix_set (coeffs, l, k, 1.0 / (1.0 + k + 2.0 * l));
  }

  return coeffs;
}

/* The series of coeffs on the panel [ra, rb] x [ta, tb] at (r, theta), term by term */
static gdouble
_series (NcmMatrix *coeffs, gdouble ra, gdouble rb, gdouble ta, gdouble tb, gdouble r, gdouble theta)
{
  const gdouble sr = (2.0 * r - ra - rb) / (rb - ra);
  const gdouble st = (2.0 * theta - ta - tb) / (tb - ta);
  gdouble sum      = 0.0;
  guint k, l;

  for (l = 0; l < ncm_matrix_nrows (coeffs); l++)
  {
    for (k = 0; k < ncm_matrix_ncols (coeffs); k++)
      sum += ncm_matrix_get (coeffs, l, k) * _cheb (k, sr) * _cheb (l, st);
  }

  return sum;
}

/* Builds kind 1 with a Gaussian of width s in a new panel set and evaluates it at one point */
static gpointer
_thread_build (gpointer data)
{
  TestThread *th              = data;
  TestFunc f                  = { 1, 1.0, 0.4 * cos (0.7), 0.4 * sin (0.7), th->s, 3.0, 0 };
  NcmSpectralPolarPanels *spp = ncm_spectral_polar_panels_new ();
  GError *error               = NULL;

  ncm_spectral_polar_build (th->sp, &_func, &f, th->r_edges, th->theta_edges, 0.0, 1.0e-8, spp, NULL, &error);
  g_assert_no_error (error);

  th->value = ncm_spectral_polar_panels_eval (spp, 0.41, 0.69);
  ncm_spectral_polar_panels_free (spp);

  return NULL;
}

/* Sets the edges r {0, 1, 2} and theta {0, 1, 2} in the panel set of the test */
static void
_set_edges_2x2 (TestNcmSpectralPolarPanels *test)
{
  const gdouble edges[] = {0.0, 1.0, 2.0};

  _set_edges (test->r_edges, edges, G_N_ELEMENTS (edges));
  _set_edges (test->theta_edges, edges, G_N_ELEMENTS (edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);
}

void
test_ncm_spectral_polar_panels_new (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  test->sp          = ncm_spectral_polar_new ();
  test->spp         = ncm_spectral_polar_panels_new ();
  test->r_edges     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->theta_edges = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->rng         = ncm_rng_seeded_new (NULL, 1234);
}

void
test_ncm_spectral_polar_panels_free (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_array_unref (test->r_edges);
  g_array_unref (test->theta_edges);
  ncm_rng_free (test->rng);
  NCM_TEST_FREE (ncm_spectral_polar_panels_free, test->spp);
  NCM_TEST_FREE (ncm_spectral_polar_free, test->sp);
}

void
test_ncm_spectral_polar_panels_empty (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmSerialize *ser           = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralPolarPanels *dup = NCM_SPECTRAL_POLAR_PANELS (ncm_serialize_dup_obj (ser, G_OBJECT (test->spp)));

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, 0);
  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (dup), ==, 0);

  ncm_spectral_polar_panels_free (dup);
  ncm_serialize_free (ser);
}

void
test_ncm_spectral_polar_panels_append (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.5, 1.0, 3.0};
  const gdouble t_edges[] = {-1.0, 0.0, 2.0};
  NcmMatrix *coeffs[5];
  const gdouble panel[5][4] = {
    {0.5, 1.0, -1.0, 0.0},
    {0.5, 1.0, 0.0, 2.0},
    {1.0, 3.0, -1.0, 0.0},
    {1.0, 2.0, 0.0, 2.0},
    {2.0, 3.0, 0.0, 2.0}
  };
  gdouble ra, rb, ta, tb;
  guint i, q;

  coeffs[0] = _coeffs_new (3, 5);
  coeffs[1] = _coeffs_new (2, 2);
  coeffs[2] = _coeffs_new (9, 3);
  coeffs[3] = _coeffs_new (5, 5);
  coeffs[4] = _coeffs_new (3, 9);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);

  for (i = 0; i < 5; i++)
    ncm_spectral_polar_panels_append (test->spp, panel[i][0], panel[i][1], panel[i][2], panel[i][3], coeffs[i]);

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, 5);

  ncm_spectral_polar_panels_get_range (test->spp, &ra, &rb, &ta, &tb);
  g_assert_cmpfloat (ra, ==, 0.5);
  g_assert_cmpfloat (rb, ==, 3.0);
  g_assert_cmpfloat (ta, ==, -1.0);
  g_assert_cmpfloat (tb, ==, 2.0);

  for (i = 0; i < 5; i++)
  {
    guint l_r, l_theta;
    const gdouble *c = ncm_spectral_polar_panels_peek_panel (test->spp, i, &ra, &rb, &ta, &tb, &l_r, &l_theta);
    NcmMatrix *copy  = ncm_spectral_polar_panels_get_panel (test->spp, i, &ra, &rb, &ta, &tb);

    g_assert_cmpfloat (ra, ==, panel[i][0]);
    g_assert_cmpfloat (rb, ==, panel[i][1]);
    g_assert_cmpfloat (ta, ==, panel[i][2]);
    g_assert_cmpfloat (tb, ==, panel[i][3]);
    g_assert_cmpuint ((1u << l_r) + 1, ==, ncm_matrix_ncols (coeffs[i]));
    g_assert_cmpuint ((1u << l_theta) + 1, ==, ncm_matrix_nrows (coeffs[i]));
    g_assert_cmpfloat (c[1], ==, ncm_matrix_get (coeffs[i], 0, 1));
    g_assert_cmpfloat (ncm_matrix_cmp (copy, coeffs[i], 0.0), ==, 0.0);

    /* The first corner and the center of a panel belong to it */
    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, ra, ta), ==, i);
    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 0.5 * (ra + rb), 0.5 * (ta + tb)), ==, i);

    for (q = 0; q < 200; q++)
    {
      const gdouble r     = ncm_rng_uniform_gen (test->rng, ra, rb);
      const gdouble theta = ncm_rng_uniform_gen (test->rng, ta, tb);
      const gdouble ref   = _series (coeffs[i], ra, rb, ta, tb, r, theta);

      g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, r, theta), ==, i);
      ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval (test->spp, r, theta), ==, ref, 1.0e-13, 1.0e-13);
      ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval_on (test->spp, i, r, theta), ==, ref, 1.0e-13, 1.0e-13);
    }

    ncm_matrix_free (copy);
  }

  /* The last edges belong to the last panels, and points outside to the nearest ones */
  g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 3.0, 2.0), ==, 4);
  g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 3.0, -0.5), ==, 2);
  g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 1.5, 2.5), ==, 3);
  g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 0.1, -2.0), ==, 0);
  g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 5.0, 3.0), ==, 4);

  ncm_spectral_polar_panels_reset (test->spp);
  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, 0);

  for (i = 0; i < 5; i++)
    ncm_matrix_free (coeffs[i]);
}

void
test_ncm_spectral_polar_panels_append_submatrix (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, 1.0};
  NcmMatrix *full         = _coeffs_new (7, 8);
  NcmMatrix *sub          = ncm_matrix_get_submatrix (full, 1, 2, 5, 3);
  guint q;

  g_assert_cmpuint (ncm_matrix_tda (sub), !=, ncm_matrix_ncols (sub));

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);
  ncm_spectral_polar_panels_append (test->spp, 0.0, 1.0, 0.0, 1.0, sub);

  for (q = 0; q < 200; q++)
  {
    const gdouble r     = ncm_rng_uniform_gen (test->rng, 0.0, 1.0);
    const gdouble theta = ncm_rng_uniform_gen (test->rng, 0.0, 1.0);

    ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval (test->spp, r, theta), ==,
                            _series (sub, 0.0, 1.0, 0.0, 1.0, r, theta), 1.0e-13, 1.0e-13);
  }

  ncm_matrix_free (sub);
  ncm_matrix_free (full);
}

void
test_ncm_spectral_polar_panels_build_smooth (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.2, 1.5};
  const gdouble t_edges[] = {-1.0, 0.5, 2.0};
  const gdouble abstol    = 1.0e-10;
  TestFunc f              = { 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };
  GError *error           = NULL;
  guint64 n_eval          = 0;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, &n_eval, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, 2);
  g_assert_cmpuint (n_eval, ==, f.n_values);
  g_assert_cmpuint (n_eval, >, 0);
  g_assert_cmpfloat (_max_error (test, test->spp, &f), <=, abstol);
}

void
test_ncm_spectral_polar_panels_build_peak (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble abstol = 1.0e-8;
  TestFunc f           = { 1, 1.0, 0.4 * cos (0.7), 0.4 * sin (0.7), 0.02, 3.0, 0 };
  GError *error        = NULL;
  gdouble ra, rb, ta, tb;

  ncm_spectral_polar_peak_edges (test->sp, 1.0, 0.4, 0.7, 0.02, 0.05, test->r_edges, test->theta_edges);
  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, NULL, &error);
  g_assert_no_error (error);

  ncm_spectral_polar_panels_get_range (test->spp, &ra, &rb, &ta, &tb);
  g_assert_cmpfloat (ra, ==, 0.0);
  g_assert_cmpfloat (rb, ==, 1.0);
  ncm_assert_cmpdouble_e (tb - ta, ==, 2.0 * M_PI, 1.0e-15, 0.0);

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), >=,
                    (test->r_edges->len - 1) * (test->theta_edges->len - 1));
  g_assert_cmpfloat (_max_error (test, test->spp, &f), <=, abstol);
}

void
test_ncm_spectral_polar_panels_build_reltol (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.2, 1.5};
  const gdouble t_edges[] = {-1.0, 2.0};
  const gdouble reltol    = 1.0e-9;
  TestFunc f              = { 0, 1.0e6, 0.0, 0.0, 0.0, 0.0, 0 };
  TestFunc f1             = { 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };
  GError *error           = NULL;
  guint64 n_eval, n_eval1;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, reltol, 0.0, test->spp, &n_eval, &error);
  g_assert_no_error (error);

  /* Largest |F| on the range: exp (0.75) 1e6 */
  g_assert_cmpfloat (_max_error (test, test->spp, &f), <=, reltol * 1.0e6 * exp (0.75));

  /* The relative test does not depend on the scale of F */
  ncm_spectral_polar_build (test->sp, &_func, &f1, test->r_edges, test->theta_edges, reltol, 0.0, test->spp, &n_eval1, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (n_eval1, ==, n_eval);
}

void
test_ncm_spectral_polar_panels_build_bisect (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {-M_PI, 0.0, M_PI};
  const gdouble abstol    = 1.0e-7;
  TestFunc f              = { 1, 1.0, 0.4 * cos (0.7), 0.4 * sin (0.7), 0.05, 3.0, 0 };
  GError *error           = NULL;
  guint n_panels, i;

  ncm_spectral_polar_set_level_cap (test->sp, 4);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, NULL, &error);
  g_assert_no_error (error);

  n_panels = ncm_spectral_polar_panels_get_n_panels (test->spp);
  g_assert_cmpuint (n_panels, >, 2);
  g_assert_cmpfloat (_max_error (test, test->spp, &f), <=, abstol);

  /* Every panel is found from its first corner, its center and its last corner moved inside */
  for (i = 0; i < n_panels; i++)
  {
    gdouble ra, rb, ta, tb;
    guint l_r, l_theta;

    ncm_spectral_polar_panels_peek_panel (test->spp, i, &ra, &rb, &ta, &tb, &l_r, &l_theta);

    g_assert_cmpuint (l_r, <=, 4);
    g_assert_cmpuint (l_theta, <=, 4);
    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, ra, ta), ==, i);
    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, 0.5 * (ra + rb), 0.5 * (ta + tb)), ==, i);
    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (test->spp, rb - 1.0e-9 * (rb - ra), tb - 1.0e-9 * (tb - ta)), ==, i);
  }
}

void
test_ncm_spectral_polar_panels_build_reuse (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.2, 1.5};
  const gdouble t_edges[] = {-1.0, 0.5, 2.0};
  const gdouble abstol    = 1.0e-8;
  TestFunc peak           = { 1, 1.0, 0.4 * cos (0.7), 0.4 * sin (0.7), 0.02, 3.0, 0 };
  TestFunc smooth         = { 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };
  GError *error           = NULL;
  gdouble v_first;
  guint n_first;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &smooth, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, NULL, &error);
  g_assert_no_error (error);

  n_first = ncm_spectral_polar_panels_get_n_panels (test->spp);
  v_first = ncm_spectral_polar_panels_eval (test->spp, 0.7, 0.3);

  ncm_spectral_polar_peak_edges (test->sp, 1.0, 0.4, 0.7, 0.02, 0.05, test->r_edges, test->theta_edges);
  ncm_spectral_polar_build (test->sp, &_func, &peak, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpfloat (_max_error (test, test->spp, &peak), <=, abstol);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &smooth, test->r_edges, test->theta_edges, 0.0, abstol, test->spp, NULL, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, n_first);
  g_assert_cmpfloat (ncm_spectral_polar_panels_eval (test->spp, 0.7, 0.3), ==, v_first);
}

void
test_ncm_spectral_polar_panels_build_threads (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble s[TEST_N_THREADS] = {0.01, 0.02, 0.05, 0.2};
  TestThread th[TEST_N_THREADS];
  GThread *threads[TEST_N_THREADS];
  gdouble serial[TEST_N_THREADS];
  guint i;

  ncm_spectral_polar_peak_edges (test->sp, 1.0, 0.4, 0.7, 0.02, 0.05, test->r_edges, test->theta_edges);

  for (i = 0; i < TEST_N_THREADS; i++)
  {
    th[i].sp          = test->sp;
    th[i].r_edges     = test->r_edges;
    th[i].theta_edges = test->theta_edges;
    th[i].s           = s[i];
    th[i].value       = 0.0;

    _thread_build (&th[i]);
    serial[i]   = th[i].value;
    th[i].value = 0.0;
  }

  for (i = 0; i < TEST_N_THREADS; i++)
    threads[i] = g_thread_new ("build", &_thread_build, &th[i]);

  for (i = 0; i < TEST_N_THREADS; i++)
  {
    g_thread_join (threads[i]);
    g_assert_cmpfloat (th[i].value, ==, serial[i]);
  }
}

void
test_ncm_spectral_polar_panels_outside (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.5, 1.0, 2.0};
  const gdouble t_edges[] = {-1.0, 0.0, 1.0};
  TestFunc f              = { 2, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };
  GError *error           = NULL;

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, 1.0e-12, test->spp, NULL, &error);
  g_assert_no_error (error);

  /* F is a polynomial of degree 1 in each variable, so its series extrapolates it */
  ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval (test->spp, 2.05, 1.05), ==, _func_at (&f, 2.05, 1.05), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval (test->spp, 0.48, -1.05), ==, _func_at (&f, 0.48, -1.05), 1.0e-13, 0.0);
  ncm_assert_cmpdouble_e (ncm_spectral_polar_panels_eval (test->spp, 1.5, 1.05), ==, _func_at (&f, 1.5, 1.05), 1.0e-13, 0.0);
}

void
test_ncm_spectral_polar_panels_serialize (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {-M_PI, 0.0, M_PI};
  NcmSerialize *ser       = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  TestFunc f              = { 1, 1.0, 0.4 * cos (0.7), 0.4 * sin (0.7), 0.05, 3.0, 0 };
  GError *error           = NULL;
  NcmSpectralPolarPanels *dup;
  gdouble ra, rb, ta, tb, dra, drb, dta, dtb;
  guint i;

  /* A level cap of 4 makes the build bisect, so the copy has initial panels with several panels */
  ncm_spectral_polar_set_level_cap (test->sp, 4);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, 1.0e-7, test->spp, NULL, &error);
  g_assert_no_error (error);

  dup = NCM_SPECTRAL_POLAR_PANELS (ncm_serialize_dup_obj (ser, G_OBJECT (test->spp)));

  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (dup), ==, ncm_spectral_polar_panels_get_n_panels (test->spp));
  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (dup), >, 2);

  ncm_spectral_polar_panels_get_range (test->spp, &ra, &rb, &ta, &tb);
  ncm_spectral_polar_panels_get_range (dup, &dra, &drb, &dta, &dtb);
  g_assert_cmpfloat (dra, ==, ra);
  g_assert_cmpfloat (drb, ==, rb);
  g_assert_cmpfloat (dta, ==, ta);
  g_assert_cmpfloat (dtb, ==, tb);

  for (i = 0; i < TEST_N_POINTS; i++)
  {
    const gdouble r     = ncm_rng_uniform_gen (test->rng, ra, rb);
    const gdouble theta = ncm_rng_uniform_gen (test->rng, ta, tb);

    g_assert_cmpuint (ncm_spectral_polar_panels_find_panel (dup, r, theta), ==,
                      ncm_spectral_polar_panels_find_panel (test->spp, r, theta));
    g_assert_cmpfloat (ncm_spectral_polar_panels_eval (dup, r, theta), ==,
                       ncm_spectral_polar_panels_eval (test->spp, r, theta));
  }

  ncm_spectral_polar_panels_free (dup);
  ncm_serialize_free (ser);
}

void
test_ncm_spectral_polar_panels_max_depth (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, M_PI, 2.0 * M_PI};
  TestFunc f              = { 3, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };
  GError *error           = NULL;

  ncm_spectral_polar_set_max_depth (test->sp, 3);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 1.0e-8, 0.0, test->spp, NULL, &error);
  g_assert_error (error, NCM_SPECTRAL_POLAR_ERROR, NCM_SPECTRAL_POLAR_ERROR_MAX_DEPTH);
  g_assert_cmpuint (ncm_spectral_polar_panels_get_n_panels (test->spp), ==, 0);
  g_clear_error (&error);
}

void
test_ncm_spectral_polar_panels_max_depth_no_error (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/max_depth/no_error/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*does not pass the test after 3 bisections*");
}

void
test_ncm_spectral_polar_panels_max_depth_no_error_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, M_PI, 2.0 * M_PI};
  TestFunc f              = { 3, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };

  ncm_spectral_polar_set_max_depth (test->sp, 3);

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 1.0e-8, 0.0, test->spp, NULL, NULL);
}

void
test_ncm_spectral_polar_panels_build_tolerances (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/build/tolerances/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*must be non-negative and not both zero*");
}

void
test_ncm_spectral_polar_panels_build_tolerances_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, 1.0};
  TestFunc f              = { 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, 0.0, test->spp, NULL, NULL);
}

void
test_ncm_spectral_polar_panels_build_negative_radius (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/build/negative_radius/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*ncm_spectral_polar_build: the first radial edge -0.1 is negative*");
}

void
test_ncm_spectral_polar_panels_build_negative_radius_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {-0.1, 1.0};
  const gdouble t_edges[] = {0.0, 1.0};
  TestFunc f              = { 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0 };

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));

  ncm_spectral_polar_build (test->sp, &_func, &f, test->r_edges, test->theta_edges, 0.0, 1.0e-8, test->spp, NULL, NULL);
}

void
test_ncm_spectral_polar_panels_edges_few (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/edges/few/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*two or more of each are needed*");
}

void
test_ncm_spectral_polar_panels_edges_few_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0};
  const gdouble t_edges[] = {0.0, 1.0};

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);
}

void
test_ncm_spectral_polar_panels_edges_order (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/edges/order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*radial edges are not ascending*");
}

void
test_ncm_spectral_polar_panels_edges_order_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 0.6, 0.5, 1.0};
  const gdouble t_edges[] = {0.0, 1.0};

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);
}

void
test_ncm_spectral_polar_panels_edges_order_theta (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/edges/order_theta/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*angular edges are not ascending*");
}

void
test_ncm_spectral_polar_panels_edges_order_theta_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble r_edges[] = {0.0, 1.0};
  const gdouble t_edges[] = {0.0, 1.0, 1.0};

  _set_edges (test->r_edges, r_edges, G_N_ELEMENTS (r_edges));
  _set_edges (test->theta_edges, t_edges, G_N_ELEMENTS (t_edges));
  ncm_spectral_polar_panels_set_edges (test->spp, test->r_edges, test->theta_edges);
}

void
test_ncm_spectral_polar_panels_append_no_edges (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/no_edges/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the edges are not set*");
}

void
test_ncm_spectral_polar_panels_append_no_edges_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmMatrix *coeffs = _coeffs_new (3, 3);

  ncm_spectral_polar_panels_append (test->spp, 0.0, 1.0, 0.0, 1.0, coeffs);
}

void
test_ncm_spectral_polar_panels_append_empty (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/empty/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is empty*");
}

void
test_ncm_spectral_polar_panels_append_empty_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmMatrix *coeffs = _coeffs_new (3, 3);

  _set_edges_2x2 (test);
  ncm_spectral_polar_panels_append (test->spp, 0.5, 0.5, 0.0, 1.0, coeffs);
}

void
test_ncm_spectral_polar_panels_append_size (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*each must be a power of two plus one*");
}

void
test_ncm_spectral_polar_panels_append_size_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmMatrix *coeffs = _coeffs_new (3, 4);

  _set_edges_2x2 (test);
  ncm_spectral_polar_panels_append (test->spp, 0.0, 1.0, 0.0, 1.0, coeffs);
}

void
test_ncm_spectral_polar_panels_append_level (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/level/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*the levels 21 and 1 are above 20*");
}

void
test_ncm_spectral_polar_panels_append_level_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble c[9] = {0.0};

  _set_edges_2x2 (test);
  ncm_spectral_polar_panels_append_data (test->spp, 0.0, 1.0, 0.0, 1.0, 21, 1, c);
}

void
test_ncm_spectral_polar_panels_append_across (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/across/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*is not inside one initial panel*");
}

void
test_ncm_spectral_polar_panels_append_across_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmMatrix *coeffs = _coeffs_new (3, 3);

  _set_edges_2x2 (test);
  ncm_spectral_polar_panels_append (test->spp, 0.5, 1.5, 0.0, 1.0, coeffs);
}

void
test_ncm_spectral_polar_panels_append_order (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/append/order/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*a later initial panel has panels already*");
}

void
test_ncm_spectral_polar_panels_append_order_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  NcmMatrix *coeffs = _coeffs_new (3, 3);

  _set_edges_2x2 (test);
  ncm_spectral_polar_panels_append (test->spp, 1.0, 2.0, 0.0, 1.0, coeffs);
  ncm_spectral_polar_panels_append (test->spp, 0.0, 1.0, 1.0, 2.0, coeffs);
}

void
test_ncm_spectral_polar_panels_variant_size (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  g_test_trap_subprocess ("/ncm/spectral_polar_panels/variant/size/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*panel 0 has 8 coefficients at the levels 1 and 1*");
}

void
test_ncm_spectral_polar_panels_variant_size_subprocess (TestNcmSpectralPolarPanels *test, gconstpointer pdata)
{
  const gdouble edges[2] = {0.0, 1.0};
  const gdouble c[8]     = {0.0};
  GVariantBuilder builder;
  GVariant *panels;

  g_variant_builder_init (&builder, G_VARIANT_TYPE ("a(dddduuad)"));
  g_variant_builder_add (&builder, "(dddduu@ad)", 0.0, 1.0, 0.0, 1.0, 1, 1,
                         g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), c, 8, sizeof (gdouble)));

  panels = g_variant_new ("(@ad@ad@a(dddduuad))",
                          g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), edges, 2, sizeof (gdouble)),
                          g_variant_new_fixed_array (G_VARIANT_TYPE ("d"), edges, 2, sizeof (gdouble)),
                          g_variant_builder_end (&builder));

  g_object_new (NCM_TYPE_SPECTRAL_POLAR_PANELS, "panels", panels, NULL);
}

