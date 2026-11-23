/***************************************************************************
 *            test_ncm_spectral_tensor_builder.c
 *
 *  Thu Nov 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spectral_tensor_builder.c
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
#include <string.h>
#include <glib.h>
#include <glib-object.h>

#define TEST_N_POINTS (2000)
#define TEST_N_THREADS (4)
#define TEST_MAX_N_VARS (3)

/*
 * The functions tabulated by the tests, of x_0, ..., x_{n-1}: kind 0 is
 * exp (0.5 x_0 cos (x_1)) sin (x_0 sin (x_1) + 0.3), the function of
 * test_ncm_spectral_polar_panels.c in the variables (r, theta), with x_1 = 0 when n = 1;
 * kind 1 is exp (sum_i (i + 1) x_i / 4) cos (sum_i x_i); kind 2 is a Gaussian of width
 * 0.05 at x_i = 0.3 plus 1; kind 3 is 1 for x_0 < 0.4 and 2 otherwise; kind 4 is T_4 of
 * x_0 on [0, 1] rounded to its values 1 and -1 at the five Chebyshev-Lobatto points, exactly
 * 1, -1, 1, -1, 1 there; kind 5 is -1.6e4 exp (sum_i (x_i - 1)), with |F| up to 1.6e4.
 * n_values counts the values asked for.
 */
typedef struct _TestFunc
{
  guint kind;
  guint64 n_values;
} TestFunc;

typedef struct _TestNcmSpectralTensorBuilder
{
  NcmSpectralTensorBuilder *builder;
  NcmRNG *rng;
} TestNcmSpectralTensorBuilder;

typedef struct _TestThread
{
  NcmSpectralTensorBuilder *builder;
  NcmSpectralTensor *st;
} TestThread;

/* A test that aborts in a subprocess with a message that matches stderr */
typedef struct _TestTrap
{
  const gchar *name;

  void (*func) (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);

  const gchar *stderr_pattern;
} TestTrap;

void test_ncm_spectral_tensor_builder_new (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_free (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);

void test_ncm_spectral_tensor_builder_properties (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_polar (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_accuracy (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_large_values (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_trim (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_fixed (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_levels (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_reuse (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_threads (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_max_depth (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_predict (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);

void test_ncm_spectral_tensor_builder_traps (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_tolerances_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_trim_frac_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_no_domain_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_bad_levels_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_all_fixed_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);
void test_ncm_spectral_tensor_builder_max_depth_no_error_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata);

static const TestTrap traps[] = {
  {"tolerances", &test_ncm_spectral_tensor_builder_tolerances_subprocess, "*must be non-negative and not both zero*"},
  {"trim_frac", &test_ncm_spectral_tensor_builder_trim_frac_subprocess, "*ncm_spectral_tensor_builder_set_trim_frac: 1 is outside [0, 1)*"},
  {"no_domain", &test_ncm_spectral_tensor_builder_no_domain_subprocess, "*the domain of the table is not set*"},
  {"bad_levels", &test_ncm_spectral_tensor_builder_bad_levels_subprocess, "*the levels 4 and 3 of variable 0 must satisfy*"},
  {"all_fixed", &test_ncm_spectral_tensor_builder_all_fixed_subprocess, "*every variable has fixed nodes*"},
  {"max_depth_no_error", &test_ncm_spectral_tensor_builder_max_depth_no_error_subprocess, "*the subdomain*does not pass the test after 2 bisections*"},
};

#define TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD(name, func)                                      \
        g_test_add ("/ncm/spectral_tensor_builder/" name, TestNcmSpectralTensorBuilder, NULL, \
                    &test_ncm_spectral_tensor_builder_new, &func, &test_ncm_spectral_tensor_builder_free)

gint
main (gint argc, gchar *argv[])
{
  guint t;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("properties", test_ncm_spectral_tensor_builder_properties);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("polar", test_ncm_spectral_tensor_builder_polar);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("accuracy", test_ncm_spectral_tensor_builder_accuracy);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("large_values", test_ncm_spectral_tensor_builder_large_values);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("trim", test_ncm_spectral_tensor_builder_trim);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("fixed", test_ncm_spectral_tensor_builder_fixed);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("levels", test_ncm_spectral_tensor_builder_levels);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("reuse", test_ncm_spectral_tensor_builder_reuse);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("threads", test_ncm_spectral_tensor_builder_threads);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("max_depth", test_ncm_spectral_tensor_builder_max_depth);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("predict", test_ncm_spectral_tensor_builder_predict);
  TEST_NCM_SPECTRAL_TENSOR_BUILDER_ADD ("traps", test_ncm_spectral_tensor_builder_traps);

  for (t = 0; t < G_N_ELEMENTS (traps); t++)
  {
    gchar *path = g_strdup_printf ("/ncm/spectral_tensor_builder/traps/%s/subprocess", traps[t].name);

    g_test_add (path, TestNcmSpectralTensorBuilder, NULL,
                &test_ncm_spectral_tensor_builder_new, traps[t].func, &test_ncm_spectral_tensor_builder_free);
    g_free (path);
  }

  g_test_run ();

  return 0;
}

/* The function f at the point x of n_vars values */
static gdouble
_func_at (const TestFunc *f, guint n_vars, const gdouble *x)
{
  switch (f->kind)
  {
    case 0:
    {
      const gdouble r = x[0];
      const gdouble t = (n_vars > 1) ? x[1] : 0.0;

      return exp (0.5 * (r * cos (t))) * sin (r * sin (t) + 0.3);
    }
    case 1:
    {
      gdouble e = 0.0;
      gdouble s = 0.0;
      guint i;

      for (i = 0; i < n_vars; i++)
      {
        e += 0.25 * (i + 1) * x[i];
        s += x[i];
      }

      return exp (e) * cos (s);
    }
    case 2:
    {
      gdouble d2 = 0.0;
      guint i;

      for (i = 0; i < n_vars; i++)
        d2 += gsl_pow_2 (x[i] - 0.3);

      return 1.0 + exp (-0.5 * d2 / gsl_pow_2 (0.05));
    }
    case 3:

      return (x[0] < 0.4) ? 1.0 : 2.0;

    case 4:

      return (lround (4.0 * acos (CLAMP (2.0 * x[0] - 1.0, -1.0, 1.0)) / M_PI) % 2 == 0) ? 1.0 : -1.0;

    case 5:
    {
      gdouble e = 0.0;
      guint i;

      for (i = 0; i < n_vars; i++)
        e += x[i] - 1.0;

      return -1.6e4 * exp (e);
    }
    default: /* LCOV_EXCL_LINE */
      g_assert_not_reached ();

      return 0.0;
  }
}

/* The NcmSpectralTensorF of f: its values on the tensor product, the last variable fastest */
static void
_func (guint n_vars, const guint *n, const gdouble * const *x, gdouble *fval, gpointer user_data)
{
  TestFunc *f = user_data;
  guint k[TEST_MAX_N_VARS];
  gdouble p[TEST_MAX_N_VARS] = {0.0};
  gsize total                = 1;
  gsize q;
  guint i;

  for (i = 0; i < n_vars; i++)
  {
    k[i]   = 0;
    total *= n[i];
  }

  for (q = 0; q < total; q++)
  {
    for (i = 0; i < n_vars; i++)
      p[i] = x[i][k[i]];

    fval[q] = _func_at (f, n_vars, p);

    for (i = n_vars; i > 0; i--)
    {
      if (++k[i - 1] < n[i - 1])
        break;

      k[i - 1] = 0;
    }
  }

  f->n_values += total;
}

/* The NcmSpectralPolarF of f, kind 0, with the cosines and sines of theta computed as _func_at() does */
static void
_polar_func (gpointer user_data, gdouble r, guint n, const gdouble *theta, const gdouble *cos_theta, const gdouble *sin_theta, gdouble *fval)
{
  TestFunc *f = user_data;
  guint j;

  for (j = 0; j < n; j++)
  {
    const gdouble x[2] = {r, theta[j]};

    fval[j] = _func_at (f, 2, x);
  }

  f->n_values += n;
}

/*
 * Sets the domain of st from the first and the last of the edges e[i][0..n_e[i] - 1] of
 * each variable and divides it as the grid of those edges, leaving every leaf empty.
 */
static void
_table_set (NcmSpectralTensor *st, guint n_vars, const gdouble * const *e, const guint *n_e)
{
  gdouble lower[TEST_MAX_N_VARS];
  gdouble upper[TEST_MAX_N_VARS];
  guint i;

  for (i = 0; i < n_vars; i++)
  {
    lower[i] = e[i][0];
    upper[i] = e[i][n_e[i] - 1];
  }

  ncm_spectral_tensor_set_domain (st, n_vars, lower, n_vars, upper);

  for (i = 0; i < n_vars; i++)
    ncm_spectral_tensor_split_empty (st, i, n_e[i], e[i]);
}

/* A table of n_vars variables divided as the grid of the edges e[i][0..n_e[i] - 1] */
static NcmSpectralTensor *
_table_new (guint n_vars, const gdouble * const *e, const guint *n_e)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (n_vars);

  _table_set (st, n_vars, e, n_e);

  return st;
}

/* The total number of coefficients stored in st */
static gsize
_n_coeffs (NcmSpectralTensor *st)
{
  const guint n_vars = ncm_spectral_tensor_get_n_vars (st);
  gsize total        = 0;
  guint p, i;

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
  {
    const gdouble *lower, *upper;
    const guint *K;
    gsize len = 1;

    ncm_spectral_tensor_peek_panel (st, p, &lower, &upper, &K);

    for (i = 0; i < n_vars; i++)
      len *= K[i];

    total += len;
  }

  return total;
}

/* Asserts that the two tables hold the same panels, bit for bit */
static void
_assert_same_tables (NcmSpectralTensor *st1, NcmSpectralTensor *st2)
{
  const guint n_vars = ncm_spectral_tensor_get_n_vars (st1);
  guint p, i;

  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st1), ==, ncm_spectral_tensor_get_n_panels (st2));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st1), ==, ncm_spectral_tensor_get_n_nodes (st2));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st1), ==, 0);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st2), ==, 0);

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st1); p++)
  {
    const gdouble *lower1, *upper1, *lower2, *upper2;
    const guint *K1, *K2;
    const gdouble *c1 = ncm_spectral_tensor_peek_panel (st1, p, &lower1, &upper1, &K1);
    const gdouble *c2 = ncm_spectral_tensor_peek_panel (st2, p, &lower2, &upper2, &K2);
    gsize len         = 1;

    for (i = 0; i < n_vars; i++)
      len *= K1[i];

    g_assert_cmpmem (lower1, sizeof (gdouble) * n_vars, lower2, sizeof (gdouble) * n_vars);
    g_assert_cmpmem (upper1, sizeof (gdouble) * n_vars, upper2, sizeof (gdouble) * n_vars);
    g_assert_cmpmem (K1, sizeof (guint) * n_vars, K2, sizeof (guint) * n_vars);
    g_assert_cmpmem (c1, sizeof (gdouble) * len, c2, sizeof (gdouble) * len);
  }
}

/* The largest |table - f| at n random points of the range of st */
static gdouble
_max_error (TestNcmSpectralTensorBuilder *test, NcmSpectralTensor *st, const TestFunc *f, guint n)
{
  const guint n_vars = ncm_spectral_tensor_get_n_vars (st);
  gdouble err        = 0.0;
  gdouble x[TEST_MAX_N_VARS];
  guint q, i;

  for (q = 0; q < n; q++)
  {
    for (i = 0; i < n_vars; i++)
    {
      gdouble a, b;

      ncm_spectral_tensor_get_range (st, i, &a, &b);
      x[i] = ncm_rng_uniform_gen (test->rng, a, b);
    }

    err = GSL_MAX (err, fabs (ncm_spectral_tensor_eval (st, n_vars, x) - _func_at (f, n_vars, x)));
  }

  return err;
}

void
test_ncm_spectral_tensor_builder_new (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  test->builder = ncm_spectral_tensor_builder_new ();
  test->rng     = ncm_rng_seeded_new (NULL, 2027);
}

void
test_ncm_spectral_tensor_builder_free (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_spectral_tensor_builder_free, test->builder);
  ncm_rng_free (test->rng);
}

void
test_ncm_spectral_tensor_builder_properties (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  g_assert_cmpuint (ncm_spectral_tensor_builder_get_level_min (test->builder), ==, 3);
  g_assert_cmpuint (ncm_spectral_tensor_builder_get_level_cap (test->builder), ==, 7);
  g_assert_cmpuint (ncm_spectral_tensor_builder_get_max_depth (test->builder), ==, 30);
  g_assert_cmpfloat (ncm_spectral_tensor_builder_get_trim_frac (test->builder), ==, 0.0);

  ncm_spectral_tensor_builder_set_level_min (test->builder, 4);
  ncm_spectral_tensor_builder_set_level_cap (test->builder, 8);
  ncm_spectral_tensor_builder_set_max_depth (test->builder, 12);
  ncm_spectral_tensor_builder_set_trim_frac (test->builder, 1.0e-4);

  g_assert_cmpuint (ncm_spectral_tensor_builder_get_level_min (test->builder), ==, 4);
  g_assert_cmpuint (ncm_spectral_tensor_builder_get_level_cap (test->builder), ==, 8);
  g_assert_cmpuint (ncm_spectral_tensor_builder_get_max_depth (test->builder), ==, 12);
  g_assert_cmpfloat (ncm_spectral_tensor_builder_get_trim_frac (test->builder), ==, 1.0e-4);

  {
    guint level_min, level_cap, max_depth;
    gdouble trim_frac;

    g_object_get (test->builder, "level-min", &level_min, "level-cap", &level_cap, "max-depth", &max_depth, "trim-frac", &trim_frac, NULL);

    g_assert_cmpuint (level_min, ==, 4);
    g_assert_cmpuint (level_cap, ==, 8);
    g_assert_cmpuint (max_depth, ==, 12);
    g_assert_cmpfloat (trim_frac, ==, 1.0e-4);
  }
}

/*
 * With two variables and no trim the builder repeats ncm_spectral_polar_build(): the same
 * panels, levels and number of values, and the same coefficients up to 4 ulp of the largest
 * one of each panel. The coefficients are not bit-identical: GCC computes part of the node
 * table of ncm_spectral_polar.c with the vector cosine of libmvec, 1 to 2 ulp from cos().
 */
void
test_ncm_spectral_tensor_builder_polar (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble r_edges[]     = {0.0, 0.5, 1.2};
  const gdouble t_edges[]     = {0.0, 1.0, 2.5, 3.1};
  const gdouble *e[2]         = {r_edges, t_edges};
  const guint n_e[2]          = {G_N_ELEMENTS (r_edges), G_N_ELEMENTS (t_edges)};
  const gdouble abstol[]      = {1.0e-6, 1.0e-10};
  NcmSpectralPolar *sp        = ncm_spectral_polar_new ();
  NcmSpectralPolarPanels *spp = ncm_spectral_polar_panels_new ();
  GArray *r_arr               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *t_arr               = g_array_new (FALSE, FALSE, sizeof (gdouble));
  guint a;

  g_array_append_vals (r_arr, r_edges, G_N_ELEMENTS (r_edges));
  g_array_append_vals (t_arr, t_edges, G_N_ELEMENTS (t_edges));

  for (a = 0; a < G_N_ELEMENTS (abstol); a++)
  {
    NcmSpectralTensor *st = _table_new (2, e, n_e);
    TestFunc f_polar      = {0, 0};
    TestFunc f_tensor     = {0, 0};
    GError *error         = NULL;
    guint64 n_polar, n_tensor;
    guint p;

    ncm_spectral_polar_build (sp, &_polar_func, &f_polar, r_arr, t_arr, 0.0, abstol[a], spp, &n_polar, &error);
    g_assert_no_error (error);
    ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, abstol[a], &f_tensor, &n_tensor, &error);
    g_assert_no_error (error);

    g_assert_cmpuint (n_tensor, ==, n_polar);
    g_assert_cmpuint (f_tensor.n_values, ==, f_polar.n_values);
    g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), ==, ncm_spectral_polar_panels_get_n_panels (spp));
    g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), >=, (n_e[0] - 1) * (n_e[1] - 1));

    for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
    {
      const gdouble *lower, *upper;
      const guint *K;
      const gdouble *c = ncm_spectral_tensor_peek_panel (st, p, &lower, &upper, &K);
      gdouble ra, rb, ta, tb;
      guint l_r, l_t;
      const gdouble *cp = ncm_spectral_polar_panels_peek_panel (spp, p, &ra, &rb, &ta, &tb, &l_r, &l_t);

      g_assert_cmpfloat (lower[0], ==, ra);
      g_assert_cmpfloat (upper[0], ==, rb);
      g_assert_cmpfloat (lower[1], ==, ta);
      g_assert_cmpfloat (upper[1], ==, tb);
      g_assert_cmpuint (K[0], ==, (1u << l_r) + 1);
      g_assert_cmpuint (K[1], ==, (1u << l_t) + 1);
      {
        gdouble cmax = 0.0;
        gdouble dmax = 0.0;
        guint q;

        for (q = 0; q < K[0] * K[1]; q++)
        {
          cmax = GSL_MAX (cmax, fabs (cp[q]));
          dmax = GSL_MAX (dmax, fabs (c[q] - cp[q]));
        }

        g_assert_cmpfloat (dmax, <=, 4.0 * GSL_DBL_EPSILON * cmax);
      }
    }

    ncm_spectral_tensor_free (st);
  }

  g_array_unref (r_arr);
  g_array_unref (t_arr);
  ncm_spectral_polar_panels_free (spp);
  ncm_spectral_polar_free (sp);
}

/*
 * The table agrees with F within the tolerance at random points, for one, two and three
 * variables and for the smooth function and the peak.
 */
void
test_ncm_spectral_tensor_builder_accuracy (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]   = {-1.0, 0.0, 1.0};
  const gdouble e1[]   = {0.0, 0.7};
  const gdouble e2[]   = {-0.5, 0.5, 1.0};
  const gdouble *e[3]  = {e0, e1, e2};
  const guint n_e[3]   = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1), G_N_ELEMENTS (e2)};
  const guint kinds[]  = {1, 2};
  const gdouble abstol = 1.0e-8;
  guint n_vars, k;

  for (n_vars = 1; n_vars <= TEST_MAX_N_VARS; n_vars++)
  {
    for (k = 0; k < G_N_ELEMENTS (kinds); k++)
    {
      NcmSpectralTensor *st = _table_new (n_vars, e, n_e);
      TestFunc f            = {kinds[k], 0};
      GError *error         = NULL;
      guint64 n_eval;

      ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, abstol, &f, &n_eval, &error);
      g_assert_no_error (error);
      g_assert_cmpuint (n_eval, ==, f.n_values);
      g_assert_cmpfloat (_max_error (test, st, &f, TEST_N_POINTS), <=, abstol);

      ncm_spectral_tensor_free (st);
    }
  }
}

/*
 * A function with |F| up to 1.6e4, kind 5, on a domain with no splits, for 1 to 3
 * variables. The tolerance 1e-12 is below the spacing of the doubles at 1.6e4, 1.8e-12. As a
 * relative tolerance it is multiplied by V, the largest |F| on the panel: the build ends and
 * the error at random points is at most reltol V. As an absolute tolerance it is not
 * reached: the build ends in NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH.
 */
void
test_ncm_spectral_tensor_builder_large_values (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {-1.0, 1.0};
  const gdouble e1[]    = {0.0, 0.7};
  const gdouble e2[]    = {-0.5, 1.0};
  const gdouble *e[3]   = {e0, e1, e2};
  const guint n_e[3]    = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1), G_N_ELEMENTS (e2)};
  const gdouble upper[] = {1.0, 0.7, 1.0};
  const gdouble tol     = 1.0e-12;
  guint n_vars;

  for (n_vars = 1; n_vars <= TEST_MAX_N_VARS; n_vars++)
  {
    NcmSpectralTensor *st = _table_new (n_vars, e, n_e);
    TestFunc f            = {5, 0};
    const gdouble V       = fabs (_func_at (&f, n_vars, upper));
    GError *error         = NULL;

    ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, tol, 0.0, &f, NULL, &error);
    g_assert_no_error (error);
    g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), ==, 1);
    g_assert_cmpfloat (_max_error (test, st, &f, TEST_N_POINTS), <=, tol * V);

    _table_set (st, n_vars, e, n_e);
    ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, tol, &f, NULL, &error);
    g_assert_error (error, NCM_SPECTRAL_TENSOR_BUILDER_ERROR, NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH);
    g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), ==, 0);
    g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, 1);

    g_clear_error (&error);
    ncm_spectral_tensor_free (st);
  }
}

/*
 * With trim-frac f the table stores fewer coefficients, differs from the full one by at
 * most f times the tolerance at random points, and stays within the tolerance of F.
 */
void
test_ncm_spectral_tensor_builder_trim (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]     = {-1.0, 0.0, 1.0};
  const gdouble e1[]     = {0.0, 0.7};
  const gdouble *e[2]    = {e0, e1};
  const guint n_e[2]     = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  const gdouble abstol   = 1.0e-6;
  const gdouble frac     = 0.5;
  NcmSpectralTensor *st  = _table_new (2, e, n_e);
  NcmSpectralTensor *stt = _table_new (2, e, n_e);
  TestFunc f             = {1, 0};
  GError *error          = NULL;
  gdouble x[2];
  guint q;

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, abstol, &f, NULL, &error);
  g_assert_no_error (error);

  ncm_spectral_tensor_builder_set_trim_frac (test->builder, frac);
  ncm_spectral_tensor_builder_build (test->builder, &_func, stt, NULL, NULL, 0.0, abstol, &f, NULL, &error);
  g_assert_no_error (error);

  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (stt), ==, ncm_spectral_tensor_get_n_panels (st));
  g_assert_cmpuint (_n_coeffs (stt), <, _n_coeffs (st));

  for (q = 0; q < TEST_N_POINTS; q++)
  {
    x[0] = ncm_rng_uniform_gen (test->rng, e0[0], e0[2]);
    x[1] = ncm_rng_uniform_gen (test->rng, e1[0], e1[1]);

    g_assert_cmpfloat (fabs (ncm_spectral_tensor_eval (stt, 2, x) - ncm_spectral_tensor_eval (st, 2, x)), <=, frac * abstol);
  }

  g_assert_cmpfloat (_max_error (test, stt, &f, TEST_N_POINTS), <=, abstol);

  ncm_spectral_tensor_free (st);
  ncm_spectral_tensor_free (stt);
}

/*
 * A variable with level-min = level-cap keeps its number of points: every panel has 2^l + 1
 * coefficients in it and spans one interval of the grid in it, and the table is still
 * within the tolerance when that number of points resolves F.
 */
void
test_ncm_spectral_tensor_builder_fixed (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {-1.0, 0.0, 1.0};
  const gdouble e1[]    = {0.0, 0.35, 0.7};
  const gdouble *e[2]   = {e0, e1};
  const guint n_e[2]    = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  const guint l_min[2]  = {3, 5};
  const guint l_cap[2]  = {7, 5};
  const gdouble abstol  = 1.0e-8;
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {1, 0};
  GError *error         = NULL;
  guint p;

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, l_min, l_cap, 0.0, abstol, &f, NULL, &error);
  g_assert_no_error (error);

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
  {
    const gdouble *lower, *upper;
    const guint *K;

    ncm_spectral_tensor_peek_panel (st, p, &lower, &upper, &K);

    g_assert_cmpuint (K[1], ==, 33);
    g_assert_true (((lower[1] == e1[0]) && (upper[1] == e1[1])) || ((lower[1] == e1[1]) && (upper[1] == e1[2])));
  }

  g_assert_cmpfloat (_max_error (test, st, &f, TEST_N_POINTS), <=, abstol);

  ncm_spectral_tensor_free (st);
}

/* Per-variable levels: every panel has at least 2^level-min + 1 coefficients in each variable */
void
test_ncm_spectral_tensor_builder_levels (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {-1.0, 1.0};
  const gdouble e1[]    = {0.0, 0.7};
  const gdouble *e[2]   = {e0, e1};
  const guint n_e[2]    = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  const guint l_min[2]  = {2, 5};
  const guint l_cap[2]  = {8, 6};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {1, 0};
  GError *error         = NULL;
  guint p;

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, l_min, l_cap, 0.0, 1.0e-10, &f, NULL, &error);
  g_assert_no_error (error);

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
  {
    const gdouble *lower, *upper;
    const guint *K;

    ncm_spectral_tensor_peek_panel (st, p, &lower, &upper, &K);

    g_assert_cmpuint (K[0], >=, 5);
    g_assert_cmpuint (K[0], <=, 257);
    g_assert_cmpuint (K[1], >=, 33);
    g_assert_cmpuint (K[1], <=, 65);
  }

  ncm_spectral_tensor_free (st);
}

/* A table set again and built again gives the same table; a build on a full table changes nothing */
void
test_ncm_spectral_tensor_builder_reuse (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]     = {-1.0, 0.0, 1.0};
  const gdouble e1[]     = {0.0, 0.7};
  const gdouble *e[2]    = {e0, e1};
  const guint n_e[2]     = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  NcmSpectralTensor *st  = _table_new (2, e, n_e);
  NcmSpectralTensor *st2 = _table_new (2, e, n_e);
  TestFunc f             = {2, 0};
  GError *error          = NULL;

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);
  g_assert_no_error (error);
  ncm_spectral_tensor_builder_build (test->builder, &_func, st2, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);
  g_assert_no_error (error);
  _table_set (st, 2, e, n_e);
  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);
  g_assert_no_error (error);

  _assert_same_tables (st, st2);

  f.n_values = 0;
  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);
  g_assert_no_error (error);
  g_assert_cmpuint (f.n_values, ==, 0);

  _assert_same_tables (st, st2);

  ncm_spectral_tensor_free (st);
  ncm_spectral_tensor_free (st2);
}

static gpointer
_thread_build (gpointer data)
{
  TestThread *th = data;
  TestFunc f     = {2, 0};
  GError *error  = NULL;

  ncm_spectral_tensor_builder_build (th->builder, &_func, th->st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);
  g_assert_no_error (error);

  return NULL;
}

/* Builds in several threads with one builder give the table of a build in one thread */
void
test_ncm_spectral_tensor_builder_threads (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {-1.0, 0.0, 1.0};
  const gdouble e1[]    = {0.0, 0.7};
  const gdouble *e[2]   = {e0, e1};
  const guint n_e[2]    = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestThread th[TEST_N_THREADS];
  GThread *threads[TEST_N_THREADS];
  TestThread serial = {test->builder, st};
  guint t;

  _thread_build (&serial);

  for (t = 0; t < TEST_N_THREADS; t++)
  {
    th[t].builder = test->builder;
    th[t].st      = _table_new (2, e, n_e);
    threads[t]    = g_thread_new ("build", &_thread_build, &th[t]);
  }

  for (t = 0; t < TEST_N_THREADS; t++)
  {
    g_thread_join (threads[t]);
    _assert_same_tables (st, th[t].st);
    ncm_spectral_tensor_free (th[t].st);
  }

  ncm_spectral_tensor_free (st);
}

/* A discontinuity exhausts max-depth: the error, and the table left with its domain and an empty root */
void
test_ncm_spectral_tensor_builder_max_depth (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {0.0, 1.0};
  const gdouble e1[]    = {0.0, 1.0};
  const gdouble *e[2]   = {e0, e1};
  const guint n_e[2]    = {G_N_ELEMENTS (e0), G_N_ELEMENTS (e1)};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {3, 0};
  GError *error         = NULL;
  gdouble a, b;

  ncm_spectral_tensor_builder_set_max_depth (test->builder, 2);
  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, &error);

  g_assert_error (error, NCM_SPECTRAL_TENSOR_BUILDER_ERROR, NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), ==, 0);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, 1);
  ncm_spectral_tensor_get_range (st, 0, &a, &b);
  g_assert_cmpfloat (a, ==, 0.0);
  g_assert_cmpfloat (b, ==, 1.0);

  g_clear_error (&error);
  ncm_spectral_tensor_free (st);
}

/*
 * The two ends of the prediction of the next change, one level below the level cap.
 *
 * Kind 3 does not depend on x_1, so the coefficients above index 0 in that variable are
 * exactly zero: with x_0 fixed and x_1 at level 3 of 4 the predicted change is 0, the level
 * is raised instead of the node being split, and the build still fails on the step in x_0.
 *
 * Kind 4 at level 2 in x_0 has the coefficient of T_4 only: the band (1, 2] is empty and
 * the band (2, 4] is not, so no decay is measured, the predicted change is infinite and the
 * node is split at once, which max-depth 0 does not allow.
 */
void
test_ncm_spectral_tensor_builder_predict (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]  = {0.0, 1.0};
  const gdouble *e[2] = {e0, e0};
  const guint n_e[2]  = {2, 2};

  {
    const guint l_min[2]  = {3, 3};
    const guint l_cap[2]  = {3, 4};
    NcmSpectralTensor *st = _table_new (2, e, n_e);
    TestFunc f            = {3, 0};
    GError *error         = NULL;

    ncm_spectral_tensor_builder_set_max_depth (test->builder, 1);
    ncm_spectral_tensor_builder_build (test->builder, &_func, st, l_min, l_cap, 0.0, 1.0e-8, &f, NULL, &error);

    g_assert_error (error, NCM_SPECTRAL_TENSOR_BUILDER_ERROR, NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH);
    g_assert_nonnull (strstr (error->message, "after 1 bisections"));
    /* Both levels of x_1 are sampled on the root and on its first child, where the build stops */
    g_assert_cmpuint (f.n_values, ==, 2 * 9 * 17);

    g_clear_error (&error);
    ncm_spectral_tensor_free (st);
  }

  {
    const guint l_min[2]  = {2, 2};
    const guint l_cap[2]  = {3, 2};
    NcmSpectralTensor *st = _table_new (2, e, n_e);
    TestFunc f            = {4, 0};
    GError *error         = NULL;

    ncm_spectral_tensor_builder_set_max_depth (test->builder, 0);
    ncm_spectral_tensor_builder_build (test->builder, &_func, st, l_min, l_cap, 0.0, 1.0e-8, &f, NULL, &error);

    g_assert_error (error, NCM_SPECTRAL_TENSOR_BUILDER_ERROR, NCM_SPECTRAL_TENSOR_BUILDER_ERROR_MAX_DEPTH);
    g_assert_nonnull (strstr (error->message, "after 0 bisections"));
    /* The node is given up on its starting grid, without the level 3 of x_0 */
    g_assert_cmpuint (f.n_values, ==, 5 * 5);

    g_clear_error (&error);
    ncm_spectral_tensor_free (st);
  }
}

void
test_ncm_spectral_tensor_builder_traps (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  guint t;

  for (t = 0; t < G_N_ELEMENTS (traps); t++)
  {
    gchar *path = g_strdup_printf ("/ncm/spectral_tensor_builder/traps/%s/subprocess", traps[t].name);

    g_test_trap_subprocess (path, 0, 0);
    g_test_trap_assert_failed ();
    g_test_trap_assert_stderr (traps[t].stderr_pattern);
    g_free (path);
  }
}

/*
 * The functions below abort in a subprocess of the traps test, so the lines after the
 * aborting call never run.
 */
/* LCOV_EXCL_START */

void
test_ncm_spectral_tensor_builder_tolerances_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {0.0, 1.0};
  const gdouble *e[2]   = {e0, e0};
  const guint n_e[2]    = {2, 2};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {1, 0};

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 0.0, &f, NULL, NULL);
}

void
test_ncm_spectral_tensor_builder_trim_frac_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  ncm_spectral_tensor_builder_set_trim_frac (test->builder, 1.0);
}

void
test_ncm_spectral_tensor_builder_no_domain_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (2);
  TestFunc f            = {1, 0};

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, NULL);
}

void
test_ncm_spectral_tensor_builder_bad_levels_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {0.0, 1.0};
  const gdouble *e[2]   = {e0, e0};
  const guint n_e[2]    = {2, 2};
  const guint l_min[2]  = {4, 3};
  const guint l_cap[2]  = {3, 7};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {1, 0};

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, l_min, l_cap, 0.0, 1.0e-8, &f, NULL, NULL);
}

void
test_ncm_spectral_tensor_builder_all_fixed_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {0.0, 1.0};
  const gdouble *e[2]   = {e0, e0};
  const guint n_e[2]    = {2, 2};
  const guint l[2]      = {3, 4};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {1, 0};

  ncm_spectral_tensor_builder_build (test->builder, &_func, st, l, l, 0.0, 1.0e-8, &f, NULL, NULL);
}

void
test_ncm_spectral_tensor_builder_max_depth_no_error_subprocess (TestNcmSpectralTensorBuilder *test, gconstpointer pdata)
{
  const gdouble e0[]    = {0.0, 1.0};
  const gdouble *e[2]   = {e0, e0};
  const guint n_e[2]    = {2, 2};
  NcmSpectralTensor *st = _table_new (2, e, n_e);
  TestFunc f            = {3, 0};

  ncm_spectral_tensor_builder_set_max_depth (test->builder, 2);
  ncm_spectral_tensor_builder_build (test->builder, &_func, st, NULL, NULL, 0.0, 1.0e-8, &f, NULL, NULL);
}

/* LCOV_EXCL_STOP */

