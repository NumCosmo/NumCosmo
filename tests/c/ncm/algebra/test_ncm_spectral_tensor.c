/***************************************************************************
 *            test_ncm_spectral_tensor.c
 *
 *  Thu Nov 19 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_spectral_tensor.c
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

#define TEST_N_POINTS (1000)
#define TEST_K_MAX (6)
#define TEST_MAX_N_VARS (5)

/*
 * A table filled by hand and the panels it was filled with. The domain is the product of
 * [a_i, b_i], divided as a grid at m_i in every variable, 2^n subdomains. The subdomain
 * that starts at (a_0, ..., a_{n-1}) is split at its midpoint in the last variable, and its
 * lower half is split again at its midpoint in the first variable, so the table has
 * 2^n + 2 panels. Each panel has random numbers of coefficients in [1, TEST_K_MAX] and
 * random coefficients in [-1, 1].
 */
typedef struct _TestNcmSpectralTensor
{
  guint n_vars;
  NcmSpectralTensor *st;
  NcmRNG *rng;
  gdouble a[TEST_MAX_N_VARS];
  gdouble m[TEST_MAX_N_VARS];
  gdouble b[TEST_MAX_N_VARS];
  GArray *node;   /* guint, per panel */
  GArray *lower;  /* gdouble, n_vars per panel */
  GArray *upper;  /* gdouble, n_vars per panel */
  GArray *K;      /* guint, n_vars per panel */
  GArray *offset; /* gsize, per panel */
  GArray *c;      /* gdouble, the coefficients of every panel */
} TestNcmSpectralTensor;

/* A test that aborts in a subprocess with a message that matches stderr */
typedef struct _TestTrap
{
  const gchar *name;

  void (*func) (TestNcmSpectralTensor *test, gconstpointer pdata);

  const gchar *stderr_pattern;
} TestTrap;

void test_ncm_spectral_tensor_new (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_free (TestNcmSpectralTensor *test, gconstpointer pdata);

void test_ncm_spectral_tensor_basic (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_tree (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_eval (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_eval_long (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_eval_rows (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_eval_array (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_edges_rule (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_serialize (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_reset (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_domain (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_incomplete (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_empty (TestNcmSpectralTensor *test, gconstpointer pdata);

void test_ncm_spectral_tensor_traps (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_domain_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_domain_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_node_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_not_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_variable_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_no_edges_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_outside_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_order_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_empty_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_split_empty_order_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_panel_not_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_panel_n_K_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_panel_n_coeffs_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_set_panel_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_find_panel_n_vars_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_eval_array_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_find_panel_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_find_panel_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_node_not_split_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_node_not_panel_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_variant_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);
void test_ncm_spectral_tensor_variant_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata);

static const TestTrap traps[] = {
  {"set_domain/length", &test_ncm_spectral_tensor_set_domain_length_subprocess, "*1 lower ends and 2 upper ends for 2 variables*"},
  {"set_domain/empty", &test_ncm_spectral_tensor_set_domain_empty_subprocess, "*of variable 1 is empty*"},
  {"split/node", &test_ncm_spectral_tensor_split_node_subprocess, "*ncm_spectral_tensor_split: node 1000, the table has*"},
  {"split/not_empty", &test_ncm_spectral_tensor_split_not_empty_subprocess, "*ncm_spectral_tensor_split: node 0 is not empty*"},
  {"split/variable", &test_ncm_spectral_tensor_split_variable_subprocess, "*variable 2, the table has 2 variables*"},
  {"split/no_edges", &test_ncm_spectral_tensor_split_no_edges_subprocess, "*a split needs one interior edge or more*"},
  {"split/outside", &test_ncm_spectral_tensor_split_outside_subprocess, "*must be ascending and strictly inside*edge 0 is*"},
  {"split/order", &test_ncm_spectral_tensor_split_order_subprocess, "*must be ascending and strictly inside*edge 1 is*"},
  {"split_empty/no_domain", &test_ncm_spectral_tensor_split_empty_no_domain_subprocess, "*ncm_spectral_tensor_split_empty: the domain is not set*"},
  {"split_empty/order", &test_ncm_spectral_tensor_split_empty_order_subprocess, "*the edges of variable 0 are not ascending at index 0*"},
  {"set_panel/not_empty", &test_ncm_spectral_tensor_set_panel_not_empty_subprocess, "*ncm_spectral_tensor_set_panel: node 0 is not empty*"},
  {"set_panel/n_K", &test_ncm_spectral_tensor_set_panel_n_K_subprocess, "*3 numbers of coefficients for 2 variables*"},
  {"set_panel/n_coeffs", &test_ncm_spectral_tensor_set_panel_n_coeffs_subprocess, "*0 coefficients in variable 1, outside*"},
  {"set_panel/length", &test_ncm_spectral_tensor_set_panel_length_subprocess, "*5 coefficients, the numbers of coefficients give 6*"},
  {"find_panel/n_vars", &test_ncm_spectral_tensor_find_panel_n_vars_subprocess, "*a point of 3 values for 2 variables*"},
  {"eval_array/length", &test_ncm_spectral_tensor_eval_array_length_subprocess, "*ncm_spectral_tensor_eval_array: 3 values are not points of 2 variables*"},
  {"find_panel/no_domain", &test_ncm_spectral_tensor_find_panel_no_domain_subprocess, "*ncm_spectral_tensor_find_panel: the domain is not set*"},
  {"find_panel/empty", &test_ncm_spectral_tensor_find_panel_empty_subprocess, "*the point is in the empty node 0*"},
  {"node/not_split", &test_ncm_spectral_tensor_node_not_split_subprocess, "*is not a split*"},
  {"node/not_panel", &test_ncm_spectral_tensor_node_not_panel_subprocess, "*node 0 is not a panel*"},
  {"variant/length", &test_ncm_spectral_tensor_variant_length_subprocess, "*3 coefficients, the numbers of coefficients give 4*"},
  {"variant/no_domain", &test_ncm_spectral_tensor_variant_no_domain_subprocess, "*0 lower ends and 0 upper ends for 2 variables*"},
};

static const struct
{
  const gchar *name;

  void (*func) (TestNcmSpectralTensor *test, gconstpointer pdata);
} tests_per_n_vars[] = {
  {"basic", &test_ncm_spectral_tensor_basic},
  {"tree", &test_ncm_spectral_tensor_tree},
  {"eval", &test_ncm_spectral_tensor_eval},
  {"eval_long", &test_ncm_spectral_tensor_eval_long},
  {"eval_rows", &test_ncm_spectral_tensor_eval_rows},
  {"eval_array", &test_ncm_spectral_tensor_eval_array},
  {"edges_rule", &test_ncm_spectral_tensor_edges_rule},
  {"serialize", &test_ncm_spectral_tensor_serialize},
  {"reset", &test_ncm_spectral_tensor_reset},
  {"set_domain", &test_ncm_spectral_tensor_set_domain},
  {"incomplete", &test_ncm_spectral_tensor_incomplete},
  {"split_empty", &test_ncm_spectral_tensor_split_empty},
};

gint
main (gint argc, gchar *argv[])
{
  guint n, t;

  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  for (n = 1; n <= TEST_MAX_N_VARS; n++)
  {
    for (t = 0; t < G_N_ELEMENTS (tests_per_n_vars); t++)
    {
      gchar *path = g_strdup_printf ("/ncm/spectral_tensor/n%u/%s", n, tests_per_n_vars[t].name);

      g_test_add (path, TestNcmSpectralTensor, GUINT_TO_POINTER (n),
                  &test_ncm_spectral_tensor_new, tests_per_n_vars[t].func, &test_ncm_spectral_tensor_free);
      g_free (path);
    }
  }

  g_test_add ("/ncm/spectral_tensor/traps", TestNcmSpectralTensor, GUINT_TO_POINTER (2),
              &test_ncm_spectral_tensor_new, &test_ncm_spectral_tensor_traps, &test_ncm_spectral_tensor_free);

  for (t = 0; t < G_N_ELEMENTS (traps); t++)
  {
    gchar *path = g_strdup_printf ("/ncm/spectral_tensor/traps/%s/subprocess", traps[t].name);

    g_test_add (path, TestNcmSpectralTensor, GUINT_TO_POINTER (2),
                &test_ncm_spectral_tensor_new, traps[t].func, &test_ncm_spectral_tensor_free);
    g_free (path);
  }

  g_test_run ();

  return 0;
}

/* The index of the empty node of st whose subdomain starts at the lower ends lower */
static guint
_empty_node_at (TestNcmSpectralTensor *test, NcmSpectralTensor *st, const gdouble *lower)
{
  guint node;

  for (node = 0; node < ncm_spectral_tensor_get_n_nodes (st); node++)
  {
    const gdouble *lo, *up;

    ncm_spectral_tensor_peek_node_bounds (st, node, &lo, &up);

    if ((ncm_spectral_tensor_get_node_kind (st, node) == NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY) &&
        (memcmp (lo, lower, sizeof (gdouble) * test->n_vars) == 0))
      return node;
  }

  g_assert_not_reached ();

  return 0;
}

/* Sets the domain of st and divides it as described in TestNcmSpectralTensor, leaving every leaf empty */
static void
_make_tree (TestNcmSpectralTensor *test, NcmSpectralTensor *st)
{
  const guint il = test->n_vars - 1;
  const gdouble *lo, *up;
  gdouble mid;
  guint node, i;

  ncm_spectral_tensor_set_domain (st, test->n_vars, test->a, test->n_vars, test->b);

  /* The ends of the domain are passed too: they are not inside any node */
  for (i = 0; i < test->n_vars; i++)
  {
    const gdouble e[3] = {test->a[i], test->m[i], test->b[i]};

    g_assert_cmpuint (ncm_spectral_tensor_split_empty (st, i, 3, e), ==, 1u << i);
  }

  node = _empty_node_at (test, st, test->a);
  ncm_spectral_tensor_peek_node_bounds (st, node, &lo, &up);
  mid  = 0.5 * (lo[il] + up[il]);
  node = ncm_spectral_tensor_split (st, node, il, 1, &mid);

  ncm_spectral_tensor_peek_node_bounds (st, node, &lo, &up);
  mid = 0.5 * (lo[0] + up[0]);
  ncm_spectral_tensor_split (st, node, 0, 1, &mid);
}

/* Makes the node of each recorded panel a panel node of st, tree from _make_tree() */
static void
_refill (TestNcmSpectralTensor *test, NcmSpectralTensor *st)
{
  const guint n_vars = test->n_vars;
  guint p, i;

  for (p = 0; p < test->offset->len; p++)
  {
    const guint *K = &g_array_index (test->K, guint, n_vars * p);
    gsize len      = 1;

    for (i = 0; i < n_vars; i++)
      len *= K[i];

    g_assert_cmpuint (ncm_spectral_tensor_set_panel (st, g_array_index (test->node, guint, p), n_vars, K, len,
                                                     &g_array_index (test->c, gdouble, g_array_index (test->offset, gsize, p))), ==, p);
  }
}

/*
 * The panel holding x by the rule of the class, found by testing every recorded panel: a
 * coordinate is clamped to [a_i, b_i], belongs to [lower, upper), and the upper edge b_i
 * belongs to the panels that end there.
 */
static guint
_ref_panel (TestNcmSpectralTensor *test, const gdouble *x)
{
  const guint n_vars = test->n_vars;
  guint p, i;

  for (p = 0; p < test->offset->len; p++)
  {
    const gdouble *lower = &g_array_index (test->lower, gdouble, n_vars * p);
    const gdouble *upper = &g_array_index (test->upper, gdouble, n_vars * p);
    gboolean inside      = TRUE;

    for (i = 0; inside && (i < n_vars); i++)
    {
      const gdouble xc = GSL_MIN (GSL_MAX (x[i], test->a[i]), test->b[i]);

      inside = (xc >= lower[i]) && ((xc < upper[i]) || ((xc == test->b[i]) && (upper[i] == test->b[i])));
    }

    if (inside)
      return p;
  }

  g_assert_not_reached ();

  return 0;
}

/*
 * The series of recorded panel p at x by the direct sum over every coefficient, with T_k
 * from the three-term recurrence; *l1 receives the l1 norm of its coefficients.
 */
static gdouble
_ref_eval (TestNcmSpectralTensor *test, guint p, const gdouble *x, gdouble *l1)
{
  const guint n_vars   = test->n_vars;
  const gdouble *lower = &g_array_index (test->lower, gdouble, n_vars * p);
  const gdouble *upper = &g_array_index (test->upper, gdouble, n_vars * p);
  const guint *K       = &g_array_index (test->K, guint, n_vars * p);
  const gdouble *c     = &g_array_index (test->c, gdouble, g_array_index (test->offset, gsize, p));
  gdouble T[TEST_MAX_N_VARS][TEST_K_MAX];
  gdouble sum = 0.0;
  gsize len   = 1;
  gsize q;
  guint i, k;

  *l1 = 0.0;

  for (i = 0; i < n_vars; i++)
  {
    const gdouble s = (2.0 * x[i] - lower[i] - upper[i]) / (upper[i] - lower[i]);

    T[i][0] = 1.0;

    if (K[i] > 1)
      T[i][1] = s;

    /* Inside the panel T_k (s) = cos (k acos (s)), which does not share the recurrence of the class */
    for (k = 2; k < K[i]; k++)
      T[i][k] = (fabs (s) <= 1.0) ? cos (k * acos (s)) : 2.0 * s * T[i][k - 1] - T[i][k - 2];

    len *= K[i];
  }

  for (q = 0; q < len; q++)
  {
    gdouble prod = c[q];
    gsize r      = q;

    for (i = 0; i < n_vars; i++)
    {
      prod *= T[i][r % K[i]];
      r    /= K[i];
    }

    sum += prod;
    *l1 += fabs (c[q]);
  }

  return sum;
}

/* A random point of the domain of the table */
static void
_random_point (TestNcmSpectralTensor *test, gdouble *x)
{
  guint i;

  for (i = 0; i < test->n_vars; i++)
    x[i] = ncm_rng_uniform_gen (test->rng, test->a[i], test->b[i]);
}

/* Asserts that st2 has the nodes and the panels of st, and the same values at random points */
static void
_assert_same_tables (TestNcmSpectralTensor *test, NcmSpectralTensor *st, NcmSpectralTensor *st2)
{
  const guint n_vars = test->n_vars;
  gdouble x[TEST_MAX_N_VARS];
  guint node, q, i, p;

  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (st2), ==, n_vars);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st2), ==, ncm_spectral_tensor_get_n_nodes (st));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st2), ==, ncm_spectral_tensor_get_n_empty (st));
  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st2), ==, ncm_spectral_tensor_get_n_panels (st));

  for (node = 0; node < ncm_spectral_tensor_get_n_nodes (st); node++)
  {
    const gdouble *lo, *up, *lo2, *up2;

    ncm_spectral_tensor_peek_node_bounds (st, node, &lo, &up);
    ncm_spectral_tensor_peek_node_bounds (st2, node, &lo2, &up2);

    g_assert_cmpint (ncm_spectral_tensor_get_node_kind (st2, node), ==, ncm_spectral_tensor_get_node_kind (st, node));
    g_assert_cmpmem (lo2, sizeof (gdouble) * n_vars, lo, sizeof (gdouble) * n_vars);
    g_assert_cmpmem (up2, sizeof (gdouble) * n_vars, up, sizeof (gdouble) * n_vars);
  }

  for (p = 0; p < ncm_spectral_tensor_get_n_panels (st); p++)
  {
    const gdouble *lower, *upper, *lower2, *upper2;
    const guint *K, *K2;
    const gdouble *c  = ncm_spectral_tensor_peek_panel (st, p, &lower, &upper, &K);
    const gdouble *c2 = ncm_spectral_tensor_peek_panel (st2, p, &lower2, &upper2, &K2);
    gsize len         = 1;

    for (i = 0; i < n_vars; i++)
      len *= K[i];

    g_assert_cmpmem (lower2, sizeof (gdouble) * n_vars, lower, sizeof (gdouble) * n_vars);
    g_assert_cmpmem (upper2, sizeof (gdouble) * n_vars, upper, sizeof (gdouble) * n_vars);
    g_assert_cmpmem (K2, sizeof (guint) * n_vars, K, sizeof (guint) * n_vars);
    g_assert_cmpmem (c2, sizeof (gdouble) * len, c, sizeof (gdouble) * len);
  }

  if (ncm_spectral_tensor_get_n_empty (st) > 0)
    return;

  for (q = 0; q < TEST_N_POINTS; q++)
  {
    _random_point (test, x);
    g_assert_cmpfloat (ncm_spectral_tensor_eval (st2, n_vars, x), ==, ncm_spectral_tensor_eval (st, n_vars, x));
  }
}

/* A table of two variables with the domain of the fixture and the root empty */
static NcmSpectralTensor *
_domain_table (TestNcmSpectralTensor *test)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (2);

  ncm_spectral_tensor_set_domain (st, 2, test->a, 2, test->b);

  return st;
}

void
test_ncm_spectral_tensor_new (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars = GPOINTER_TO_UINT (pdata);
  guint node, i;

  g_assert_cmpuint (n_vars, >=, 1);
  g_assert_cmpuint (n_vars, <=, TEST_MAX_N_VARS);

  test->n_vars = n_vars;
  test->st     = ncm_spectral_tensor_new (n_vars);
  test->rng    = ncm_rng_seeded_new (NULL, 2026 + n_vars);
  test->node   = g_array_new (FALSE, FALSE, sizeof (guint));
  test->lower  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->upper  = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->K      = g_array_new (FALSE, FALSE, sizeof (guint));
  test->offset = g_array_new (FALSE, FALSE, sizeof (gsize));
  test->c      = g_array_new (FALSE, FALSE, sizeof (gdouble));

  for (i = 0; i < n_vars; i++)
  {
    test->a[i] = -1.0 + 0.1 * i;
    test->m[i] = 0.25 + 0.05 * i;
    test->b[i] = 1.5 + 0.2 * i;
  }

  _make_tree (test, test->st);

  /* Every empty node becomes a panel node, in the order of the nodes */
  for (node = 0; node < ncm_spectral_tensor_get_n_nodes (test->st); node++)
  {
    const gsize offset = test->c->len;
    const gdouble *lo, *up;
    guint K[TEST_MAX_N_VARS];
    gsize len = 1;
    gsize q;

    if (ncm_spectral_tensor_get_node_kind (test->st, node) != NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY)
      continue;

    for (i = 0; i < n_vars; i++)
    {
      K[i] = 1 + ncm_rng_uniform_int_gen (test->rng, TEST_K_MAX);
      len *= K[i];
    }

    for (q = 0; q < len; q++)
    {
      const gdouble v = ncm_rng_uniform_gen (test->rng, -1.0, 1.0);

      g_array_append_val (test->c, v);
    }

    ncm_spectral_tensor_peek_node_bounds (test->st, node, &lo, &up);

    g_array_append_val (test->node, node);
    g_array_append_vals (test->lower, lo, n_vars);
    g_array_append_vals (test->upper, up, n_vars);
    g_array_append_vals (test->K, K, n_vars);
    g_array_append_val (test->offset, offset);

    g_assert_cmpuint (ncm_spectral_tensor_set_panel (test->st, node, n_vars, K, len, &g_array_index (test->c, gdouble, offset)), ==, test->offset->len - 1);
  }
}

void
test_ncm_spectral_tensor_free (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_spectral_tensor_free, test->st);
  ncm_rng_free (test->rng);
  g_array_unref (test->node);
  g_array_unref (test->lower);
  g_array_unref (test->upper);
  g_array_unref (test->K);
  g_array_unref (test->offset);
  g_array_unref (test->c);
}

void
test_ncm_spectral_tensor_basic (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars = test->n_vars;
  guint i, p;

  g_assert_cmpuint (ncm_spectral_tensor_get_n_vars (test->st), ==, n_vars);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (test->st), ==, (1u << n_vars) + 2);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (test->st), ==, (2u << n_vars) + 3);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (test->st), ==, 0);

  for (i = 0; i < n_vars; i++)
  {
    gdouble a, b;

    ncm_spectral_tensor_get_range (test->st, i, &a, &b);
    g_assert_cmpfloat (a, ==, test->a[i]);
    g_assert_cmpfloat (b, ==, test->b[i]);
  }

  for (p = 0; p < test->offset->len; p++)
  {
    const guint node = g_array_index (test->node, guint, p);
    const gdouble *lower, *upper;
    const guint *K;
    const gdouble *c = ncm_spectral_tensor_peek_panel (test->st, p, &lower, &upper, &K);
    gsize len        = 1;

    g_assert_cmpint (ncm_spectral_tensor_get_node_kind (test->st, node), ==, NCM_SPECTRAL_TENSOR_NODE_KIND_PANEL);
    g_assert_cmpuint (ncm_spectral_tensor_get_node_panel (test->st, node), ==, p);

    g_assert_cmpmem (lower, sizeof (gdouble) * n_vars, &g_array_index (test->lower, gdouble, n_vars * p), sizeof (gdouble) * n_vars);
    g_assert_cmpmem (upper, sizeof (gdouble) * n_vars, &g_array_index (test->upper, gdouble, n_vars * p), sizeof (gdouble) * n_vars);
    g_assert_cmpmem (K, sizeof (guint) * n_vars, &g_array_index (test->K, guint, n_vars * p), sizeof (guint) * n_vars);

    for (i = 0; i < n_vars; i++)
      len *= K[i];

    g_assert_cmpmem (c, sizeof (gdouble) * len, &g_array_index (test->c, gdouble, g_array_index (test->offset, gsize, p)), sizeof (gdouble) * len);

    /* The copies of ncm_spectral_tensor_get_panel() */
    {
      GArray *lower_a, *upper_a, *K_a;
      GArray *c_a = ncm_spectral_tensor_get_panel (test->st, p, &lower_a, &upper_a, &K_a);

      g_assert_cmpmem (lower_a->data, sizeof (gdouble) * lower_a->len, lower, sizeof (gdouble) * n_vars);
      g_assert_cmpmem (upper_a->data, sizeof (gdouble) * upper_a->len, upper, sizeof (gdouble) * n_vars);
      g_assert_cmpmem (K_a->data, sizeof (guint) * K_a->len, K, sizeof (guint) * n_vars);
      g_assert_cmpmem (c_a->data, sizeof (gdouble) * c_a->len, c, sizeof (gdouble) * len);

      g_array_unref (lower_a);
      g_array_unref (upper_a);
      g_array_unref (K_a);
      g_array_unref (c_a);
    }
  }
}

/*
 * The tree read back node by node: the root splits the first variable at m_0 with the
 * children 1 and 2, the children of every split node are consecutive and divide its
 * interval at its interior edges, and the other intervals of a child are those of its
 * parent.
 */
void
test_ncm_spectral_tensor_tree (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars  = test->n_vars;
  const guint n_nodes = ncm_spectral_tensor_get_n_nodes (test->st);
  guint n_children    = 1;
  guint node, i, k;

  {
    guint var, first;
    GArray *edges = ncm_spectral_tensor_get_node_split (test->st, 0, &var, &first);

    g_assert_cmpuint (var, ==, 0);
    g_assert_cmpuint (first, ==, 1);
    g_assert_cmpuint (edges->len, ==, 1);
    g_assert_cmpfloat (g_array_index (edges, gdouble, 0), ==, test->m[0]);
    g_array_unref (edges);
  }

  for (node = 0; node < n_nodes; node++)
  {
    GArray *lower, *upper, *edges;
    guint var, first;

    if (ncm_spectral_tensor_get_node_kind (test->st, node) != NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT)
      continue;

    ncm_spectral_tensor_get_node_bounds (test->st, node, &lower, &upper);
    edges = ncm_spectral_tensor_get_node_split (test->st, node, &var, &first);

    g_assert_cmpuint (first + edges->len, <, n_nodes);

    for (k = 0; k <= edges->len; k++)
    {
      const gdouble *lo, *up;

      ncm_spectral_tensor_peek_node_bounds (test->st, first + k, &lo, &up);

      for (i = 0; i < n_vars; i++)
      {
        const gdouble e_lo = ((i != var) || (k == 0)) ? g_array_index (lower, gdouble, i) : g_array_index (edges, gdouble, k - 1);
        const gdouble e_up = ((i != var) || (k == edges->len)) ? g_array_index (upper, gdouble, i) : g_array_index (edges, gdouble, k);

        g_assert_cmpfloat (lo[i], ==, e_lo);
        g_assert_cmpfloat (up[i], ==, e_up);
      }
    }

    n_children += edges->len + 1;

    g_array_unref (lower);
    g_array_unref (upper);
    g_array_unref (edges);
  }

  /* Every node but the root is the child of one split node */
  g_assert_cmpuint (n_children, ==, n_nodes);
}

void
test_ncm_spectral_tensor_eval (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  gdouble x[TEST_MAX_N_VARS];
  guint q;

  for (q = 0; q < TEST_N_POINTS; q++)
  {
    guint p;
    gdouble l1, ref, val;

    _random_point (test, x);

    p   = _ref_panel (test, x);
    ref = _ref_eval (test, p, x, &l1);
    val = ncm_spectral_tensor_eval (test->st, test->n_vars, x);

    g_assert_cmpuint (ncm_spectral_tensor_find_panel (test->st, test->n_vars, x), ==, p);
    g_assert_cmpfloat (fabs (val - ref), <=, 1.0e-14 * l1);
  }
}

/*
 * A panel with more coefficients in a variable than the sums with the values T_k (s_i) take
 * (1025): the series is summed by the Clenshaw recurrences, compared with the direct sum.
 * The coefficients decay as 1 / (1 + k_0)^2.
 */
void
test_ncm_spectral_tensor_eval_long (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars    = test->n_vars;
  const guint K0        = 1026;
  const guint K1        = (n_vars > 1) ? 2 : 1;
  NcmSpectralTensor *st = ncm_spectral_tensor_new (n_vars);
  GArray *c             = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble lower[TEST_MAX_N_VARS], upper[TEST_MAX_N_VARS], x[TEST_MAX_N_VARS];
  guint K[TEST_MAX_N_VARS];
  gdouble l1 = 0.0;
  guint i, k0, k1, q;

  for (i = 0; i < n_vars; i++)
  {
    lower[i] = -1.0;
    upper[i] = 1.0;
    K[i]     = (i == 0) ? K0 : ((i == 1) ? K1 : 1);
  }

  for (k1 = 0; k1 < K1; k1++)
  {
    for (k0 = 0; k0 < K0; k0++)
    {
      const gdouble ck = ncm_rng_uniform_gen (test->rng, -1.0, 1.0) / gsl_pow_2 (1.0 + k0);

      g_array_append_val (c, ck);
      l1 += fabs (ck);
    }
  }

  ncm_spectral_tensor_set_domain (st, n_vars, lower, n_vars, upper);
  ncm_spectral_tensor_set_panel (st, 0, n_vars, K, c->len, &g_array_index (c, gdouble, 0));

  for (q = 0; q < 50; q++)
  {
    gdouble ref = 0.0;

    for (i = 0; i < n_vars; i++)
      x[i] = ncm_rng_uniform_gen (test->rng, -1.0, 1.0);

    for (k1 = 0; k1 < K1; k1++)
    {
      const gdouble T1 = (k1 == 0) ? 1.0 : x[1];

      for (k0 = 0; k0 < K0; k0++)
        ref += g_array_index (c, gdouble, k1 * K0 + k0) * cos (k0 * acos (x[0])) * T1;
    }

    g_assert_cmpfloat (fabs (ncm_spectral_tensor_eval (st, n_vars, x) - ref), <=, 1.0e-13 * l1);
  }

  g_array_unref (c);
  ncm_spectral_tensor_free (st);
}

/*
 * Every row length from 1 to 40: the sums of two to four variables have one copy of their
 * code for each row length up to 33 and one for the longer ones. The second variable has 5
 * or 40 coefficients, below and above the 33 values that select how its T_k are computed.
 * One panel on [-1, 1]^n with random coefficients, compared with the direct sum with
 * T_k (s) = cos (k acos (s)).
 */
void
test_ncm_spectral_tensor_eval_rows (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars   = test->n_vars;
  const guint K_rest[] = {5, 3, 2, 2};
  GArray *c            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble lower[TEST_MAX_N_VARS], upper[TEST_MAX_N_VARS], x[TEST_MAX_N_VARS];
  gdouble T[TEST_MAX_N_VARS][40];
  guint K[TEST_MAX_N_VARS];
  guint K0, v, i, k, q;

  for (i = 0; i < n_vars; i++)
  {
    lower[i] = -1.0;
    upper[i] = 1.0;
  }

  for (K0 = 1; K0 <= 40; K0++)
  {
    for (v = 0; v < 2; v++)
    {
      NcmSpectralTensor *st = ncm_spectral_tensor_new (n_vars);
      gdouble l1            = 0.0;
      gsize len             = 1;
      gsize r;

      for (i = 0; i < n_vars; i++)
      {
        K[i] = (i == 0) ? K0 : (((i == 1) && (v == 1)) ? 40 : K_rest[i - 1]);
        len *= K[i];
      }

      g_array_set_size (c, len);

      for (r = 0; r < len; r++)
      {
        g_array_index (c, gdouble, r) = ncm_rng_uniform_gen (test->rng, -1.0, 1.0);
        l1                           += fabs (g_array_index (c, gdouble, r));
      }

      ncm_spectral_tensor_set_domain (st, n_vars, lower, n_vars, upper);
      ncm_spectral_tensor_set_panel (st, 0, n_vars, K, len, &g_array_index (c, gdouble, 0));

      for (q = 0; q < 8; q++)
      {
        gdouble ref = 0.0;

        for (i = 0; i < n_vars; i++)
        {
          x[i] = ncm_rng_uniform_gen (test->rng, -1.0, 1.0);

          for (k = 0; k < K[i]; k++)
            T[i][k] = cos (k * acos (x[i]));
        }

        for (r = 0; r < len; r++)
        {
          gdouble prod = g_array_index (c, gdouble, r);
          gsize rr     = r;

          for (i = 0; i < n_vars; i++)
          {
            prod *= T[i][rr % K[i]];
            rr   /= K[i];
          }

          ref += prod;
        }

        g_assert_cmpfloat (fabs (ncm_spectral_tensor_eval (st, n_vars, x) - ref), <=, 1.0e-14 * l1);
      }

      ncm_spectral_tensor_free (st);
    }
  }

  g_array_unref (c);
}

/*
 * ncm_spectral_tensor_eval_array() against the direct sum, for points in random order (runs
 * of one point on a panel) and for the same points grouped by panel (runs of every length,
 * which take the passes of 8 points, of 4 points and the points left), with a new array and
 * with one reused. Then 20 points on a panel of 1026 coefficients, which the sums of
 * several points do not take.
 */
void
test_ncm_spectral_tensor_eval_array (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars   = test->n_vars;
  const guint n_points = 500;
  const guint n_panels = ncm_spectral_tensor_get_n_panels (test->st);
  GArray *x            = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *xs           = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *fval         = NULL;
  GArray *fval_s       = g_array_new (FALSE, FALSE, sizeof (gdouble));
  gdouble pt[TEST_MAX_N_VARS];
  guint p, q;

  for (q = 0; q < n_points; q++)
  {
    _random_point (test, pt);
    g_array_append_vals (x, pt, n_vars);
  }

  for (p = 0; p < n_panels; p++)
  {
    for (q = 0; q < n_points; q++)
    {
      const gdouble *xq = &g_array_index (x, gdouble, n_vars * q);

      if (_ref_panel (test, xq) == p)
        g_array_append_vals (xs, xq, n_vars);
    }
  }

  g_assert_cmpuint (xs->len, ==, x->len);

  ncm_spectral_tensor_eval_array (test->st, x, &fval);
  g_assert_nonnull (fval);
  g_assert_cmpuint (fval->len, ==, n_points);

  g_array_set_size (fval_s, 3);
  ncm_spectral_tensor_eval_array (test->st, xs, &fval_s);
  g_assert_cmpuint (fval_s->len, ==, n_points);

  for (q = 0; q < n_points; q++)
  {
    gdouble l1, ref;

    ref = _ref_eval (test, _ref_panel (test, &g_array_index (x, gdouble, n_vars * q)), &g_array_index (x, gdouble, n_vars * q), &l1);
    g_assert_cmpfloat (fabs (g_array_index (fval, gdouble, q) - ref), <=, 1.0e-14 * l1);

    ref = _ref_eval (test, _ref_panel (test, &g_array_index (xs, gdouble, n_vars * q)), &g_array_index (xs, gdouble, n_vars * q), &l1);
    g_assert_cmpfloat (fabs (g_array_index (fval_s, gdouble, q) - ref), <=, 1.0e-14 * l1);
  }

  /* No point: an empty array */
  g_array_set_size (x, 0);
  ncm_spectral_tensor_eval_array (test->st, x, &fval);
  g_assert_cmpuint (fval->len, ==, 0);

  {
    const guint K0        = 1026;
    NcmSpectralTensor *st = ncm_spectral_tensor_new (n_vars);
    GArray *c             = g_array_new (FALSE, FALSE, sizeof (gdouble));
    gdouble lower[TEST_MAX_N_VARS], upper[TEST_MAX_N_VARS];
    guint K[TEST_MAX_N_VARS];
    guint i, k0;

    for (i = 0; i < n_vars; i++)
    {
      lower[i] = -1.0;
      upper[i] = 1.0;
      K[i]     = (i == 0) ? K0 : 1;
    }

    for (k0 = 0; k0 < K0; k0++)
    {
      const gdouble ck = ncm_rng_uniform_gen (test->rng, -1.0, 1.0) / gsl_pow_2 (1.0 + k0);

      g_array_append_val (c, ck);
    }

    ncm_spectral_tensor_set_domain (st, n_vars, lower, n_vars, upper);
    ncm_spectral_tensor_set_panel (st, 0, n_vars, K, c->len, &g_array_index (c, gdouble, 0));

    for (q = 0; q < 20; q++)
    {
      for (i = 0; i < n_vars; i++)
        pt[i] = ncm_rng_uniform_gen (test->rng, -1.0, 1.0);

      g_array_append_vals (x, pt, n_vars);
    }

    ncm_spectral_tensor_eval_array (st, x, &fval);
    g_assert_cmpuint (fval->len, ==, 20);

    for (q = 0; q < 20; q++)
      g_assert_cmpfloat (g_array_index (fval, gdouble, q), ==, ncm_spectral_tensor_eval (st, n_vars, &g_array_index (x, gdouble, n_vars * q)));

    g_array_unref (c);
    ncm_spectral_tensor_free (st);
  }

  g_array_unref (x);
  g_array_unref (xs);
  g_array_unref (fval);
  g_array_unref (fval_s);
}

void
test_ncm_spectral_tensor_edges_rule (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint n_vars = test->n_vars;
  gdouble x[TEST_MAX_N_VARS];
  guint q, i;

  /*
   * Coordinates on the ends of the domain, on the interior edges of the split nodes and
   * outside the domain, against the rule: a point on an interior edge belongs to the child
   * whose interval starts there, the upper end of the domain to the panels that end there,
   * and a point outside to the panel that contains it clamped to the domain.
   */
  for (q = 0; q < TEST_N_POINTS; q++)
  {
    guint p;

    for (i = 0; i < n_vars; i++)
    {
      const gdouble half      = 0.5 * (test->a[i] + test->m[i]);
      const gdouble choices[] = {
        test->a[i], test->m[i], test->b[i], half, 0.5 * (test->a[i] + half),
        test->a[i] - 0.3, test->b[i] + 0.3,
        ncm_rng_uniform_gen (test->rng, test->a[i], test->b[i]),
      };

      x[i] = choices[ncm_rng_uniform_int_gen (test->rng, G_N_ELEMENTS (choices))];
    }

    p = _ref_panel (test, x);

    g_assert_cmpuint (ncm_spectral_tensor_find_panel (test->st, n_vars, x), ==, p);
    g_assert_cmpfloat (ncm_spectral_tensor_eval (test->st, n_vars, x), ==, ncm_spectral_tensor_eval_on (test->st, p, n_vars, x));
  }
}

void
test_ncm_spectral_tensor_serialize (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NcmSerialize *ser      = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralTensor *dup = NCM_SPECTRAL_TENSOR (ncm_serialize_dup_obj (ser, G_OBJECT (test->st)));

  _assert_same_tables (test, test->st, dup);

  ncm_spectral_tensor_free (dup);
  ncm_serialize_free (ser);
}

void
test_ncm_spectral_tensor_reset (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (test->n_vars);
  guint i;

  _make_tree (test, st);
  _refill (test, st);
  _assert_same_tables (test, test->st, st);

  /* The domain stays, and the same tree can be made again */
  ncm_spectral_tensor_reset (st);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (st), ==, 0);
  g_assert_cmpint (ncm_spectral_tensor_get_node_kind (st, 0), ==, NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY);

  for (i = 0; i < test->n_vars; i++)
  {
    gdouble a, b;

    ncm_spectral_tensor_get_range (st, i, &a, &b);
    g_assert_cmpfloat (a, ==, test->a[i]);
    g_assert_cmpfloat (b, ==, test->b[i]);
  }

  {
    const gdouble *lo, *up;
    gdouble mid;

    /* After the reset the root can be split again; its children are the nodes 1 and 2 */
    ncm_spectral_tensor_peek_node_bounds (st, 0, &lo, &up);
    mid = 0.5 * (lo[0] + up[0]);
    g_assert_cmpuint (ncm_spectral_tensor_split (st, 0, 0, 1, &mid), ==, 1);
    g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 3);
    ncm_spectral_tensor_reset (st);
  }

  _make_tree (test, st);
  _refill (test, st);
  _assert_same_tables (test, test->st, st);

  ncm_spectral_tensor_free (st);
}

/* A new domain removes the tree; a table with no domain has no node and reset leaves it so */
void
test_ncm_spectral_tensor_set_domain (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (test->n_vars);
  gdouble lower[TEST_MAX_N_VARS], upper[TEST_MAX_N_VARS];
  guint i;

  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 0);
  ncm_spectral_tensor_reset (st);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 0);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, 0);

  for (i = 0; i < test->n_vars; i++)
  {
    lower[i] = -2.0 - i;
    upper[i] = 3.0 + i;
  }

  ncm_spectral_tensor_set_domain (test->st, test->n_vars, lower, test->n_vars, upper);

  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (test->st), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (test->st), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_panels (test->st), ==, 0);

  for (i = 0; i < test->n_vars; i++)
  {
    gdouble a, b;

    ncm_spectral_tensor_get_range (test->st, i, &a, &b);
    g_assert_cmpfloat (a, ==, lower[i]);
    g_assert_cmpfloat (b, ==, upper[i]);
  }

  ncm_spectral_tensor_free (st);
}

/* A table with empty nodes is serialized with them, and its copy can be completed */
void
test_ncm_spectral_tensor_incomplete (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (test->n_vars);
  NcmSerialize *ser     = ncm_serialize_new (NCM_SERIALIZE_OPT_NONE);
  NcmSpectralTensor *dup;

  _make_tree (test, st);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, (1u << test->n_vars) + 2);

  dup = NCM_SPECTRAL_TENSOR (ncm_serialize_dup_obj (ser, G_OBJECT (st)));
  _assert_same_tables (test, st, dup);

  _refill (test, dup);
  _assert_same_tables (test, test->st, dup);

  ncm_spectral_tensor_free (dup);
  ncm_spectral_tensor_free (st);
  ncm_serialize_free (ser);
}

/* ncm_spectral_tensor_split_empty() splits the empty nodes that hold one of the edges, and no other */
void
test_ncm_spectral_tensor_split_empty (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  NcmSpectralTensor *st = ncm_spectral_tensor_new (test->n_vars);
  const gdouble inside  = 0.5 * (test->a[0] + test->m[0]);
  const gdouble out[2]  = {test->a[0] - 1.0, test->b[0] + 1.0};

  ncm_spectral_tensor_set_domain (st, test->n_vars, test->a, test->n_vars, test->b);
  ncm_spectral_tensor_split (st, 0, 0, 1, &test->m[0]);

  g_assert_cmpuint (ncm_spectral_tensor_split_empty (st, 0, 2, out), ==, 0);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 3);

  /* Only the child [a_0, m_0] holds the edge */
  g_assert_cmpuint (ncm_spectral_tensor_split_empty (st, 0, 1, &inside), ==, 1);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_nodes (st), ==, 5);
  g_assert_cmpint (ncm_spectral_tensor_get_node_kind (st, 1), ==, NCM_SPECTRAL_TENSOR_NODE_KIND_SPLIT);
  g_assert_cmpint (ncm_spectral_tensor_get_node_kind (st, 2), ==, NCM_SPECTRAL_TENSOR_NODE_KIND_EMPTY);
  g_assert_cmpuint (ncm_spectral_tensor_get_n_empty (st), ==, 3);

  ncm_spectral_tensor_free (st);
}

void
test_ncm_spectral_tensor_traps (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  guint t;

  for (t = 0; t < G_N_ELEMENTS (traps); t++)
  {
    gchar *path = g_strdup_printf ("/ncm/spectral_tensor/traps/%s/subprocess", traps[t].name);

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
test_ncm_spectral_tensor_set_domain_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_set_domain (test->st, 1, test->a, 2, test->b);
}

void
test_ncm_spectral_tensor_set_domain_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble upper[2] = {test->b[0], test->a[1]};

  ncm_spectral_tensor_set_domain (test->st, 2, test->a, 2, upper);
}

void
test_ncm_spectral_tensor_split_node_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_split (test->st, 1000, 0, 1, &test->m[0]);
}

void
test_ncm_spectral_tensor_split_not_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_split (test->st, 0, 0, 1, &test->m[0]);
}

void
test_ncm_spectral_tensor_split_variable_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_split (_domain_table (test), 0, 2, 1, &test->m[0]);
}

void
test_ncm_spectral_tensor_split_no_edges_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_split (_domain_table (test), 0, 0, 0, &test->m[0]);
}

void
test_ncm_spectral_tensor_split_outside_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  /* An edge on the lower end of the node is not inside it */
  ncm_spectral_tensor_split (_domain_table (test), 0, 0, 1, &test->a[0]);
}

void
test_ncm_spectral_tensor_split_order_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble e[2] = {test->m[0], test->m[0]};

  ncm_spectral_tensor_split (_domain_table (test), 0, 0, 2, e);
}

void
test_ncm_spectral_tensor_split_empty_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_split_empty (ncm_spectral_tensor_new (2), 0, 1, &test->m[0]);
}

void
test_ncm_spectral_tensor_split_empty_order_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble e[2] = {test->m[0], test->a[0]};

  ncm_spectral_tensor_split_empty (_domain_table (test), 0, 2, e);
}

void
test_ncm_spectral_tensor_set_panel_not_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint K[2]   = {1, 1};
  const gdouble c[1] = {1.0};

  ncm_spectral_tensor_set_panel (test->st, 0, 2, K, 1, c);
}

void
test_ncm_spectral_tensor_set_panel_n_K_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint K[3]   = {1, 1, 1};
  const gdouble c[1] = {1.0};

  ncm_spectral_tensor_set_panel (_domain_table (test), 0, 3, K, 1, c);
}

void
test_ncm_spectral_tensor_set_panel_n_coeffs_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint K[2]   = {1, 0};
  const gdouble c[1] = {1.0};

  ncm_spectral_tensor_set_panel (_domain_table (test), 0, 2, K, 1, c);
}

void
test_ncm_spectral_tensor_set_panel_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const guint K[2]   = {2, 3};
  const gdouble c[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

  ncm_spectral_tensor_set_panel (_domain_table (test), 0, 2, K, 5, c);
}

void
test_ncm_spectral_tensor_eval_array_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble v[3] = {0.0, 0.0, 0.0};
  GArray *x          = g_array_new (FALSE, FALSE, sizeof (gdouble));
  GArray *fval       = NULL;

  g_array_append_vals (x, v, 3);
  ncm_spectral_tensor_eval_array (test->st, x, &fval);
}

void
test_ncm_spectral_tensor_find_panel_n_vars_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble x[3] = {0.0, 0.0, 0.0};

  ncm_spectral_tensor_eval (test->st, 3, x);
}

void
test_ncm_spectral_tensor_find_panel_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble x[2] = {0.0, 0.0};

  ncm_spectral_tensor_eval (ncm_spectral_tensor_new (2), 2, x);
}

void
test_ncm_spectral_tensor_find_panel_empty_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  const gdouble x[2] = {0.0, 0.0};

  ncm_spectral_tensor_eval (_domain_table (test), 2, x);
}

void
test_ncm_spectral_tensor_node_not_split_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  guint var, first;

  ncm_spectral_tensor_get_node_split (test->st, g_array_index (test->node, guint, 0), &var, &first);
}

void
test_ncm_spectral_tensor_node_not_panel_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  ncm_spectral_tensor_get_node_panel (test->st, 0);
}

void
test_ncm_spectral_tensor_variant_length_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  GVariant *tree = g_variant_new_parsed ("([0.0, 0.0], [1.0, 1.0], @a(uuad) [], [(@u 0, @au [2, 2], [1.0, 2.0, 3.0])])");
  GObject *obj   = g_object_new (NCM_TYPE_SPECTRAL_TENSOR, "n-vars", 2, "tree", tree, NULL);

  g_object_unref (obj);
}

void
test_ncm_spectral_tensor_variant_no_domain_subprocess (TestNcmSpectralTensor *test, gconstpointer pdata)
{
  GVariant *tree = g_variant_new_parsed ("(@ad [], @ad [], [(@u 0, @u 0, [0.5])], @a(uauad) [])");
  GObject *obj   = g_object_new (NCM_TYPE_SPECTRAL_TENSOR, "n-vars", 2, "tree", tree, NULL);

  g_object_unref (obj);
}

/* LCOV_EXCL_STOP */

