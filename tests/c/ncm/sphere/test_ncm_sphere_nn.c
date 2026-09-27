/***************************************************************************
 *            test_ncm_sphere_nn.c
 *
 *  Sun Sep 27 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_sphere_nn.c
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

#define TEST_NPOINTS 2000

typedef struct _TestNcmSphereNN
{
  NcmSphereNN *snn;
  GArray *r;
  GArray *theta;
  GArray *phi;
} TestNcmSphereNN;

void test_ncm_sphere_nn_new (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_free (TestNcmSphereNN *test, gconstpointer pdata);

void test_ncm_sphere_nn_get (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_insert_array (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_brute_force (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_searches_agree (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_rebuild (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_repeated (TestNcmSphereNN *test, gconstpointer pdata);
void test_ncm_sphere_nn_get_near_pole (void);
void test_ncm_sphere_nn_dump (void);
void test_ncm_sphere_nn_dump_subprocess (void);

void test_ncm_sphere_nn_traps (void);
void test_ncm_sphere_nn_invalid_not_built (void);
void test_ncm_sphere_nn_invalid_stale (void);
void test_ncm_sphere_nn_invalid_k_zero (void);
void test_ncm_sphere_nn_invalid_k_large (void);
void test_ncm_sphere_nn_invalid_get (void);
void test_ncm_sphere_nn_invalid_lengths (void);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_add ("/ncm/sphere_nn/get", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_get, &test_ncm_sphere_nn_free);
  g_test_add ("/ncm/sphere_nn/insert_array", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_insert_array, &test_ncm_sphere_nn_free);
  g_test_add ("/ncm/sphere_nn/brute_force", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_brute_force, &test_ncm_sphere_nn_free);
  g_test_add ("/ncm/sphere_nn/searches_agree", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_searches_agree, &test_ncm_sphere_nn_free);
  g_test_add ("/ncm/sphere_nn/rebuild", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_rebuild, &test_ncm_sphere_nn_free);

  g_test_add ("/ncm/sphere_nn/repeated", TestNcmSphereNN, NULL,
              &test_ncm_sphere_nn_new, &test_ncm_sphere_nn_repeated, &test_ncm_sphere_nn_free);
  g_test_add_func ("/ncm/sphere_nn/get_near_pole", &test_ncm_sphere_nn_get_near_pole);
  g_test_add_func ("/ncm/sphere_nn/dump", &test_ncm_sphere_nn_dump);
  g_test_add_func ("/ncm/sphere_nn/dump/subprocess", &test_ncm_sphere_nn_dump_subprocess);
  g_test_add_func ("/ncm/sphere_nn/traps", &test_ncm_sphere_nn_traps);
  g_test_add_func ("/ncm/sphere_nn/invalid/not_built/subprocess", &test_ncm_sphere_nn_invalid_not_built);
  g_test_add_func ("/ncm/sphere_nn/invalid/stale/subprocess", &test_ncm_sphere_nn_invalid_stale);
  g_test_add_func ("/ncm/sphere_nn/invalid/k_zero/subprocess", &test_ncm_sphere_nn_invalid_k_zero);
  g_test_add_func ("/ncm/sphere_nn/invalid/k_large/subprocess", &test_ncm_sphere_nn_invalid_k_large);
  g_test_add_func ("/ncm/sphere_nn/invalid/get/subprocess", &test_ncm_sphere_nn_invalid_get);
  g_test_add_func ("/ncm/sphere_nn/invalid/lengths/subprocess", &test_ncm_sphere_nn_invalid_lengths);

  g_test_run ();

  return 0;
}

/* Points uniform on the sphere with radii in [0.5, 2], from a fixed seed. */
void
test_ncm_sphere_nn_new (TestNcmSphereNN *test, gconstpointer pdata)
{
  NcmRNG *rng = ncm_rng_seeded_new (NULL, 20260927);
  guint i;

  test->snn   = ncm_sphere_nn_new ();
  test->r     = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->theta = g_array_new (FALSE, FALSE, sizeof (gdouble));
  test->phi   = g_array_new (FALSE, FALSE, sizeof (gdouble));

  for (i = 0; i < TEST_NPOINTS; i++)
  {
    const gdouble r     = ncm_rng_uniform_gen (rng, 0.5, 2.0);
    const gdouble theta = acos (ncm_rng_uniform_gen (rng, -1.0, 1.0));
    const gdouble phi   = ncm_rng_uniform_gen (rng, -M_PI, M_PI);

    g_array_append_val (test->r, r);
    g_array_append_val (test->theta, theta);
    g_array_append_val (test->phi, phi);
  }

  ncm_sphere_nn_insert_array (test->snn, test->r, test->theta, test->phi);
  ncm_sphere_nn_rebuild (test->snn);

  ncm_rng_free (rng);
}

void
test_ncm_sphere_nn_free (TestNcmSphereNN *test, gconstpointer pdata)
{
  g_array_unref (test->r);
  g_array_unref (test->theta);
  g_array_unref (test->phi);
  NCM_TEST_FREE (ncm_sphere_nn_free, test->snn);
}

static void
_test_cartesian (const gdouble r, const gdouble theta, const gdouble phi, gdouble x[3])
{
  x[0] = r * sin (theta) * cos (phi);
  x[1] = r * sin (theta) * sin (phi);
  x[2] = r * cos (theta);
}

/* The stored point comes back as inserted, up to rounding. */
void
test_ncm_sphere_nn_get (TestNcmSphereNN *test, gconstpointer pdata)
{
  guint i;

  g_assert_cmpint (ncm_sphere_nn_get_n (test->snn), ==, TEST_NPOINTS);

  for (i = 0; i < TEST_NPOINTS; i++)
  {
    gdouble r, theta, phi;

    ncm_sphere_nn_get (test->snn, i, &r, &theta, &phi);
    ncm_assert_cmpdouble_e (r, ==, g_array_index (test->r, gdouble, i), 1.0e-14, 0.0);
    ncm_assert_cmpdouble_e (theta, ==, g_array_index (test->theta, gdouble, i), 1.0e-13, 1.0e-15);
    ncm_assert_cmpdouble_e (phi, ==, g_array_index (test->phi, gdouble, i), 1.0e-13, 1.0e-15);
  }
}

/* Inserting one by one and as an array store the same points. */
void
test_ncm_sphere_nn_insert_array (TestNcmSphereNN *test, gconstpointer pdata)
{
  NcmSphereNN *one = ncm_sphere_nn_new ();
  guint i;

  for (i = 0; i < TEST_NPOINTS; i++)
    ncm_sphere_nn_insert (one, g_array_index (test->r, gdouble, i), g_array_index (test->theta, gdouble, i), g_array_index (test->phi, gdouble, i));

  g_assert_cmpint (ncm_sphere_nn_get_n (one), ==, ncm_sphere_nn_get_n (test->snn));

  for (i = 0; i < TEST_NPOINTS; i++)
  {
    gdouble r1, t1, p1, r2, t2, p2;

    ncm_sphere_nn_get (one, i, &r1, &t1, &p1);
    ncm_sphere_nn_get (test->snn, i, &r2, &t2, &p2);
    g_assert_cmpfloat (r1, ==, r2);
    g_assert_cmpfloat (t1, ==, t2);
    g_assert_cmpfloat (p1, ==, p2);
  }

  ncm_sphere_nn_free (one);
}

static gint
_test_cmp_double (gconstpointer a, gconstpointer b)
{
  const gdouble x = *(const gdouble *) a;
  const gdouble y = *(const gdouble *) b;

  return (x < y) ? -1 : (x > y);
}

/* The k nearest by squared Euclidean distance, against sorting all distances, for
 * targets off the stored points and for each k. */
void
test_ncm_sphere_nn_brute_force (TestNcmSphereNN *test, gconstpointer pdata)
{
  NcmRNG *rng       = ncm_rng_seeded_new (NULL, 7);
  const gint64 ks[] = {1, 7, 50};
  gdouble *d2       = g_new (gdouble, TEST_NPOINTS);
  guint t, j, i;

  for (t = 0; t < 40; t++)
  {
    const gdouble r     = ncm_rng_uniform_gen (rng, 0.5, 2.0);
    const gdouble theta = acos (ncm_rng_uniform_gen (rng, -1.0, 1.0));
    const gdouble phi   = ncm_rng_uniform_gen (rng, -M_PI, M_PI);
    gdouble x[3];

    _test_cartesian (r, theta, phi, x);

    for (i = 0; i < TEST_NPOINTS; i++)
    {
      gdouble y[3];

      _test_cartesian (g_array_index (test->r, gdouble, i), g_array_index (test->theta, gdouble, i), g_array_index (test->phi, gdouble, i), y);
      d2[i] = gsl_pow_2 (x[0] - y[0]) + gsl_pow_2 (x[1] - y[1]) + gsl_pow_2 (x[2] - y[2]);
    }

    for (j = 0; j < G_N_ELEMENTS (ks); j++)
    {
      gdouble *sorted = g_memdup2 (d2, sizeof (gdouble) * TEST_NPOINTS);
      GArray *dist, *idx;

      qsort (sorted, TEST_NPOINTS, sizeof (gdouble), _test_cmp_double);
      ncm_sphere_nn_knn_search_distances (test->snn, r, theta, phi, ks[j], &dist, &idx);

      g_assert_cmpuint (dist->len, ==, ks[j]);
      g_assert_cmpuint (idx->len, ==, ks[j]);

      for (i = 0; i < ks[j]; i++)
      {
        const glong index = g_array_index (idx, glong, i);

        /* Nearest first, the same squared distances as the sorted list, and each index
         * at its own distance. */
        ncm_assert_cmpdouble_e (g_array_index (dist, gdouble, i), ==, sorted[i], 1.0e-12, 0.0);
        ncm_assert_cmpdouble_e (g_array_index (dist, gdouble, i), ==, d2[index], 1.0e-12, 0.0);

        if (i > 0)
          g_assert_cmpfloat (g_array_index (dist, gdouble, i), >=, g_array_index (dist, gdouble, i - 1));
      }

      g_array_unref (dist);
      g_array_unref (idx);
      g_free (sorted);
    }
  }

  g_free (d2);
  ncm_rng_free (rng);
}

/* A stored point is its own nearest neighbour, and the three searches return the same. */
void
test_ncm_sphere_nn_searches_agree (TestNcmSphereNN *test, gconstpointer pdata)
{
  const gint64 k = 5;
  GArray *bdist, *bidx;
  guint i, j;

  ncm_sphere_nn_knn_search_distances_batch (test->snn, test->r, test->theta, test->phi, k, &bdist, &bidx);
  g_assert_cmpuint (bdist->len, ==, k * TEST_NPOINTS);

  for (i = 0; i < TEST_NPOINTS; i += 37)
  {
    const gdouble r     = g_array_index (test->r, gdouble, i);
    const gdouble theta = g_array_index (test->theta, gdouble, i);
    const gdouble phi   = g_array_index (test->phi, gdouble, i);
    GArray *idx         = ncm_sphere_nn_knn_search (test->snn, r, theta, phi, k);
    GArray *dist, *idx2;

    ncm_sphere_nn_knn_search_distances (test->snn, r, theta, phi, k, &dist, &idx2);

    g_assert_cmpint (g_array_index (idx, glong, 0), ==, i);
    g_assert_cmpfloat (g_array_index (dist, gdouble, 0), ==, 0.0);

    for (j = 0; j < k; j++)
    {
      g_assert_cmpint (g_array_index (idx, glong, j), ==, g_array_index (idx2, glong, j));
      g_assert_cmpint (g_array_index (bidx, glong, i * k + j), ==, g_array_index (idx2, glong, j));
      g_assert_cmpfloat (g_array_index (bdist, gdouble, i * k + j), ==, g_array_index (dist, gdouble, j));
    }

    g_array_unref (idx);
    g_array_unref (idx2);
    g_array_unref (dist);
  }

  g_array_unref (bdist);
  g_array_unref (bidx);
}

/* Points inserted after a rebuild enter the searches at the next one. */
void
test_ncm_sphere_nn_rebuild (TestNcmSphereNN *test, gconstpointer pdata)
{
  const gdouble r = 1.0, theta = 0.3, phi = 1.1;
  GArray *idx;

  ncm_sphere_nn_insert (test->snn, r, theta, phi);
  ncm_sphere_nn_rebuild (test->snn);

  idx = ncm_sphere_nn_knn_search (test->snn, r, theta, phi, 1);
  g_assert_cmpint (g_array_index (idx, glong, 0), ==, TEST_NPOINTS);
  g_array_unref (idx);
}

static NcmSphereNN *
_test_two_points (gboolean build)
{
  NcmSphereNN *snn = ncm_sphere_nn_new ();

  ncm_sphere_nn_insert (snn, 1.0, 0.5, 0.5);
  ncm_sphere_nn_insert (snn, 1.0, 1.5, 2.5);

  if (build)
    ncm_sphere_nn_rebuild (snn);

  return snn;
}

/* Repeated points tie at the same distance: a stored point that was inserted twice has
 * both copies, at distance zero, as its two nearest neighbours. */
void
test_ncm_sphere_nn_repeated (TestNcmSphereNN *test, gconstpointer pdata)
{
  const guint nrepeat = 500;
  guint i;

  for (i = 0; i < nrepeat; i++)
    ncm_sphere_nn_insert (test->snn, g_array_index (test->r, gdouble, i), g_array_index (test->theta, gdouble, i), g_array_index (test->phi, gdouble, i));

  ncm_sphere_nn_rebuild (test->snn);

  for (i = 0; i < nrepeat; i += 7)
  {
    GArray *dist, *idx;
    glong a, b;

    ncm_sphere_nn_knn_search_distances (test->snn, g_array_index (test->r, gdouble, i), g_array_index (test->theta, gdouble, i), g_array_index (test->phi, gdouble, i), 3, &dist, &idx);

    a = GSL_MIN (g_array_index (idx, glong, 0), g_array_index (idx, glong, 1));
    b = GSL_MAX (g_array_index (idx, glong, 0), g_array_index (idx, glong, 1));

    g_assert_cmpint (a, ==, i);
    g_assert_cmpint (b, ==, TEST_NPOINTS + i);
    g_assert_cmpfloat (g_array_index (dist, gdouble, 1), ==, 0.0);
    g_assert_cmpfloat (g_array_index (dist, gdouble, 2), >, 0.0);

    g_array_unref (dist);
    g_array_unref (idx);
  }
}

void
test_ncm_sphere_nn_dump (void)
{
  g_test_trap_subprocess ("/ncm/sphere_nn/dump/subprocess", 0, G_TEST_SUBPROCESS_INHERIT_STDERR);
  g_test_trap_assert_passed ();
  g_test_trap_assert_stdout ("*+-------*");
}

void
test_ncm_sphere_nn_dump_subprocess (void)
{
  NcmSphereNN *snn = _test_two_points (TRUE);

  ncm_sphere_nn_dump_tree (snn);
  ncm_sphere_nn_free (snn);
}

void
test_ncm_sphere_nn_traps (void)
{
  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/not_built/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*call ncm_sphere_nn_rebuild() after inserting*");

  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/stale/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*3 points inserted but the tree was built with 2*");

  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/k_zero/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*k = 0 is out of range*");

  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/k_large/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*k = 3 is out of range, the tree holds 2 points*");

  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/get/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*index 2 out of range*");

  g_test_trap_subprocess ("/ncm/sphere_nn/invalid/lengths/subprocess", 0, 0);
  g_test_trap_assert_failed ();
  g_test_trap_assert_stderr ("*different lengths*");
}

/* A search before any rebuild used to dereference an empty result. */
void
test_ncm_sphere_nn_invalid_not_built (void)
{
  NcmSphereNN *snn = _test_two_points (FALSE);

  ncm_sphere_nn_knn_search (snn, 1.0, 0.5, 0.5, 1);
}

/* Points inserted after the rebuild used to be silently missing from the searches. */
void
test_ncm_sphere_nn_invalid_stale (void)
{
  NcmSphereNN *snn = _test_two_points (TRUE);

  ncm_sphere_nn_insert (snn, 1.0, 2.5, 0.1);
  ncm_sphere_nn_knn_search (snn, 1.0, 0.5, 0.5, 1);
}

void
test_ncm_sphere_nn_invalid_k_zero (void)
{
  NcmSphereNN *snn = _test_two_points (TRUE);
  GArray *dist, *idx;

  ncm_sphere_nn_knn_search_distances (snn, 1.0, 0.5, 0.5, 0, &dist, &idx);
}

/* The single searches returned fewer than k, the batch one aborted; both now abort. */
void
test_ncm_sphere_nn_invalid_k_large (void)
{
  NcmSphereNN *snn = _test_two_points (TRUE);

  ncm_sphere_nn_knn_search (snn, 1.0, 0.5, 0.5, 3);
}

void
test_ncm_sphere_nn_invalid_get (void)
{
  NcmSphereNN *snn = _test_two_points (TRUE);
  gdouble r, theta, phi;

  ncm_sphere_nn_get (snn, 2, &r, &theta, &phi);
}

void
test_ncm_sphere_nn_invalid_lengths (void)
{
  NcmSphereNN *snn = ncm_sphere_nn_new ();
  GArray *a        = g_array_new (FALSE, TRUE, sizeof (gdouble));
  GArray *b        = g_array_new (FALSE, TRUE, sizeof (gdouble));

  g_array_set_size (a, 3);
  g_array_set_size (b, 2);

  ncm_sphere_nn_insert_array (snn, a, a, b);
}

/* The polar angle of a stored point comes back at full precision near the poles:
 * acos (z / r) lost it as 1 / theta^2 (4.1e-8 at theta = 1e-5, 4.0e-4 at 1e-7). */
void
test_ncm_sphere_nn_get_near_pole (void)
{
  NcmSphereNN *snn        = ncm_sphere_nn_new ();
  const gdouble theta_a[] = {1.0e-3, 1.0e-5, 1.0e-7, M_PI - 1.0e-7};
  guint i;

  for (i = 0; i < G_N_ELEMENTS (theta_a); i++)
    ncm_sphere_nn_insert (snn, 1.0, theta_a[i], 0.3);

  for (i = 0; i < G_N_ELEMENTS (theta_a); i++)
  {
    gdouble r, theta, phi;

    ncm_sphere_nn_get (snn, i, &r, &theta, &phi);
    ncm_assert_cmpdouble_e (theta, ==, theta_a[i], 1.0e-15, 0.0);
  }

  ncm_sphere_nn_free (snn);
}

