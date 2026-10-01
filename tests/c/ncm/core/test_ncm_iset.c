/***************************************************************************
 *            test_ncm_iset.c
 *
 *  Thu September 25 12:00:00 2026
 *  Copyright  2026  Sandro Dias Pinto Vitenti
 *  <vitenti@uel.br>
 ****************************************************************************/
/*
 * test_ncm_iset.c
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

#define TEST_ISET_N 10

/* The set {2, 4, 7}, added out of order */
static NcmISet *
_test_iset_new (void)
{
  NcmISet *iset = ncm_iset_new (TEST_ISET_N);

  ncm_iset_add (iset, 7);
  ncm_iset_add_range (iset, 2, 5);
  ncm_iset_del (iset, 3);

  return iset;
}

static const gint _test_iset_idx[] = {2, 4, 7};

static void
_test_iset_assert_equal (NcmISet *iset, const gint *idx, guint len)
{
  NcmVector *v = ncm_vector_new (TEST_ISET_N);
  NcmVector *sub;
  guint i;

  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, i);

  g_assert_cmpuint (ncm_iset_get_len (iset), ==, len);

  if (len > 0)
  {
    sub = ncm_iset_get_subvector (iset, v, NULL);

    for (i = 0; i < len; i++)
      g_assert_cmpfloat (ncm_vector_get (sub, i), ==, idx[i]);

    ncm_vector_free (sub);
  }

  ncm_vector_free (v);
}

static void
test_ncm_iset_add_del (void)
{
  NcmISet *iset = _test_iset_new ();
  guint max_index;

  g_assert_cmpuint (ncm_iset_get_max_size (iset), ==, TEST_ISET_N);
  g_object_get (iset, "max-index", &max_index, NULL);
  g_assert_cmpuint (max_index, ==, TEST_ISET_N);

  _test_iset_assert_equal (iset, _test_iset_idx, 3);

  ncm_iset_reset (iset);
  _test_iset_assert_equal (iset, NULL, 0);

  ncm_iset_free (iset);
}

static void
test_ncm_iset_copy_complement (void)
{
  NcmISet *iset          = _test_iset_new ();
  NcmISet *target        = ncm_iset_new (TEST_ISET_N);
  const gint cmplm_idx[] = {0, 1, 3, 5, 6, 8, 9};

  ncm_iset_add (target, 9);
  ncm_iset_copy (iset, target);
  _test_iset_assert_equal (target, _test_iset_idx, 3);

  ncm_iset_set_complement (iset, target);
  _test_iset_assert_equal (target, cmplm_idx, G_N_ELEMENTS (cmplm_idx));

  /* Removing the complement from the full set leaves the set */
  ncm_iset_reset (iset);
  ncm_iset_add_range (iset, 0, TEST_ISET_N);
  ncm_iset_remove_subset (target, iset);
  _test_iset_assert_equal (iset, _test_iset_idx, 3);

  ncm_iset_free (iset);
  ncm_iset_free (target);
}

static void
test_ncm_iset_subvector (void)
{
  NcmISet *iset   = _test_iset_new ();
  NcmVector *v    = ncm_vector_new (TEST_ISET_N);
  NcmVector *w    = ncm_vector_new (TEST_ISET_N);
  NcmVector *vdup = ncm_vector_new (TEST_ISET_N);
  NcmVector *sub, *sub_dup;
  guint i;

  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, i * i);

  sub     = ncm_iset_get_subvector (iset, v, NULL);
  sub_dup = ncm_iset_get_subvector (iset, v, vdup);

  g_assert_cmpuint (ncm_vector_len (sub), ==, 3);
  g_assert_cmpuint (ncm_vector_len (sub_dup), ==, 3);

  for (i = 0; i < 3; i++)
  {
    const gint k = _test_iset_idx[i];

    g_assert_cmpfloat (ncm_vector_get (sub, i), ==, k * k);
    g_assert_cmpfloat (ncm_vector_get (vdup, i), ==, k * k);
  }

  /* ncm_iset_set_subvector() writes the components back */
  ncm_vector_set_zero (w);
  ncm_iset_set_subvector (iset, w, sub);

  for (i = 0; i < TEST_ISET_N; i++)
    g_assert_cmpfloat (ncm_vector_get (w, i), ==, (i == 2 || i == 4 || i == 7) ? i * i : 0.0);

  ncm_vector_free (sub);
  ncm_vector_free (sub_dup);
  ncm_vector_free (vdup);
  ncm_vector_free (w);
  ncm_vector_free (v);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_subarray (void)
{
  NcmISet *iset = _test_iset_new ();
  GArray *a     = g_array_new (FALSE, FALSE, sizeof (gint));
  GArray *a_dup = g_array_new (FALSE, FALSE, sizeof (gint));
  GArray *sub;
  gint i;

  for (i = 0; i < TEST_ISET_N; i++)
  {
    const gint a_i = 3 * i;

    g_array_append_val (a, a_i);
  }

  g_array_set_size (a_dup, TEST_ISET_N);
  sub = ncm_iset_get_subarray (iset, a, a_dup);
  g_assert_true (sub == a_dup);

  for (i = 0; i < 3; i++)
    g_assert_cmpint (g_array_index (sub, gint, i), ==, 3 * _test_iset_idx[i]);

  g_array_unref (sub);

  sub = ncm_iset_get_subarray (iset, a, NULL);
  g_assert_cmpuint (sub->len, ==, 3);

  for (i = 0; i < 3; i++)
    g_assert_cmpint (g_array_index (sub, gint, i), ==, 3 * _test_iset_idx[i]);

  g_array_unref (sub);

  g_array_unref (a_dup);
  g_array_unref (a);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_submatrix (void)
{
  NcmISet *iset   = _test_iset_new ();
  NcmMatrix *M    = ncm_matrix_new (TEST_ISET_N, TEST_ISET_N);
  NcmMatrix *Mdup = ncm_matrix_new (TEST_ISET_N, TEST_ISET_N);
  NcmMatrix *S, *S_dup, *U, *L;
  guint i, j;

  for (i = 0; i < TEST_ISET_N; i++)
    for (j = 0; j < TEST_ISET_N; j++)
      ncm_matrix_set (M, i, j, 10.0 * i + j);

  S     = ncm_iset_get_submatrix (iset, M, NULL);
  S_dup = ncm_iset_get_submatrix (iset, M, Mdup);
  U     = ncm_iset_get_sym_submatrix (iset, 'U', M, NULL);
  L     = ncm_iset_get_sym_submatrix (iset, 'L', M, NULL);

  g_assert_cmpuint (ncm_matrix_nrows (S), ==, 3);
  g_assert_cmpuint (ncm_matrix_ncols (S), ==, 3);

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      const gdouble M_kl = ncm_matrix_get (M, _test_iset_idx[i], _test_iset_idx[j]);

      g_assert_cmpfloat (ncm_matrix_get (S, i, j), ==, M_kl);
      g_assert_cmpfloat (ncm_matrix_get (S_dup, i, j), ==, M_kl);

      if (j >= i)
        g_assert_cmpfloat (ncm_matrix_get (U, i, j), ==, M_kl);

      if (j <= i)
        g_assert_cmpfloat (ncm_matrix_get (L, i, j), ==, M_kl);
    }
  }

  ncm_matrix_free (S);
  ncm_matrix_free (S_dup);
  ncm_matrix_free (U);
  ncm_matrix_free (L);
  ncm_matrix_free (Mdup);
  ncm_matrix_free (M);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_submatrix_cols (void)
{
  NcmISet *iset   = _test_iset_new ();
  NcmMatrix *M    = ncm_matrix_new (4, TEST_ISET_N);
  NcmMatrix *Mdup = ncm_matrix_new (4, TEST_ISET_N);
  NcmMatrix *S, *S_dup, *C;
  guint i, j;

  for (i = 0; i < 4; i++)
    for (j = 0; j < TEST_ISET_N; j++)
      ncm_matrix_set (M, i, j, 10.0 * i + j);

  S     = ncm_iset_get_submatrix_cols (iset, M, NULL);
  S_dup = ncm_iset_get_submatrix_cols (iset, M, Mdup);
  C     = ncm_iset_get_submatrix_colmajor_cols (iset, M, NULL);

  g_assert_cmpuint (ncm_matrix_nrows (S), ==, 4);
  g_assert_cmpuint (ncm_matrix_ncols (S), ==, 3);

  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 3; j++)
    {
      const gdouble M_ik = ncm_matrix_get (M, i, _test_iset_idx[j]);

      g_assert_cmpfloat (ncm_matrix_get (S, i, j), ==, M_ik);
      g_assert_cmpfloat (ncm_matrix_get (S_dup, i, j), ==, M_ik);
      g_assert_cmpfloat (ncm_matrix_get_colmajor (C, i, j), ==, M_ik);
    }
  }

  ncm_matrix_free (S);
  ncm_matrix_free (S_dup);
  ncm_matrix_free (C);
  ncm_matrix_free (Mdup);
  ncm_matrix_free (M);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_vector_max (void)
{
  NcmISet *iset = _test_iset_new ();
  NcmVector *v  = ncm_vector_new (TEST_ISET_N);
  gint max_i;
  guint i;

  /* The global maximum, at 9, is outside the set */
  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, (i == 4) ? 50.0 : i);

  g_assert_cmpfloat (ncm_iset_get_vector_max (iset, v, &max_i), ==, 50.0);
  g_assert_cmpint (max_i, ==, 4);

  ncm_iset_reset (iset);
  g_assert_cmpfloat (ncm_iset_get_vector_max (iset, v, &max_i), ==, -INFINITY);
  g_assert_cmpint (max_i, ==, -1);

  ncm_vector_free (v);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_subset_vec_lt (void)
{
  NcmISet *iset       = _test_iset_new ();
  NcmISet *out        = ncm_iset_new (TEST_ISET_N);
  NcmVector *v        = ncm_vector_new (TEST_ISET_N);
  const gint lt_idx[] = {2, 4};
  guint i;

  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, i);

  ncm_iset_add (out, 0);
  ncm_iset_get_subset_vec_lt (iset, out, v, 5.0);
  _test_iset_assert_equal (out, lt_idx, G_N_ELEMENTS (lt_idx));

  ncm_vector_free (v);
  ncm_iset_free (out);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_remove_smallest_subset (void)
{
  NcmISet *iset   = _test_iset_new ();
  NcmISet *target = ncm_iset_new (TEST_ISET_N);
  NcmVector *v    = ncm_vector_new (TEST_ISET_N);
  guint i;

  /* Decreasing values: the smallest components of the set are at 7 and 4 */
  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, TEST_ISET_N - i);

  ncm_iset_add_range (target, 0, TEST_ISET_N);
  g_assert_cmpuint (ncm_iset_remove_smallest_subset (iset, target, v, 2), ==, 2);
  {
    const gint left2_idx[] = {0, 1, 2, 3, 5, 6, 8, 9};

    _test_iset_assert_equal (target, left2_idx, G_N_ELEMENTS (left2_idx));
  }

  /* At least as many as the set has: removes all of them */
  ncm_iset_reset (target);
  ncm_iset_add_range (target, 0, TEST_ISET_N);
  g_assert_cmpuint (ncm_iset_remove_smallest_subset (iset, target, v, 5), ==, 3);
  g_assert_cmpuint (ncm_iset_get_len (target), ==, TEST_ISET_N - 3);


  ncm_vector_free (v);
  ncm_iset_free (target);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_add_largest_subset (void)
{
  NcmISet *iset = _test_iset_new ();
  NcmVector *v  = ncm_vector_new (TEST_ISET_N);
  guint i;

  for (i = 0; i < TEST_ISET_N; i++)
    ncm_vector_set (v, i, i);

  /* Candidates with v_i > 0.5 outside the set: {1, 3, 5, 6, 8, 9}, k = 6 */
  g_assert_cmpuint (ncm_iset_add_largest_subset (iset, v, 0.5, 0.5), ==, 3);
  {
    const gint idx[] = {2, 4, 6, 7, 8, 9};

    _test_iset_assert_equal (iset, idx, G_N_ELEMENTS (idx));
  }

  /* Candidates {1, 3, 5}: k f < 1, so one index is added */
  g_assert_cmpuint (ncm_iset_add_largest_subset (iset, v, 0.5, 0.01), ==, 1);
  {
    const gint idx[] = {2, 4, 5, 6, 7, 8, 9};

    _test_iset_assert_equal (iset, idx, G_N_ELEMENTS (idx));
  }

  /* No candidate */
  g_assert_cmpuint (ncm_iset_add_largest_subset (iset, v, 100.0, 1.0), ==, 0);

  /* Full set */
  ncm_iset_add_range (iset, 0, 2);
  ncm_iset_add (iset, 3);
  g_assert_cmpuint (ncm_iset_add_largest_subset (iset, v, 0.5, 1.0), ==, 0);

  ncm_vector_free (v);
  ncm_iset_free (iset);
}

static void
test_ncm_iset_vector_inv_cmp (void)
{
  NcmISet *iset = _test_iset_new ();
  NcmVector *u  = ncm_vector_new (TEST_ISET_N);
  NcmVector *v  = ncm_vector_new (TEST_ISET_N);
  NcmVector *cmp;
  guint i;

  for (i = 0; i < TEST_ISET_N; i++)
  {
    ncm_vector_set (u, i, 1.0 + i);
    ncm_vector_set (v, i, 0.5 * (1.0 + i));
  }

  cmp = ncm_iset_get_vector_inv_cmp (iset, u, v, NULL);
  g_assert_cmpuint (ncm_vector_len (cmp), ==, 3);

  for (i = 0; i < 3; i++)
    g_assert_cmpfloat (ncm_vector_get (cmp, i), ==, 2.0);

  ncm_vector_free (cmp);
  ncm_vector_free (u);
  ncm_vector_free (v);
  ncm_iset_free (iset);
}

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add_func ("/ncm/iset/add_del", &test_ncm_iset_add_del);
  g_test_add_func ("/ncm/iset/copy_complement", &test_ncm_iset_copy_complement);
  g_test_add_func ("/ncm/iset/subvector", &test_ncm_iset_subvector);
  g_test_add_func ("/ncm/iset/subarray", &test_ncm_iset_subarray);
  g_test_add_func ("/ncm/iset/submatrix", &test_ncm_iset_submatrix);
  g_test_add_func ("/ncm/iset/submatrix_cols", &test_ncm_iset_submatrix_cols);
  g_test_add_func ("/ncm/iset/vector_max", &test_ncm_iset_vector_max);
  g_test_add_func ("/ncm/iset/subset_vec_lt", &test_ncm_iset_subset_vec_lt);
  g_test_add_func ("/ncm/iset/remove_smallest_subset", &test_ncm_iset_remove_smallest_subset);
  g_test_add_func ("/ncm/iset/add_largest_subset", &test_ncm_iset_add_largest_subset);
  g_test_add_func ("/ncm/iset/vector_inv_cmp", &test_ncm_iset_vector_inv_cmp);

  g_test_run ();
}

