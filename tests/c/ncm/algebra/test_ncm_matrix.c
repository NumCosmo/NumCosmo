/***************************************************************************
 *            test_ncm_matrix.c
 *
 *  Sat April 21 14:30:26 2012
 *  Copyright  2012  Mariana Penna Lima
 *  <pennalima@gmail.com>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Mariana Penna Lima 2012 <pennalima@gmail.com>
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

typedef struct _TestNcmMatrix
{
  NcmMatrix *m;
  guint nrows;
  guint ncols;
  gdouble *d;
} TestNcmMatrix;

void test_ncm_matrix_new (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_gsl (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_array (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_data_slice (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_data_malloc (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_data_static (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_new_data_static_tda (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_free (TestNcmMatrix *test, gconstpointer pdata);

void test_ncm_matrix_sanity (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_operations (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_colmajor (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_add_mul (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_log_exp (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_square_to_sym (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_triang_to_sym (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_zero_triangle (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_update_vector (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_sym_update_vector (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_dtrmm_dtrsm (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_dtrmv_dtrsv (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_dsyrk (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_scale_rows_cols (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_sub_row_vector (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_chol_chi2_cols (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_is_identity (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_cholesky_decomp_nearPD (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_free (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_submatrix (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_serialization (TestNcmMatrix *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  g_test_set_nonfatal_assertions ();

  g_test_add ("/ncm/matrix/new/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_gsl/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_gsl,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_array/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_array,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_data_slice/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_data_slice,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_data_malloc/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_data_malloc,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_data_static/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_data_static,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/new_data_static_tda/sanity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new_data_static_tda,
              &test_ncm_matrix_sanity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/operations", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_operations,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/colmajor", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_colmajor,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/add_mul", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_add_mul,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/log_exp", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_log_exp,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/square_to_sym", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_square_to_sym,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/triang_to_sym", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_triang_to_sym,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/zero_triangle", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_zero_triangle,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/update_vector", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_update_vector,
              &test_ncm_matrix_free);
  g_test_add ("/ncm/matrix/sym_update_vector", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_sym_update_vector,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/dtrmm_dtrsm", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_dtrmm_dtrsm,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/dtrmv_dtrsv", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_dtrmv_dtrsv,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/dsyrk", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_dsyrk,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/scale_rows_cols", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_scale_rows_cols,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/sub_row_vector", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_sub_row_vector,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/chol_chi2_cols", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_chol_chi2_cols,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/is_identity", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_is_identity,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/cholesky_decomp_nearPD", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_cholesky_decomp_nearPD,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/submatrix", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_submatrix,
              &test_ncm_matrix_free);

  g_test_add ("/ncm/matrix/serialization", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_serialization,
              &test_ncm_matrix_free);

  g_test_run ();
}

void
test_ncm_matrix_new (TestNcmMatrix *test, gconstpointer pdata)
{
  guint i, j;

  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);

  test->m = ncm_matrix_new (test->nrows, test->ncols);
  test->d = NULL;

  for (i = 0; i < test->nrows; i++)
  {
    for (j = 0; j < test->ncols; j++)
    {
      const gdouble d = g_test_rand_double_range (-1.0, 1.0);

      ncm_matrix_set (test->m, i, j, d);
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
    }
  }
}

void
test_ncm_matrix_new_gsl (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gsl_matrix *gm = gsl_matrix_alloc (test->nrows, test->ncols);
    guint i, j;

    test->m = ncm_matrix_new_gsl (gm);
    test->d = NULL;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);

        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_assert_true (ncm_matrix_nrows (test->m) == gm->size1 && ncm_matrix_ncols (test->m) == gm->size2);
  }
}

void
test_ncm_matrix_new_array (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    GArray *ga = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), test->nrows * test->ncols);
    guint i, j;

    g_array_set_size (ga, test->nrows * test->ncols);

    test->m = ncm_matrix_new_array (ga, test->ncols);
    test->d = NULL;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);


        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_array_unref (ga);

    g_assert_cmpuint (ncm_matrix_nrows (test->m), ==, ga->len / test->ncols);

    g_assert_true (ga == ncm_matrix_get_array (test->m));
    {
      GArray *ga_dup = ncm_matrix_dup_array (test->m);

      g_assert_true (ga_dup != ga);
      g_assert_cmpuint (ga_dup->len, ==, ga->len);

      for (i = 0; i < test->nrows * test->ncols; i++)
        ncm_assert_cmpdouble (g_array_index (ga, gdouble, i), ==, g_array_index (ga_dup, gdouble, i));

      g_array_unref (ga_dup);
    }

    g_array_unref (ga);
  }
}

void
test_ncm_matrix_new_data_slice (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gdouble *d = g_slice_alloc (test->nrows * test->ncols * sizeof (gdouble));
    guint i, j;

    test->m = ncm_matrix_new_data_slice (d, test->nrows, test->ncols);
    test->d = NULL;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);

        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_assert_true ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));
  }
}

void
test_ncm_matrix_new_data_malloc (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gdouble *d = g_malloc (test->nrows * test->ncols * sizeof (gdouble));
    guint i, j;

    test->m = ncm_matrix_new_data_malloc (d, test->nrows, test->ncols);
    test->d = NULL;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);

        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_assert_cmpuint ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)), ==, (test->nrows * test->ncols));
  }
}

void
test_ncm_matrix_new_data_static (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gdouble *d = g_new (gdouble, test->nrows * test->ncols);
    guint i, j;

    test->m = ncm_matrix_new_data_static (d, test->nrows, test->ncols);
    test->d = d;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);

        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_assert_true ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));
  }

  {
    NcmVector *v      = ncm_matrix_as_vector (test->m);
    const guint ncols = ncm_matrix_ncols (test->m);
    guint i;

    g_assert_true (ncm_vector_len (v) == ncols * ncm_matrix_nrows (test->m));

    for (i = 0; i < ncm_vector_len (v); i++)
      ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_matrix_get (test->m, i / ncols, i % ncols));

    ncm_vector_free (v);
  }
}

void
test_ncm_matrix_new_data_static_tda (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gdouble *d = g_new (gdouble, test->nrows * test->ncols * 2);
    guint i, j;

    test->m = ncm_matrix_new_data_static_tda (d, test->nrows, test->ncols, test->ncols * 2);
    test->d = d;

    for (i = 0; i < test->nrows; i++)
    {
      for (j = 0; j < test->ncols; j++)
      {
        const gdouble d = g_test_rand_double_range (-1.0, 1.0);

        ncm_matrix_set (test->m, i, j, d);
        ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
      }
    }

    g_assert_true ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));

    for (i = 0; i < 10 * test->nrows; i++)
    {
      const guint nr    = g_test_rand_int_range (0, test->nrows);
      const guint nc    = g_test_rand_int_range (0, test->ncols);
      const gdouble val = g_test_rand_double ();

      ncm_matrix_set (test->m, nr, nc, val);

      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, val);
    }

    {
      NcmVector *v = ncm_matrix_get_col (test->m, test->ncols - 1);

      g_assert_true (ncm_vector_len (v) == ncm_matrix_nrows (test->m));

      for (i = 0; i < test->nrows; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_matrix_get (test->m, i, test->ncols - 1));
      }

      ncm_vector_free (v);
    }

    {
      NcmVector *v = ncm_matrix_get_row (test->m, test->nrows - 1);

      g_assert_true (ncm_vector_len (v) == ncm_matrix_ncols (test->m));

      for (i = 0; i < test->ncols; i++)
        g_assert_true (ncm_vector_get (v, i) == ncm_matrix_get (test->m, test->nrows - 1, i));

      ncm_vector_free (v);
    }
  }
}

void
test_ncm_matrix_sanity (TestNcmMatrix *test, gconstpointer pdata)
{
  guint i, j;

  g_assert_true (NCM_IS_MATRIX (test->m));

  for (i = 0; i < test->nrows; i++)
  {
    for (j = 0; j < test->ncols; j++)
    {
      const gdouble d = g_test_rand_double ();

      ncm_matrix_set (test->m, i, j, d);
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, d);
    }
  }

  for (i = 0; i < 10 * test->nrows; i++)
  {
    const guint nr  = g_test_rand_int_range (0, test->nrows);
    const guint nc  = g_test_rand_int_range (0, test->ncols);
    const gdouble d = g_test_rand_double ();

    ncm_matrix_set (test->m, nr, nc, d);

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, d);
  }
}

void
test_ncm_matrix_operations (TestNcmMatrix *test, gconstpointer pdata)
{
  NcmMatrix *cm = ncm_matrix_dup (test->m);
  guint i;

  test_ncm_matrix_sanity (test, pdata);

  for (i = 0; i < 10 * test->nrows; i++)
  {
    const guint nr  = g_test_rand_int_range (0, test->nrows);
    const guint nc  = g_test_rand_int_range (0, test->ncols);
    const gdouble d = g_test_rand_double ();

    ncm_matrix_set (test->m, nr, nc, d);

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, d);
  }

  for (i = 0; i < 10 * test->nrows; i++)
  {
    guint nr, nc;
    gdouble *d;

    nr    = g_test_rand_int_range (0, test->nrows);
    nc    = g_test_rand_int_range (0, test->ncols);
    d     = ncm_matrix_ptr (test->m, nr, nc);
    (*d) *= g_test_rand_double ();

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, *d);
  }

  {
    const gint nrows = (test->nrows <= test->ncols) ? 0 : test->nrows - test->ncols;
    const gint ncols = (test->ncols <= test->nrows) ? 0 : test->ncols - test->nrows;
    NcmMatrix *subm  = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
    NcmMatrix *sm_cp = ncm_matrix_dup (subm);

    ncm_matrix_transpose (subm);

    for (i = 0; i < 10 * test->nrows; i++)
    {
      guint nr, nc;
      gdouble d;

      nr = g_test_rand_int_range (0, test->nrows - nrows);
      nc = g_test_rand_int_range (0, test->ncols - ncols);
      d  = ncm_matrix_get (sm_cp, nr, nc);

      ncm_assert_cmpdouble (ncm_matrix_get (subm, nc, nr), ==, d);
    }

    ncm_matrix_free (sm_cp);
    ncm_matrix_free (subm);
  }

  for (i = 0; i < 10 * test->nrows; i++)
  {
    guint nr = g_test_rand_int_range (0, test->nrows);
    guint nc = g_test_rand_int_range (0, test->ncols);

    ncm_matrix_set_identity (test->m);

    if (nr == nc)
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, 1.0);
    else
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, 0.0);
  }

  {
    ncm_matrix_set_zero (test->m);

    for (i = 0; i < 10 * test->nrows; i++)
    {
      guint nr = g_test_rand_int_range (0, test->nrows);
      guint nc = g_test_rand_int_range (0, test->ncols);

      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, 0.0);
    }
  }

  for (i = 0; i < 10 * test->nrows; i++)
  {
    guint nr, nc;
    gdouble d, d1 = g_test_rand_double ();

    nr = g_test_rand_int_range (0, test->nrows);
    nc = g_test_rand_int_range (0, test->ncols);
    d  = ncm_matrix_get (test->m, nr, nc);
    ncm_matrix_scale (test->m, d1);

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, d * d1);
  }

  {
    NcmVector *v = ncm_vector_new (test->nrows);
    guint nc     = g_test_rand_int_range (0, test->ncols);

    for (i = 0; i < test->nrows; i++)
    {
      ncm_vector_set (v, i, 8.2 + i);
    }

    ncm_matrix_set_col (test->m, nc, v);

    for (i = 0; i < test->nrows; i++)
    {
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, nc), ==, ncm_vector_get (v, i));
    }

    ncm_vector_free (v);
  }

  ncm_matrix_memcpy (test->m, cm);

  for (i = 0; i < test->nrows; i++)
  {
    guint j;

    for (j = 0; j < test->ncols; j++)
    {
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, ncm_matrix_get (cm, i, j));
    }
  }

  NCM_TEST_FREE (ncm_matrix_free, cm);
}

void
test_ncm_matrix_colmajor (TestNcmMatrix *test, gconstpointer pdata)
{
  NcmMatrix *cm0 = ncm_matrix_dup (test->m);
  NcmMatrix *cm1 = ncm_matrix_dup (test->m);
  guint i;

  test_ncm_matrix_sanity (test, pdata);

  ncm_matrix_memcpy_to_colmajor (cm1, cm0);

  for (i = 0; i < test->nrows; i++)
  {
    guint j;

    for (j = 0; j < test->ncols; j++)
      ncm_assert_cmpdouble (ncm_matrix_get_colmajor (cm1, i, j), ==, ncm_matrix_get (cm0, i, j));
  }

  for (i = 0; i < test->nrows; i++)
  {
    guint j;

    for (j = 0; j < test->ncols; j++)
    {
      const gdouble m_ij = g_test_rand_double_range (-1.0, 1.0);

      ncm_matrix_set_colmajor (cm0, i, j, m_ij);

      ncm_assert_cmpdouble (ncm_matrix_get_colmajor (cm0, i, j), ==, m_ij);
    }
  }

  NCM_TEST_FREE (ncm_matrix_free, cm0);
  NCM_TEST_FREE (ncm_matrix_free, cm1);
}

void
test_ncm_matrix_submatrix (TestNcmMatrix *test, gconstpointer pdata)
{
  const gint nrows = g_test_rand_int_range (0, test->nrows);
  const gint ncols = g_test_rand_int_range (0, test->ncols);
  NcmMatrix *sm    = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
  guint ntests     = 20 * test->nrows;

  g_assert_true (ncm_matrix_nrows (sm) == (test->nrows - nrows) && ncm_matrix_ncols (sm) == (test->ncols - ncols));

  while (ntests--)
  {
    guint nr = g_test_rand_int_range (0, test->nrows - nrows);
    guint nc = g_test_rand_int_range (0, test->ncols - ncols);

    ncm_matrix_set (sm, nr, nc, g_test_rand_double ());
    ncm_assert_cmpdouble (ncm_matrix_get (sm, nr, nc), ==, ncm_matrix_get (test->m, nr + nrows, nc + ncols));
  }

  g_assert_true (NCM_IS_MATRIX (sm));

  ncm_matrix_free (test->m);
  g_assert_true (G_IS_OBJECT (test->m));
  ncm_matrix_ref (test->m);

  NCM_TEST_FREE (ncm_matrix_free, sm);
}

void
test_ncm_matrix_add_mul (TestNcmMatrix *test, gconstpointer pdata)
{
  const gint nrows = g_test_rand_int_range (0, test->nrows);
  const gint ncols = g_test_rand_int_range (0, test->ncols);
  NcmMatrix *sm    = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
  guint i, j;

  g_assert_true (ncm_matrix_nrows (sm) == (test->nrows - nrows) && ncm_matrix_ncols (sm) == (test->ncols - ncols));

  ncm_matrix_set_zero (sm);

  for (i = 0; i < ncm_matrix_nrows (sm); i++)
  {
    for (j = 0; j < ncm_matrix_ncols (sm); j++)
    {
      ncm_matrix_set (sm, i, j, g_test_rand_double ());
      ncm_assert_cmpdouble (ncm_matrix_get (sm, i, j), ==, ncm_matrix_get (test->m, i + nrows, j + ncols));
    }
  }

  {
    NcmMatrix *osm      = ncm_matrix_dup (sm);
    NcmMatrix *res      = ncm_matrix_dup (sm);
    const gdouble alpha = g_test_rand_double ();

    for (i = 0; i < ncm_matrix_nrows (osm); i++)
    {
      for (j = 0; j < ncm_matrix_ncols (osm); j++)
      {
        ncm_matrix_set (osm, i, j, g_test_rand_double ());
      }
    }

    ncm_matrix_add_mul (res, alpha, osm);

    for (i = 0; i < ncm_matrix_nrows (osm); i++)
    {
      for (j = 0; j < ncm_matrix_ncols (osm); j++)
      {
        ncm_assert_cmpdouble (ncm_matrix_get (res, i, j), ==, ncm_matrix_get (sm, i, j) + alpha * ncm_matrix_get (osm, i, j));
      }
    }

    ncm_matrix_memcpy (res, osm);

    ncm_matrix_add_mul (res, alpha, sm);

    for (i = 0; i < ncm_matrix_nrows (sm); i++)
    {
      for (j = 0; j < ncm_matrix_ncols (sm); j++)
      {
        ncm_assert_cmpdouble (ncm_matrix_get (res, i, j), ==, ncm_matrix_get (osm, i, j) + alpha * ncm_matrix_get (sm, i, j));
      }
    }

    NCM_TEST_FREE (ncm_matrix_free, osm);
    NCM_TEST_FREE (ncm_matrix_free, res);
  }

  NCM_TEST_FREE (ncm_matrix_free, sm);
}

void
test_ncm_matrix_log_exp (TestNcmMatrix *test, gconstpointer pdata)
{
  const gint nrows = (test->nrows <= test->ncols) ? 0 : test->nrows - test->ncols;
  const gint ncols = (test->ncols <= test->nrows) ? 0 : test->ncols - test->nrows;
  NcmMatrix *sm    = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
  NcmMatrix *dup   = ncm_matrix_dup (sm);
  NcmMatrix *exp_m = ncm_matrix_dup (sm);
  NcmMatrix *exp_c = ncm_matrix_dup (sm);
  NcmMatrix *log_m = ncm_matrix_dup (sm);

  ncm_matrix_transpose (dup);

  ncm_matrix_add (sm, dup);
  ncm_matrix_scale (sm, 0.5);

  ncm_matrix_memcpy (dup, sm);

  /*printf ("\n"); ncm_matrix_log_vals (sm, "LOGM: ", "% 22.15g");*/

  ncm_matrix_sym_exp_cholesky (sm, 'U', exp_c);

  /*printf ("\n"); ncm_matrix_log_vals (exp_c, "EXPC: ", "% 22.15g");*/

  ncm_matrix_triang_to_sym (exp_c, 'U', TRUE, exp_m);

  /*printf ("\n"); ncm_matrix_log_vals (exp_m, "EXPM: ", "% 22.15g");*/

  ncm_matrix_sym_posdef_log (exp_m, 'U', log_m);
  ncm_matrix_copy_triangle (log_m, 'U');

  ncm_assert_cmpdouble_e (ncm_matrix_cmp (dup, log_m, 0.0), ==, 0.0, 1.0e-7, 1.0e-7);

  NCM_TEST_FREE (ncm_matrix_free, sm);
  NCM_TEST_FREE (ncm_matrix_free, dup);
  NCM_TEST_FREE (ncm_matrix_free, exp_m);
  NCM_TEST_FREE (ncm_matrix_free, exp_c);
  NCM_TEST_FREE (ncm_matrix_free, log_m);
}

/*
 * @UL names the triangle that survives: 'U' clears the strict lower and leaves the upper
 * and the diagonal untouched, 'L' the other way round.
 */
void
test_ncm_matrix_zero_triangle (TestNcmMatrix *test, gconstpointer pdata)
{
  gint tests;

  for (tests = 0; tests < 12; tests++)
  {
    const guint n  = g_test_rand_int_range (2, 40);
    const gchar UL = (g_test_rand_int_range (0, 2) == 0) ? 'U' : 'L';
    NcmMatrix *m   = ncm_matrix_new (n, n);
    NcmMatrix *ref = ncm_matrix_new (n, n);
    guint i, j;

    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        const gdouble v = g_test_rand_double_range (-9.0, 9.0);

        ncm_matrix_set (m, i, j, v);
        ncm_matrix_set (ref, i, j, v);
      }
    }

    ncm_matrix_zero_triangle (m, UL);

    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        const gboolean cleared = (UL == 'U') ? (j < i) : (j > i);

        if (cleared)
          g_assert_cmpfloat (ncm_matrix_get (m, i, j), ==, 0.0);
        else
          g_assert_cmpfloat (ncm_matrix_get (m, i, j), ==, ncm_matrix_get (ref, i, j));
      }
    }

    /* Clearing a triangle twice changes nothing more. */
    ncm_matrix_zero_triangle (m, UL);

    for (i = 0; i < n; i++)
      for (j = i; j < n; j++)
        if ((UL == 'U'))
          g_assert_cmpfloat (ncm_matrix_get (m, i, j), ==, ncm_matrix_get (ref, i, j));

    NCM_TEST_FREE (ncm_matrix_free, m);
    NCM_TEST_FREE (ncm_matrix_free, ref);
  }
}

/*
 * U^T U for an upper triangular U, and L L^T for a lower one, both against a direct
 * triple loop. Covers each triangle and both settings of @zero: with it off the caller
 * promises the unused triangle is already clear, with it on the routine clears it, and
 * the product must come out the same either way.
 */
void
test_ncm_matrix_triang_to_sym (TestNcmMatrix *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-13;
  gint tests;

  for (tests = 0; tests < 12; tests++)
  {
    const guint n    = g_test_rand_int_range (2, 40);
    const gchar UL   = (g_test_rand_int_range (0, 2) == 0) ? 'U' : 'L';
    NcmMatrix *T     = ncm_matrix_new (n, n);
    NcmMatrix *dirty = ncm_matrix_new (n, n);
    NcmMatrix *S0    = ncm_matrix_new (n, n);
    NcmMatrix *S1    = ncm_matrix_new (n, n);
    guint i, j, k;

    /* T is triangular with the unused triangle exactly zero; dirty is the same matrix
     * with rubbish in that triangle, which @zero is supposed to clear. */
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        const gboolean used = (UL == 'U') ? (j >= i) : (j <= i);
        const gdouble v     = g_test_rand_double_range (0.5, 2.0);

        ncm_matrix_set (T, i, j, used ? v : 0.0);
        ncm_matrix_set (dirty, i, j, used ? v : g_test_rand_double_range (-9.0, 9.0));
      }
    }

    ncm_matrix_triang_to_sym (T, UL, FALSE, S0);
    ncm_matrix_triang_to_sym (dirty, UL, TRUE, S1);

    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        gdouble ref = 0.0;

        for (k = 0; k < n; k++)
          ref += (UL == 'U') ? ncm_matrix_get (T, k, i) * ncm_matrix_get (T, k, j)
                 : ncm_matrix_get (T, i, k) * ncm_matrix_get (T, j, k);

        ncm_assert_cmpdouble_e (ncm_matrix_get (S0, i, j), ==, ref, reltol, 0.0);
        ncm_assert_cmpdouble_e (ncm_matrix_get (S1, i, j), ==, ref, reltol, 0.0);

        /* @zero also clears the unused triangle of its input. */
        if ((UL == 'U') ? (j < i) : (j > i))
          g_assert_cmpfloat (ncm_matrix_get (dirty, i, j), ==, 0.0);
      }
    }

    NCM_TEST_FREE (ncm_matrix_free, T);
    NCM_TEST_FREE (ncm_matrix_free, dirty);
    NCM_TEST_FREE (ncm_matrix_free, S0);
    NCM_TEST_FREE (ncm_matrix_free, S1);
  }
}

void
test_ncm_matrix_square_to_sym (TestNcmMatrix *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-14;
  const gdouble abstol = 0.0;
  gint tests;

  for (tests = 0; tests < 10; tests++)
  {
    gint nrows    = g_test_rand_int_range (2, 100);
    gint ncols    = g_test_rand_int_range (2, 100);
    NcmMatrix *A  = ncm_matrix_new (nrows, ncols);
    NcmMatrix *C0 = ncm_matrix_new (nrows, nrows);
    NcmMatrix *C1 = ncm_matrix_new (ncols, ncols);
    gint i, j, k;

    if (g_test_rand_int_range (0, 1))
    {
      NcmMatrix *subA, *subC0, *subC1;

      nrows = g_test_rand_int_range (2, nrows);
      ncols = g_test_rand_int_range (2, ncols);
      subA  = ncm_matrix_get_submatrix (A,  0, 0, nrows, ncols);
      subC0 = ncm_matrix_get_submatrix (C0, 0, 0, nrows, nrows);
      subC1 = ncm_matrix_get_submatrix (C1, 0, 0, ncols, ncols);

      ncm_matrix_clear (&A);
      ncm_matrix_clear (&C0);
      ncm_matrix_clear (&C1);

      A  = subA;
      C0 = subC0;
      C1 = subC1;
    }

    for (i = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++)
      {
        ncm_matrix_set (A, i, j, g_test_rand_double ());
      }
    }

    ncm_matrix_square_to_sym (A, 'N', 'U', C0);
    ncm_matrix_square_to_sym (A, 'T', 'U', C1);

    for (i = 0; i < nrows; i++)
    {
      for (j = i; j < nrows; j++)
      {
        gdouble C0_ij = 0.0;

        for (k = 0; k < ncols; k++)
        {
          C0_ij += ncm_matrix_get (A, i, k) * ncm_matrix_get (A, j, k);
        }

        ncm_assert_cmpdouble_e (ncm_matrix_get (C0, i, j), ==, C0_ij, reltol, abstol);
      }
    }

    for (i = 0; i < ncols; i++)
    {
      for (j = i; j < ncols; j++)
      {
        gdouble C1_ij = 0.0;

        for (k = 0; k < nrows; k++)
        {
          C1_ij += ncm_matrix_get (A, k, i) * ncm_matrix_get (A, k, j);
        }

        ncm_assert_cmpdouble_e (ncm_matrix_get (C1, i, j), ==, C1_ij, reltol, abstol);
      }
    }

    ncm_matrix_square_to_sym (A, 'N', 'L', C0);
    ncm_matrix_square_to_sym (A, 'T', 'L', C1);

    for (i = 0; i < nrows; i++)
    {
      for (j = 0; j <= i; j++)
      {
        gdouble C0_ij = 0.0;

        for (k = 0; k < ncols; k++)
        {
          C0_ij += ncm_matrix_get (A, i, k) * ncm_matrix_get (A, j, k);
        }

        ncm_assert_cmpdouble_e (ncm_matrix_get (C0, i, j), ==, C0_ij, reltol, abstol);
      }
    }

    for (i = 0; i < ncols; i++)
    {
      for (j = 0; j <= i; j++)
      {
        gdouble C1_ij = 0.0;

        for (k = 0; k < nrows; k++)
        {
          C1_ij += ncm_matrix_get (A, k, i) * ncm_matrix_get (A, k, j);
        }

        ncm_assert_cmpdouble_e (ncm_matrix_get (C1, i, j), ==, C1_ij, reltol, abstol);
      }
    }

    ncm_matrix_clear (&A);
    ncm_matrix_clear (&C0);
    ncm_matrix_clear (&C1);
  }
}

void
test_ncm_matrix_update_vector (TestNcmMatrix *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-14;
  const gdouble abstol = 0.0;
  gint tests;

  for (tests = 0; tests < 10; tests++)
  {
    gint nrows          = g_test_rand_int_range (2, 100);
    gint ncols          = g_test_rand_int_range (2, 100);
    NcmMatrix *A        = ncm_matrix_new (nrows, ncols);
    NcmVector *v        = ncm_vector_new (ncols);
    NcmVector *v_dup    = ncm_vector_new (ncols);
    NcmVector *u        = ncm_vector_new (nrows);
    NcmVector *u_dup    = ncm_vector_new (nrows);
    const gdouble alpha = g_test_rand_double ();
    const gdouble beta  = g_test_rand_double ();
    gint i, j;

    if (g_test_rand_int_range (0, 1))
    {
      NcmMatrix *subA;

      nrows = g_test_rand_int_range (2, nrows);
      ncols = g_test_rand_int_range (2, ncols);
      subA  = ncm_matrix_get_submatrix (A,  0, 0, nrows, ncols);

      ncm_matrix_clear (&A);

      A = subA;
    }

    for (j = 0; j < ncols; j++)
    {
      ncm_vector_set (v, j, g_test_rand_double ());

      for (i = 0; i < nrows; i++)
      {
        ncm_matrix_set (A, i, j, g_test_rand_double ());

        if (j == 0)
          ncm_vector_set (u, i, g_test_rand_double ());
      }
    }

    ncm_vector_memcpy (u_dup, u);
    ncm_matrix_update_vector (A, 'N', alpha, v, beta, u_dup);

    for (i = 0; i < nrows; i++)
    {
      gdouble Av_i = 0.0;

      for (j = 0; j < ncols; j++)
      {
        Av_i += ncm_matrix_get (A, i, j) * ncm_vector_get (v, j);
      }

      Av_i = alpha * Av_i + beta * ncm_vector_get (u, i);

      ncm_assert_cmpdouble_e (ncm_vector_get (u_dup, i), ==, Av_i, reltol, abstol);
    }

    ncm_vector_memcpy (v_dup, v);
    ncm_matrix_update_vector (A, 'T', alpha, u, beta, v_dup);

    for (i = 0; i < ncols; i++)
    {
      gdouble Au_i = 0.0;

      for (j = 0; j < nrows; j++)
      {
        Au_i += ncm_matrix_get (A, j, i) * ncm_vector_get (u, j);
      }

      Au_i = alpha * Au_i + beta * ncm_vector_get (v, i);

      ncm_assert_cmpdouble_e (ncm_vector_get (v_dup, i), ==, Au_i, reltol, abstol);
    }

    ncm_matrix_clear (&A);
    ncm_vector_clear (&v);
    ncm_vector_clear (&v_dup);
    ncm_vector_clear (&u);
    ncm_vector_clear (&u_dup);
  }
}

void
test_ncm_matrix_sym_update_vector (TestNcmMatrix *test, gconstpointer pdata)
{
  const gdouble reltol = 1.0e-14;
  const gdouble abstol = 0.0;
  gint tests;

  for (tests = 0; tests < 10; tests++)
  {
    gint n              = g_test_rand_int_range (2, 100);
    NcmMatrix *A        = ncm_matrix_new (n, n);
    NcmVector *v        = ncm_vector_new (n);
    NcmVector *u        = ncm_vector_new (n);
    NcmVector *u_dup    = ncm_vector_new (n);
    const gdouble alpha = g_test_rand_double ();
    const gdouble beta  = g_test_rand_double ();
    gchar Uplo          = g_test_rand_int_range (0, 1) ? 'U' : 'L';
    gint i, j;

    if (g_test_rand_int_range (0, 1))
    {
      NcmMatrix *subA;

      n    = g_test_rand_int_range (2, n);
      subA = ncm_matrix_get_submatrix (A,  0, 0, n, n);

      ncm_matrix_clear (&A);

      A = subA;
    }

    for (j = 0; j < n; j++)
    {
      ncm_vector_set (v, j, g_test_rand_double ());

      for (i = 0; i < n; i++)
      {
        ncm_matrix_set (A, i, j, g_test_rand_double ());

        if (j == 0)
          ncm_vector_set (u, i, g_test_rand_double ());
      }
    }

    ncm_vector_memcpy (u_dup, u);
    ncm_matrix_sym_update_vector (A, Uplo, alpha, v, beta, u_dup);

    for (i = 0; i < n; i++)
    {
      gdouble Av_i = 0.0;

      if (Uplo == 'U')
      {
        for (j = 0; j < i; j++)
        {
          Av_i += ncm_matrix_get (A, j, i) * ncm_vector_get (v, j);
        }

        for (j = i; j < n; j++)
        {
          Av_i += ncm_matrix_get (A, i, j) * ncm_vector_get (v, j);
        }
      }
      else
      {
        for (j = 0; j < i; j++)
        {
          Av_i += ncm_matrix_get (A, i, j) * ncm_vector_get (v, j);
        }

        for (j = i; j < n; j++)
        {
          Av_i += ncm_matrix_get (A, j, i) * ncm_vector_get (v, j);
        }
      }

      Av_i = alpha * Av_i + beta * ncm_vector_get (u, i);

      ncm_assert_cmpdouble_e (ncm_vector_get (u_dup, i), ==, Av_i, reltol, abstol);
    }

    ncm_matrix_clear (&A);
    ncm_vector_clear (&v);
    ncm_vector_clear (&u);
    ncm_vector_clear (&u_dup);
  }
}

/* op(A) as a dense matrix: the named triangle of A, the other one zero, transposed on 'T'. */
static void
_test_ncm_matrix_triangle_op (NcmMatrix *A, gchar UL, gchar T, NcmMatrix *opA)
{
  const guint n = ncm_matrix_nrows (A);
  guint i, j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      const gboolean in_tri = (UL == 'U') ? (j >= i) : (j <= i);
      const gdouble a_ij    = in_tri ? ncm_matrix_get (A, i, j) : 0.0;

      if (T == 'N')
        ncm_matrix_set (opA, i, j, a_ij);
      else
        ncm_matrix_set (opA, j, i, a_ij);
    }
  }
}

/* Both triangles filled, so the one a routine must not read is never zero; the diagonal
 * dominates, so the solves stay well conditioned. */
static void
_test_ncm_matrix_fill_triangular (NcmMatrix *A)
{
  const guint n = ncm_matrix_nrows (A);
  guint i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      ncm_matrix_set (A, i, j, g_test_rand_double_range (-1.0, 1.0) + ((i == j) ? 2.0 * n : 0.0));
}

void
test_ncm_matrix_dtrmm_dtrsm (TestNcmMatrix *test, gconstpointer pdata)
{
  const gchar Sides[2] = {'L', 'R'};
  const gchar ULs[2]   = {'U', 'L'};
  const gchar Ts[2]    = {'N', 'T'};
  gint tests;

  for (tests = 0; tests < 6; tests++)
  {
    const guint n  = g_test_rand_int_range (2, 30);
    const guint m  = g_test_rand_int_range (2, 30);
    NcmMatrix *A   = ncm_matrix_new (n, n);
    NcmMatrix *opA = ncm_matrix_new (n, n);
    guint si, ui, ti;

    _test_ncm_matrix_fill_triangular (A);

    for (si = 0; si < 2; si++)
    {
      for (ui = 0; ui < 2; ui++)
      {
        for (ti = 0; ti < 2; ti++)
        {
          const gchar Side    = Sides[si];
          const gchar UL      = ULs[ui];
          const gchar T       = Ts[ti];
          const guint nrows   = (Side == 'L') ? n : m;
          const guint ncols   = (Side == 'L') ? m : n;
          const gdouble alpha = g_test_rand_double_range (0.5, 2.0);
          NcmMatrix *B0       = ncm_matrix_new (nrows, ncols);
          NcmMatrix *B        = ncm_matrix_new (nrows, ncols);
          guint i, j, k;

          for (i = 0; i < nrows; i++)
            for (j = 0; j < ncols; j++)
              ncm_matrix_set (B0, i, j, g_test_rand_double_range (-1.0, 1.0));

          ncm_matrix_memcpy (B, B0);
          _test_ncm_matrix_triangle_op (A, UL, T, opA);

          ncm_matrix_dtrmm (B, Side, UL, T, alpha, A);

          for (i = 0; i < nrows; i++)
          {
            for (j = 0; j < ncols; j++)
            {
              gdouble ref = 0.0;

              for (k = 0; k < n; k++)
              {
                if (Side == 'L')
                  ref += ncm_matrix_get (opA, i, k) * ncm_matrix_get (B0, k, j);
                else
                  ref += ncm_matrix_get (B0, i, k) * ncm_matrix_get (opA, k, j);
              }

              ncm_assert_cmpdouble_e (ncm_matrix_get (B, i, j), ==, alpha * ref, 1.0e-12, 1.0e-12);
            }
          }

          /* The solve undoes the product. */
          ncm_matrix_dtrsm (B, Side, UL, T, 1.0 / alpha, A);

          for (i = 0; i < nrows; i++)
            for (j = 0; j < ncols; j++)
              ncm_assert_cmpdouble_e (ncm_matrix_get (B, i, j), ==, ncm_matrix_get (B0, i, j), 1.0e-11, 1.0e-12);

          NCM_TEST_FREE (ncm_matrix_free, B0);
          NCM_TEST_FREE (ncm_matrix_free, B);
        }
      }
    }

    NCM_TEST_FREE (ncm_matrix_free, A);
    NCM_TEST_FREE (ncm_matrix_free, opA);
  }
}

void
test_ncm_matrix_dtrmv_dtrsv (TestNcmMatrix *test, gconstpointer pdata)
{
  const gchar ULs[2] = {'U', 'L'};
  const gchar Ts[2]  = {'N', 'T'};
  gint tests;

  for (tests = 0; tests < 6; tests++)
  {
    const guint n  = g_test_rand_int_range (2, 30);
    NcmMatrix *A   = ncm_matrix_new (n, n);
    NcmMatrix *opA = ncm_matrix_new (n, n);
    NcmVector *v0  = ncm_vector_new (n);
    NcmVector *v   = ncm_vector_new (n);
    guint ui, ti;

    _test_ncm_matrix_fill_triangular (A);

    for (ui = 0; ui < 2; ui++)
    {
      for (ti = 0; ti < 2; ti++)
      {
        const gchar UL = ULs[ui];
        const gchar T  = Ts[ti];
        guint i, k;

        for (i = 0; i < n; i++)
          ncm_vector_set (v0, i, g_test_rand_double_range (-1.0, 1.0));

        ncm_vector_memcpy (v, v0);
        _test_ncm_matrix_triangle_op (A, UL, T, opA);

        ncm_matrix_dtrmv (A, UL, T, v);

        for (i = 0; i < n; i++)
        {
          gdouble ref = 0.0;

          for (k = 0; k < n; k++)
            ref += ncm_matrix_get (opA, i, k) * ncm_vector_get (v0, k);

          ncm_assert_cmpdouble_e (ncm_vector_get (v, i), ==, ref, 1.0e-12, 1.0e-12);
        }

        ncm_matrix_dtrsv (A, UL, T, v);

        for (i = 0; i < n; i++)
          ncm_assert_cmpdouble_e (ncm_vector_get (v, i), ==, ncm_vector_get (v0, i), 1.0e-11, 1.0e-12);
      }
    }

    NCM_TEST_FREE (ncm_matrix_free, A);
    NCM_TEST_FREE (ncm_matrix_free, opA);
    NCM_TEST_FREE (ncm_vector_free, v0);
    NCM_TEST_FREE (ncm_vector_free, v);
  }
}

void
test_ncm_matrix_dsyrk (TestNcmMatrix *test, gconstpointer pdata)
{
  const gchar ULs[2] = {'U', 'L'};
  const gchar Ts[2]  = {'N', 'T'};
  gint tests;

  for (tests = 0; tests < 6; tests++)
  {
    const guint n = g_test_rand_int_range (2, 30);
    const guint k = g_test_rand_int_range (1, 30);
    NcmMatrix *C0 = ncm_matrix_new (n, n);
    NcmMatrix *C  = ncm_matrix_new (n, n);
    guint ui, ti;

    for (ui = 0; ui < 2; ui++)
    {
      for (ti = 0; ti < 2; ti++)
      {
        const gchar UL      = ULs[ui];
        const gchar T       = Ts[ti];
        const gdouble alpha = g_test_rand_double_range (0.5, 2.0);
        const gdouble beta  = g_test_rand_double_range (-1.0, 1.0);
        NcmMatrix *A        = (T == 'N') ? ncm_matrix_new (n, k) : ncm_matrix_new (k, n);
        guint i, j, l;

        for (i = 0; i < ncm_matrix_nrows (A); i++)
          for (j = 0; j < ncm_matrix_ncols (A); j++)
            ncm_matrix_set (A, i, j, g_test_rand_double_range (-1.0, 1.0));

        for (i = 0; i < n; i++)
          for (j = 0; j < n; j++)
            ncm_matrix_set (C0, i, j, g_test_rand_double_range (-1.0, 1.0));

        ncm_matrix_memcpy (C, C0);
        ncm_matrix_dsyrk (C, UL, T, alpha, A, beta);

        for (i = 0; i < n; i++)
        {
          for (j = 0; j < n; j++)
          {
            const gboolean in_tri = (UL == 'U') ? (j >= i) : (j <= i);

            if (in_tri)
            {
              gdouble ref = 0.0;

              for (l = 0; l < k; l++)
              {
                const gdouble a_il = (T == 'N') ? ncm_matrix_get (A, i, l) : ncm_matrix_get (A, l, i);
                const gdouble a_jl = (T == 'N') ? ncm_matrix_get (A, j, l) : ncm_matrix_get (A, l, j);

                ref += a_il * a_jl;
              }

              ref = alpha * ref + beta * ncm_matrix_get (C0, i, j);
              ncm_assert_cmpdouble_e (ncm_matrix_get (C, i, j), ==, ref, 1.0e-12, 1.0e-12);
            }
            else
            {
              /* The other triangle is not touched. */
              g_assert_cmpfloat (ncm_matrix_get (C, i, j), ==, ncm_matrix_get (C0, i, j));
            }
          }
        }

        NCM_TEST_FREE (ncm_matrix_free, A);
      }
    }

    NCM_TEST_FREE (ncm_matrix_free, C0);
    NCM_TEST_FREE (ncm_matrix_free, C);
  }
}

void
test_ncm_matrix_scale_rows_cols (TestNcmMatrix *test, gconstpointer pdata)
{
  gint tests;

  for (tests = 0; tests < 8; tests++)
  {
    const guint nrows = g_test_rand_int_range (1, 30);
    const guint ncols = g_test_rand_int_range (1, 30);
    NcmMatrix *M0     = ncm_matrix_new (nrows, ncols);
    NcmMatrix *M      = ncm_matrix_new (nrows, ncols);
    NcmVector *r      = ncm_vector_new (nrows);
    NcmVector *c      = ncm_vector_new (ncols);
    guint i, j;

    for (i = 0; i < nrows; i++)
    {
      ncm_vector_set (r, i, g_test_rand_double_range (-2.0, 2.0));

      for (j = 0; j < ncols; j++)
        ncm_matrix_set (M0, i, j, g_test_rand_double_range (-1.0, 1.0));
    }

    for (j = 0; j < ncols; j++)
      ncm_vector_set (c, j, g_test_rand_double_range (-2.0, 2.0));

    ncm_matrix_memcpy (M, M0);
    ncm_matrix_scale_rows (M, r);

    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        g_assert_cmpfloat (ncm_matrix_get (M, i, j), ==, ncm_vector_get (r, i) * ncm_matrix_get (M0, i, j));

    ncm_matrix_memcpy (M, M0);
    ncm_matrix_scale_cols (M, c);

    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        g_assert_cmpfloat (ncm_matrix_get (M, i, j), ==, ncm_matrix_get (M0, i, j) * ncm_vector_get (c, j));

    NCM_TEST_FREE (ncm_matrix_free, M0);
    NCM_TEST_FREE (ncm_matrix_free, M);
    NCM_TEST_FREE (ncm_vector_free, r);
    NCM_TEST_FREE (ncm_vector_free, c);
  }
}

void
test_ncm_matrix_sub_row_vector (TestNcmMatrix *test, gconstpointer pdata)
{
  gint tests;

  for (tests = 0; tests < 8; tests++)
  {
    const guint nrows = g_test_rand_int_range (1, 30);
    const guint ncols = g_test_rand_int_range (1, 30);
    NcmMatrix *M0     = ncm_matrix_new (nrows, ncols);
    NcmMatrix *M      = ncm_matrix_new (nrows, ncols);
    NcmVector *v      = ncm_vector_new (ncols);
    NcmVector *v2     = ncm_vector_new (2 * ncols);
    NcmVector *vs     = ncm_vector_get_subvector_stride (v2, 0, ncols, 2);
    guint i, j;

    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        ncm_matrix_set (M0, i, j, g_test_rand_double_range (-1.0, 1.0));

    for (j = 0; j < ncols; j++)
    {
      ncm_vector_set (v, j, g_test_rand_double_range (-2.0, 2.0));
      ncm_vector_set (v2, 2 * j, ncm_vector_get (v, j));
      ncm_vector_set (v2, 2 * j + 1, 1.0e3);
    }

    ncm_matrix_memcpy (M, M0);
    ncm_matrix_sub_row_vector (M, v);

    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        g_assert_cmpfloat (ncm_matrix_get (M, i, j), ==, ncm_matrix_get (M0, i, j) - ncm_vector_get (v, j));

    /* A strided vector is read through its stride. */
    ncm_matrix_memcpy (M, M0);
    ncm_matrix_sub_row_vector (M, vs);

    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        g_assert_cmpfloat (ncm_matrix_get (M, i, j), ==, ncm_matrix_get (M0, i, j) - ncm_vector_get (v, j));

    NCM_TEST_FREE (ncm_matrix_free, M0);
    NCM_TEST_FREE (ncm_matrix_free, M);
    NCM_TEST_FREE (ncm_vector_free, vs);
    NCM_TEST_FREE (ncm_vector_free, v2);
    NCM_TEST_FREE (ncm_vector_free, v);
  }
}

void
test_ncm_matrix_chol_chi2_cols (TestNcmMatrix *test, gconstpointer pdata)
{
  const guint d_a[4]  = { 1, 3, 8, 50 };
  const guint np_a[3] = { 1, 7, 256 };
  const guint nb_a[4] = { 1, 7, 32, 256 };
  NcmRNG *rng         = ncm_rng_seeded_new (NULL, 20260923);
  guint i_d, i_np, i_nb;

  for (i_d = 0; i_d < 4; i_d++)
  {
    const guint d    = d_a[i_d];
    NcmMatrix *cov   = ncm_matrix_new (d, d);
    NcmMatrix *U     = ncm_matrix_new (d, d);
    NcmVector *mu    = ncm_vector_new (d);
    NcmVector *theta = ncm_vector_new (d);
    guint i, j;

    /* @mu scales the covariance, so it is an input here; a large correlation level keeps
     * the factor well conditioned up to d = 50, where a small one is singular. */
    ncm_vector_set_all (mu, 1.0);
    ncm_matrix_fill_rand_cov2 (cov, mu, 0.1, 2.0, 10.0, rng);
    ncm_matrix_memcpy (U, cov);
    g_assert_cmpint (ncm_matrix_cholesky_decomp (U, 'U'), ==, 0);

    /* Only the upper triangle may be read: the lower one is filled with garbage. */
    for (i = 1; i < d; i++)
      for (j = 0; j < i; j++)
        ncm_matrix_set (U, i, j, 1.0e10 * ncm_rng_uniform_gen (rng, -1.0, 1.0));

    for (j = 0; j < d; j++)
      ncm_vector_set (theta, j, ncm_rng_uniform_gen (rng, -2.0, 2.0));

    for (i_np = 0; i_np < 3; i_np++)
    {
      const guint np  = np_a[i_np];
      NcmMatrix *X    = ncm_matrix_new (d, np);
      NcmMatrix *Xr   = ncm_matrix_new (np, d);
      NcmVector *ref  = ncm_vector_new (np);
      NcmVector *chi2 = ncm_vector_new (np);
      NcmVector *dq   = ncm_vector_new_data_static (ncm_matrix_ptr (Xr, 0, 0), d, 1);
      guint p;

      for (p = 0; p < np; p++)
        for (j = 0; j < d; j++)
        {
          const gdouble x_pj = ncm_rng_uniform_gen (rng, -3.0, 3.0);

          ncm_matrix_set (X, j, p, x_pj);
          ncm_matrix_set (Xr, p, j, x_pj);
        }

      /* Reference: the sequence the VKDE evaluator used, on points as rows. */
      ncm_matrix_sub_row_vector (Xr, theta);
      ncm_matrix_dtrsm (Xr, 'R', 'U', 'N', 1.0, U);

      for (p = 0; p < np; p++)
      {
        ncm_vector_replace_data (dq, ncm_matrix_ptr (Xr, p, 0));
        ncm_vector_set (ref, p, ncm_vector_dot (dq, dq));
      }

      for (i_nb = 0; i_nb < 4; i_nb++)
      {
        NcmMatrix *work = ncm_matrix_new (d, nb_a[i_nb]);

        ncm_vector_set_all (chi2, -1.0);
        ncm_matrix_chol_chi2_cols (X, theta, U, work, chi2);

        for (p = 0; p < np; p++)
          ncm_assert_cmpdouble_e (ncm_vector_get (chi2, p), ==, ncm_vector_get (ref, p), 1.0e-13, 0.0);

        NCM_TEST_FREE (ncm_matrix_free, work);
      }

      /* A strided centre is read through its stride. */
      {
        NcmMatrix *work = ncm_matrix_new (d, 32);
        NcmVector *t2   = ncm_vector_new (2 * d);
        NcmVector *ts   = ncm_vector_get_subvector_stride (t2, 0, d, 2);

        for (j = 0; j < d; j++)
        {
          ncm_vector_set (t2, 2 * j, ncm_vector_get (theta, j));
          ncm_vector_set (t2, 2 * j + 1, 1.0e3);
        }

        ncm_vector_set_all (chi2, -1.0);
        ncm_matrix_chol_chi2_cols (X, ts, U, work, chi2);

        for (p = 0; p < np; p++)
          ncm_assert_cmpdouble_e (ncm_vector_get (chi2, p), ==, ncm_vector_get (ref, p), 1.0e-13, 0.0);

        NCM_TEST_FREE (ncm_vector_free, ts);
        NCM_TEST_FREE (ncm_vector_free, t2);
        NCM_TEST_FREE (ncm_matrix_free, work);
      }

      /* Submatrix inputs, whose row stride exceeds the column count. */
      {
        NcmMatrix *Xb = ncm_matrix_new (d + 2, np + 3);
        NcmMatrix *Xs = ncm_matrix_get_submatrix (Xb, 0, 0, d, np);
        NcmMatrix *Wb = ncm_matrix_new (d + 2, 32 + 5);
        NcmMatrix *Ws = ncm_matrix_get_submatrix (Wb, 0, 0, d, 32);

        ncm_matrix_set_all (Xb, 1.0e10);
        ncm_matrix_memcpy (Xs, X);

        ncm_vector_set_all (chi2, -1.0);
        ncm_matrix_chol_chi2_cols (Xs, theta, U, Ws, chi2);

        for (p = 0; p < np; p++)
          ncm_assert_cmpdouble_e (ncm_vector_get (chi2, p), ==, ncm_vector_get (ref, p), 1.0e-13, 0.0);

        NCM_TEST_FREE (ncm_matrix_free, Ws);
        NCM_TEST_FREE (ncm_matrix_free, Wb);
        NCM_TEST_FREE (ncm_matrix_free, Xs);
        NCM_TEST_FREE (ncm_matrix_free, Xb);
      }

      NCM_TEST_FREE (ncm_vector_free, dq);
      NCM_TEST_FREE (ncm_vector_free, chi2);
      NCM_TEST_FREE (ncm_vector_free, ref);
      NCM_TEST_FREE (ncm_matrix_free, Xr);
      NCM_TEST_FREE (ncm_matrix_free, X);
    }

    NCM_TEST_FREE (ncm_vector_free, theta);
    NCM_TEST_FREE (ncm_vector_free, mu);
    NCM_TEST_FREE (ncm_matrix_free, U);
    NCM_TEST_FREE (ncm_matrix_free, cov);
  }

  ncm_rng_free (rng);
}

void
test_ncm_matrix_is_identity (TestNcmMatrix *test, gconstpointer pdata)
{
  const gdouble tol = 1.0e-12;
  gint tests;

  for (tests = 0; tests < 8; tests++)
  {
    const guint n = g_test_rand_int_range (1, 30);
    const guint i = g_test_rand_int_range (0, n);
    const guint j = g_test_rand_int_range (0, n);
    NcmMatrix *m  = ncm_matrix_new (n, n);

    ncm_matrix_set_identity (m);
    g_assert_true (ncm_matrix_is_identity (m, tol));

    /* Inside the tolerance it is still the identity; at twice the tolerance it is not,
     * on the diagonal as much as off it. */
    ncm_matrix_addto (m, i, j, 0.5 * tol);
    g_assert_true (ncm_matrix_is_identity (m, tol));

    ncm_matrix_addto (m, i, j, 1.5 * tol);
    g_assert_false (ncm_matrix_is_identity (m, tol));

    ncm_matrix_set_identity (m);
    ncm_matrix_set (m, i, j, GSL_NAN);
    g_assert_false (ncm_matrix_is_identity (m, tol));

    NCM_TEST_FREE (ncm_matrix_free, m);
  }
}

void
test_ncm_matrix_cholesky_decomp_nearPD (TestNcmMatrix *test, gconstpointer pdata)
{
  const gchar ULs[2] = {'U', 'L'};
  const gdouble eps  = 1.0e-3;
  gint tests;

  for (tests = 0; tests < 6; tests++)
  {
    const guint n   = g_test_rand_int_range (2, 20);
    const gchar UL  = ULs[tests % 2];
    NcmMatrix *B    = ncm_matrix_new (n, n - 1);
    NcmMatrix *cov  = ncm_matrix_new (n, n);
    NcmMatrix *ref  = ncm_matrix_new (n, n);
    NcmMatrix *dec  = ncm_matrix_new (n, n);
    NcmMatrix *back = ncm_matrix_new (n, n);
    gboolean repaired;
    guint i, j;

    /* Positive definite: the plain factor, no repair. */
    for (i = 0; i < n; i++)
      for (j = 0; j < n - 1; j++)
        ncm_matrix_set (B, i, j, g_test_rand_double_range (-1.0, 1.0));

    ncm_matrix_dgemm (cov, 'N', 'T', 1.0, B, B, 0.0);

    for (i = 0; i < n; i++)
      ncm_matrix_addto (cov, i, i, 1.0);

    ncm_matrix_memcpy (ref, cov);
    g_assert_cmpint (ncm_matrix_cholesky_decomp (ref, UL), ==, 0);

    g_assert_cmpint (ncm_matrix_cholesky_decomp_nearPD (cov, dec, UL, 200, &repaired), ==, 0);
    g_assert_false (repaired);

    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        if ((UL == 'U') ? (j >= i) : (j <= i))
          g_assert_cmpfloat (ncm_matrix_get (dec, i, j), ==, ncm_matrix_get (ref, i, j));

    /* Indefinite by eps along one direction (an exactly singular matrix can still pass
     * dpotrf through rounding): refused without the repair, repaired with it, the input
     * untouched, and the repaired factor reproduces the matrix up to the eps that had to
     * be added back. */
    ncm_matrix_dgemm (cov, 'N', 'T', 1.0, B, B, 0.0);

    for (i = 0; i < n; i++)
      ncm_matrix_addto (cov, i, i, -eps);

    ncm_matrix_memcpy (ref, cov);

    g_assert_cmpint (ncm_matrix_cholesky_decomp_nearPD (cov, dec, UL, 0, &repaired), !=, 0);
    g_assert_false (repaired);

    g_assert_cmpint (ncm_matrix_cholesky_decomp_nearPD (cov, dec, UL, 200, &repaired), ==, 0);
    g_assert_true (repaired);
    g_assert_cmpfloat (ncm_matrix_cmp (cov, ref, 0.0), ==, 0.0);

    ncm_matrix_triang_to_sym (dec, UL, TRUE, back);

    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
        ncm_assert_cmpdouble_e (ncm_matrix_get (back, i, j), ==, ncm_matrix_get (ref, i, j), 0.0, 2.0 * eps);

    /* NULL is accepted for the flag. */
    g_assert_cmpint (ncm_matrix_cholesky_decomp_nearPD (cov, dec, UL, 200, NULL), ==, 0);

    NCM_TEST_FREE (ncm_matrix_free, B);
    NCM_TEST_FREE (ncm_matrix_free, cov);
    NCM_TEST_FREE (ncm_matrix_free, ref);
    NCM_TEST_FREE (ncm_matrix_free, dec);
    NCM_TEST_FREE (ncm_matrix_free, back);
  }
}

void
test_ncm_matrix_serialization (TestNcmMatrix *test, gconstpointer pdata)
{
  gchar *mser      = ncm_serialize_global_to_string (G_OBJECT (test->m), TRUE);
  NcmMatrix *m_dup = NCM_MATRIX (ncm_serialize_global_from_string (mser));
  guint i, j;

  g_free (mser);
  g_assert_cmpint (ncm_matrix_nrows (test->m), ==, ncm_matrix_nrows (m_dup));
  g_assert_cmpint (ncm_matrix_ncols (test->m), ==, ncm_matrix_ncols (m_dup));

  for (i = 0; i < ncm_matrix_nrows (test->m); i++)
  {
    for (j = 0; j < ncm_matrix_ncols (test->m); j++)
    {
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, i, j), ==, ncm_matrix_get (m_dup, i, j));
    }
  }

  NCM_TEST_FREE (ncm_matrix_free, m_dup);
}

void
test_ncm_matrix_free (TestNcmMatrix *test, gconstpointer pdata)
{
  NCM_TEST_FREE (ncm_matrix_free, test->m);

  if (test->d != NULL)
    g_free (test->d);
}

