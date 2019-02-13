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
void test_ncm_matrix_add_mul (TestNcmMatrix *test, gconstpointer pdata);
void test_ncm_matrix_log_exp (TestNcmMatrix *test, gconstpointer pdata);
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
  
  g_test_add ("/ncm/matrix/add_mul", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_add_mul,
              &test_ncm_matrix_free);
  
  g_test_add ("/ncm/matrix/log_exp", TestNcmMatrix, NULL,
              &test_ncm_matrix_new,
              &test_ncm_matrix_log_exp,
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

    g_assert (ncm_matrix_nrows (test->m) == gm->size1 && ncm_matrix_ncols (test->m) == gm->size2);
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

    g_assert (ga == ncm_matrix_get_array (test->m));
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
    
    g_assert ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));
  }
}

void
test_ncm_matrix_new_data_malloc (TestNcmMatrix *test, gconstpointer pdata)
{
  test->nrows = g_test_rand_int_range (1, 100);
  test->ncols = g_test_rand_int_range (1, 100);
  {
    gdouble *d = g_malloc ( test->nrows * test->ncols * sizeof (gdouble));
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

    g_assert ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));
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

    g_assert ((ncm_matrix_nrows (test->m) * ncm_matrix_ncols (test->m)) == (test->nrows * test->ncols));

    for (i = 0; i < 10 * test->nrows; i++)
    {
      const guint nr = g_test_rand_int_range (0, test->nrows);
      const guint nc = g_test_rand_int_range (0, test->ncols);
      const gdouble val = g_test_rand_double ();
      ncm_matrix_set (test->m, nr, nc, val);

      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, val);
    }

    {
      NcmVector *v = ncm_matrix_get_col (test->m, test->ncols - 1);
      g_assert (ncm_vector_len (v) == ncm_matrix_nrows (test->m));
      for (i = 0; i < test->nrows; i++)
      {
        ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_matrix_get (test->m, i, test->ncols - 1));
      }
      ncm_vector_free (v);
    }

    {
      NcmVector *v = ncm_matrix_get_row (test->m, test->nrows - 1);
      g_assert (ncm_vector_len (v) == ncm_matrix_ncols (test->m));
      for (i = 0; i < test->ncols; i++)
        g_assert (ncm_vector_get (v, i) == ncm_matrix_get (test->m, test->nrows - 1, i));
      ncm_vector_free (v);
    }
  }
}

void
test_ncm_matrix_sanity (TestNcmMatrix *test, gconstpointer pdata)
{
  guint i, j;

  g_assert (NCM_IS_MATRIX (test->m));

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
    const guint nr = g_test_rand_int_range (0, test->nrows);
    const guint nc = g_test_rand_int_range (0, test->ncols);
    const gdouble d = g_test_rand_double ();
    ncm_matrix_set (test->m, nr, nc, d);

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, d);
  }

  for (i = 0; i < 10 * test->nrows; i++)
  {
    guint nr, nc;
    gdouble *d;

    nr = g_test_rand_int_range (0, test->nrows);
    nc = g_test_rand_int_range (0, test->ncols);
    d = ncm_matrix_ptr (test->m, nr, nc);
    (*d) *= g_test_rand_double ();

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, *d);
  }

  {
    const gint nrows = (test->nrows <= test->ncols) ? 0 : test->nrows - test->ncols;
    const gint ncols = (test->ncols <= test->nrows) ? 0 : test->ncols - test->nrows;
    NcmMatrix *subm = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
    NcmMatrix *sm_cp = ncm_matrix_dup (subm);

    ncm_matrix_transpose (subm);
    
    for (i = 0; i < 10 * test->nrows; i++)
    {
      guint nr, nc;
      gdouble d;

      nr = g_test_rand_int_range (0, test->nrows - nrows);
      nc = g_test_rand_int_range (0, test->ncols - ncols);
      d = ncm_matrix_get (sm_cp, nr, nc);

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
    {
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, 1.0);
    }
    else
    {
      ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, 0.0);
    }
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
    d = ncm_matrix_get (test->m, nr, nc);
    ncm_matrix_scale (test->m, d1);

    ncm_assert_cmpdouble (ncm_matrix_get (test->m, nr, nc), ==, d * d1);
  }

  {
    NcmVector *v = ncm_vector_new (test->nrows);
    guint nc = g_test_rand_int_range (0, test->ncols);
    
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
test_ncm_matrix_submatrix (TestNcmMatrix *test, gconstpointer pdata)
{
  const gint nrows = g_test_rand_int_range (0, test->nrows);
  const gint ncols = g_test_rand_int_range (0, test->ncols);
  NcmMatrix *sm    = ncm_matrix_get_submatrix (test->m, nrows, ncols, test->nrows - nrows, test->ncols - ncols);
  guint ntests     = 20 * test->nrows;

  g_assert (ncm_matrix_nrows (sm) == (test->nrows - nrows) && ncm_matrix_ncols (sm) == (test->ncols - ncols));

  while (ntests--)
  {
    guint nr = g_test_rand_int_range (0, test->nrows - nrows);
    guint nc = g_test_rand_int_range (0, test->ncols - ncols);
    ncm_matrix_set (sm, nr, nc, g_test_rand_double ());
    ncm_assert_cmpdouble (ncm_matrix_get (sm, nr, nc), ==, ncm_matrix_get (test->m, nr + nrows, nc + ncols));
  }

  g_assert (NCM_IS_MATRIX (sm));

  ncm_matrix_free (test->m);
  g_assert (G_IS_OBJECT (test->m));
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

  g_assert (ncm_matrix_nrows (sm) == (test->nrows - nrows) && ncm_matrix_ncols (sm) == (test->ncols - ncols));

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
    NcmMatrix *osm = ncm_matrix_dup (sm);
    NcmMatrix *res = ncm_matrix_dup (sm);
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
  NcmMatrix *log_m = ncm_matrix_dup (sm);

  ncm_matrix_transpose (dup);

  ncm_matrix_add (sm, dup);
  ncm_matrix_scale (sm, 0.5);

  ncm_matrix_memcpy (dup, sm);

  ncm_matrix_sym_exp_cholesky (sm, 'U', exp_m);
  ncm_matrix_triang_to_sym (exp_m, 'U', TRUE);

  ncm_matrix_sym_posdef_log (exp_m, 'U', log_m);
  ncm_matrix_copy_triangle (log_m, 'U');

  ncm_assert_cmpdouble_e (ncm_matrix_cmp (dup, log_m, 0.0), ==, 0.0, 1.0e-7, 1.0e-7);
  
  NCM_TEST_FREE (ncm_matrix_free, sm);
  NCM_TEST_FREE (ncm_matrix_free, dup);
  NCM_TEST_FREE (ncm_matrix_free, exp_m);
  NCM_TEST_FREE (ncm_matrix_free, log_m);
}

void
test_ncm_matrix_serialization (TestNcmMatrix *test, gconstpointer pdata)
{
  gchar *mser = ncm_serialize_global_to_string (G_OBJECT (test->m), TRUE);
  NcmMatrix *m_dup = NCM_MATRIX (ncm_serialize_global_from_string (mser));
  gint i, j;
  
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
