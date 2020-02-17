/***************************************************************************
 *            test_ncm_vector.c
 *
 *  Tue April 03 16:02:26 2012
 *  Copyright  2012  Sandro Dias Pinto Vitenti
 *  <sandro@isoftware.com.br>
 ****************************************************************************/
/*
 * numcosmo
 * Copyright (C) Sandro Dias Pinto Vitenti 2012 <sandro@isoftware.com.br>
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

#define _TEST_NCM_VECTOR_STATIC_SIZE 100
#define _TEST_NCM_VECTOR_MIN_SIZE 5

typedef struct _TestNcmVector
{
  NcmVector *v;
  const NcmVector *cv;
  guint v_size;
  guint ntests;
} TestNcmVector;

void test_ncm_vector_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_gsl_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_gsl_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_array_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_array_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_slice_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_slice_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_malloc_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_malloc_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_static_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_static_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_const_new (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_const_free (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_sanity (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_operations (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_subvector (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_variant (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_serialization (TestNcmVector *test, gconstpointer pdata);
void test_ncm_vector_data_const_sanity (TestNcmVector *test, gconstpointer pdata);

gint
main (gint argc, gchar *argv[])
{
  g_test_init (&argc, &argv, NULL);
  ncm_cfg_init_full_ptr (&argc, &argv);
  ncm_cfg_enable_gsl_err_handler ();

  /* Default vector allocation */

  g_test_add ("/ncm/vector/default/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_free);

  g_test_add ("/ncm/vector/default/operations", TestNcmVector, NULL, 
              &test_ncm_vector_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_free);

  g_test_add ("/ncm/vector/default/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_free);

  g_test_add ("/ncm/vector/default/variant", TestNcmVector, NULL, 
              &test_ncm_vector_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_free);

  g_test_add ("/ncm/vector/default/serialization", TestNcmVector, NULL, 
              &test_ncm_vector_new, 
              &test_ncm_vector_serialization,
              &test_ncm_vector_free);

  /* GSL vector allocation */

  g_test_add ("/ncm/vector/gsl/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_gsl_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_gsl_free);
  
  g_test_add ("/ncm/vector/gsl/operations", TestNcmVector, NULL, 
              &test_ncm_vector_gsl_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_gsl_free);

  g_test_add ("/ncm/vector/gsl/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_gsl_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_gsl_free);

  g_test_add ("/ncm/vector/gsl/variant", TestNcmVector, NULL, 
              &test_ncm_vector_gsl_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_gsl_free);

  g_test_add ("/ncm/vector/gsl/serialization", TestNcmVector, NULL, 
              &test_ncm_vector_gsl_new, 
              &test_ncm_vector_serialization,
              &test_ncm_vector_gsl_free);

  /* Array vector allocation */

  g_test_add ("/ncm/vector/array/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_array_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_array_free);
  
  g_test_add ("/ncm/vector/array/operations", TestNcmVector, NULL, 
              &test_ncm_vector_array_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_array_free);

  g_test_add ("/ncm/vector/array/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_array_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_array_free);

  g_test_add ("/ncm/vector/array/variant", TestNcmVector, NULL, 
              &test_ncm_vector_array_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_array_free);
  
  g_test_add ("/ncm/vector/array/serialization", TestNcmVector, NULL, 
              &test_ncm_vector_array_new, 
              &test_ncm_vector_serialization,
              &test_ncm_vector_array_free);

  /* Data slice vector allocation */

  g_test_add ("/ncm/vector/data_slice/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_data_slice_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_data_slice_free);
  
  g_test_add ("/ncm/vector/data_slice/operations", TestNcmVector, NULL, 
              &test_ncm_vector_data_slice_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_data_slice_free);

  g_test_add ("/ncm/vector/data_slice/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_data_slice_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_data_slice_free);

  g_test_add ("/ncm/vector/data_slice/variant", TestNcmVector, NULL, 
              &test_ncm_vector_data_slice_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_data_slice_free);
  
  g_test_add ("/ncm/vector/data_slice/serialization", TestNcmVector, NULL, 
              &test_ncm_vector_data_slice_new, 
              &test_ncm_vector_serialization,
              &test_ncm_vector_data_slice_free);

  /* Data malloc vector allocation */

  g_test_add ("/ncm/vector/data_malloc/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_data_malloc_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_data_malloc_free);
  
  g_test_add ("/ncm/vector/data_malloc/operations", TestNcmVector, NULL, 
              &test_ncm_vector_data_malloc_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_data_malloc_free);

  g_test_add ("/ncm/vector/data_malloc/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_data_malloc_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_data_malloc_free);

  g_test_add ("/ncm/vector/data_malloc/variant", TestNcmVector, NULL, 
              &test_ncm_vector_data_malloc_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_data_malloc_free);

  /* Data static vector allocation */

  g_test_add ("/ncm/vector/data_static/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_data_static_new, 
              &test_ncm_vector_sanity, 
              &test_ncm_vector_data_static_free);

  g_test_add ("/ncm/vector/data_static/operations", TestNcmVector, NULL, 
              &test_ncm_vector_data_static_new, 
              &test_ncm_vector_operations, 
              &test_ncm_vector_data_static_free);

  g_test_add ("/ncm/vector/data_static/subvector", TestNcmVector, NULL, 
              &test_ncm_vector_data_static_new, 
              &test_ncm_vector_subvector,
              &test_ncm_vector_data_static_free);

  g_test_add ("/ncm/vector/data_static/variant", TestNcmVector, NULL, 
              &test_ncm_vector_data_static_new, 
              &test_ncm_vector_variant,
              &test_ncm_vector_data_static_free);

  g_test_add ("/ncm/vector/data_static/serialization", TestNcmVector, NULL, 
              &test_ncm_vector_data_static_new, 
              &test_ncm_vector_serialization,
              &test_ncm_vector_data_static_free);

  /* Data static vector allocation */

  g_test_add ("/ncm/vector/data_const/sanity", TestNcmVector, NULL, 
              &test_ncm_vector_data_const_new, 
              &test_ncm_vector_data_const_sanity, 
              &test_ncm_vector_data_const_free);

  g_test_run ();
}

void
_random_fill (NcmVector *v)
{
  gint i;
  for (i = 0; i < ncm_vector_len (v); i++)
    ncm_vector_set (v, i, g_test_rand_double ());
}

void
test_ncm_vector_new (TestNcmVector *test, gconstpointer pdata)
{
  test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  test->v = ncm_vector_new (test->v_size);
  _random_fill (test->v);
}

void
test_ncm_vector_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  NCM_TEST_FREE (ncm_vector_free, v);
}

void
test_ncm_vector_gsl_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  gsl_vector *gv = gsl_vector_alloc (v_size);
  NcmVector *v = test->v = ncm_vector_new_gsl (gv);
  ncm_assert_cmpdouble (ncm_vector_len (v), ==, gv->size);
  _random_fill (test->v);
}

void
test_ncm_vector_gsl_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  NCM_TEST_FREE (ncm_vector_free, v);
}

void
test_ncm_vector_array_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  GArray *ga = g_array_sized_new (FALSE, FALSE, sizeof (gdouble), v_size);
  NcmVector *v;
  
  g_array_set_size (ga, v_size);
  v = test->v = ncm_vector_new_array (ga);
  g_array_unref (ga);

  g_assert_cmpuint (ncm_vector_len (v), ==, ga->len);

  g_assert_true (ga == ncm_vector_get_array (v));
  g_array_unref (ga);
  _random_fill (test->v);
}

void
test_ncm_vector_array_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  
  NCM_TEST_FREE (ncm_vector_free, v);
}

void
test_ncm_vector_data_slice_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  gdouble *d = g_slice_alloc (v_size * sizeof (gdouble));
  NcmVector *v = test->v = ncm_vector_new_data_slice (d, v_size, 1);

  g_assert_cmpuint (ncm_vector_len (v), ==, v_size);
  _random_fill (test->v);
}

void
test_ncm_vector_data_slice_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;

  NCM_TEST_FREE (ncm_vector_free, v);  
}

void
test_ncm_vector_data_malloc_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = g_test_rand_int_range (_TEST_NCM_VECTOR_MIN_SIZE, _TEST_NCM_VECTOR_STATIC_SIZE);
  gdouble *d = g_malloc (v_size * sizeof (gdouble));
  NcmVector *v = test->v = ncm_vector_new_data_malloc (d, v_size, 1);

  g_assert_cmpuint (ncm_vector_len (v), ==, v_size);
  _random_fill (test->v);
}

void
test_ncm_vector_data_malloc_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;  
  NCM_TEST_FREE (ncm_vector_free, v);
}

void
test_ncm_vector_data_static_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = _TEST_NCM_VECTOR_STATIC_SIZE;
  static gdouble d[_TEST_NCM_VECTOR_STATIC_SIZE];
  NcmVector *v = test->v = ncm_vector_new_data_static (d, v_size, 1);

  g_assert_cmpuint (ncm_vector_len (v), ==, v_size);
  _random_fill (test->v);
}

void
test_ncm_vector_data_static_free (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  NCM_TEST_FREE (ncm_vector_free, v);
}

void
test_ncm_vector_data_const_new (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size = _TEST_NCM_VECTOR_STATIC_SIZE;
  const gdouble d[_TEST_NCM_VECTOR_STATIC_SIZE] = {0.0, };
  const NcmVector *v = test->cv = ncm_vector_const_new_data (d, v_size, 1);

  g_assert_cmpuint (ncm_vector_len (v), ==, v_size);
}

void
test_ncm_vector_data_const_free (TestNcmVector *test, gconstpointer pdata)
{
  const NcmVector *cv = test->cv;
  ncm_vector_const_free (cv);
}

void
test_ncm_vector_sanity (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size;
  NcmVector *v = test->v;
  guint i;

  g_assert_true (NCM_IS_VECTOR (v));

  for (i = 0; i < 10 * v_size; i++)
  {
    const guint n = g_test_rand_int_range (0, v_size - 1);
    const gdouble d = g_test_rand_double ();
    ncm_vector_set (v, n, d);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d);
  }
}

void
test_ncm_vector_data_const_sanity (TestNcmVector *test, gconstpointer pdata)
{
  g_assert_true (NCM_IS_VECTOR (test->cv));
}

void
test_ncm_vector_operations (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size;
  NcmVector *v = test->v;
  NcmVector *cv = ncm_vector_dup (v);
  guint i;

  for (i = 0; i < 10 * v_size; i++)
  {
    const guint n = g_test_rand_int_range (0, v_size - 1);
    const gdouble d = g_test_rand_double ();
    ncm_vector_set (v, n, d);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble *d;

    n = g_test_rand_int_range (0, v_size - 1);
    d = ncm_vector_ptr (v, n);
    (*d) *= g_test_rand_double ();

    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, *d);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);
    d = ncm_vector_get (v, n);
    d += d1;

    ncm_vector_addto (v, n, d1);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);
    d = ncm_vector_get (v, n);
    d -= d1;

    ncm_vector_subfrom (v, n, d1);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);
    ncm_vector_set_all (v, d1);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d1);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d, d1 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);
    d = ncm_vector_get (v, n);
    ncm_vector_scale (v, d1);
    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d * d1);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_div (v, cv);

    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d1 / d2);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_add (v, cv);

    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d1 + d2);
  }

  for (i = 0; i < 10 * v_size; i++)
  {
    guint n;
    gdouble d1 = g_test_rand_double (),
    d2 = g_test_rand_double ();

    n = g_test_rand_int_range (0, v_size - 1);

    ncm_vector_set_all (v, d1);
    ncm_vector_set_all (cv, d2);

    ncm_vector_sub (v, cv);

    ncm_assert_cmpdouble (ncm_vector_get (v, n), ==, d1 - d2);
  }

  ncm_vector_set_zero (v);
  for (i = 0; i < v_size; i++)
    ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, 0.0);

  ncm_vector_memcpy (v, cv);
  for (i = 0; i < v_size; i++)
    ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_vector_get (cv, i));

  ncm_vector_set_zero (v);
  ncm_vector_memcpy2 (v, cv, 2, 0, v_size - 2);

  ncm_assert_cmpdouble (ncm_vector_get (v, 0), ==, 0.0);
  ncm_assert_cmpdouble (ncm_vector_get (v, 1), ==, 0.0);
  for (i = 0; i < v_size - 2; i++)
    ncm_assert_cmpdouble (ncm_vector_get (v, i + 2), ==, ncm_vector_get (cv, i));

	ncm_vector_clear (&cv);
}

void
test_ncm_vector_subvector (TestNcmVector *test, gconstpointer pdata)
{
  guint v_size = test->v_size;
  NcmVector *v = test->v;
  NcmVector *sv = ncm_vector_get_subvector (v, 1, v_size - 1);
  guint ntests = 20 * v_size;

  g_assert_cmpuint (ncm_vector_len (sv), ==, (v_size - 1));

  while (ntests--)
  {
    guint i = g_test_rand_int_range (0, v_size - 1);
    ncm_vector_set (sv, i, g_test_rand_double ());
    ncm_assert_cmpdouble (ncm_vector_get (sv, i), ==, ncm_vector_get (v, i + 1));
  }

  g_assert_true (NCM_IS_VECTOR (sv));

  ncm_vector_free (v);
  g_assert_true (G_IS_OBJECT (v));
  ncm_vector_ref (v);

  NCM_TEST_FREE (ncm_vector_free, sv);
}

void
test_ncm_vector_variant (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  GVariant *var = ncm_vector_get_variant (v);
  g_assert_true (!g_variant_is_floating (var));
  g_assert_true (g_variant_is_container (var));
  g_assert_cmpuint (ncm_vector_len (v), ==, g_variant_n_children (var));

  {
    NcmVector *nv = ncm_vector_new_variant (var);
    gint i;

    g_assert_cmpuint (ncm_vector_len (v), ==, ncm_vector_len (nv));
    for (i = 0; i < ncm_vector_len (v); i++)
    {
      ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_vector_get (nv, i));
    }
    NCM_TEST_FREE (ncm_vector_free, nv);
  }
  
  g_variant_unref (var);
  NCM_TEST_FAIL (g_variant_unref (var); fprintf (stderr, "fail (%s)", g_variant_get_type_string (var)));
}

void
test_ncm_vector_serialization (TestNcmVector *test, gconstpointer pdata)
{
  NcmVector *v = test->v;
  gchar *vser = ncm_serialize_global_to_string (G_OBJECT (v), TRUE);
  NcmVector *v_dup = NCM_VECTOR (ncm_serialize_global_from_string (vser));
  gint i;
  g_free (vser);
  g_assert_cmpint (ncm_vector_len (v), ==, ncm_vector_len (v_dup));

  for (i = 0; i < ncm_vector_len (v); i++)
  {
    ncm_assert_cmpdouble (ncm_vector_get (v, i), ==, ncm_vector_get (v_dup, i));
  }

  NCM_TEST_FREE (ncm_vector_free, v_dup);
}
